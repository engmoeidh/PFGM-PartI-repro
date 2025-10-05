import os, time, math, numpy as np, pandas as pd
import matplotlib.pyplot as plt

ROOT="."
g = pd.read_csv(os.path.join(ROOT,"data","grid_refinement.csv")).sort_values("h")
finest_h = float(g["h"].min())
R = float(g["R"].iloc[0]) if "R" in g.columns else 20.0
N = int(round(R/finest_h))+1

# --- Recover a finest-grid profile by rerunning one short relax (uses your semi-implicit GF)
from gf.solver1d import gf_solve_1d

x, phi, dfh, _ = gf_solve_1d("doublewell", R=R, N=N, iters=60000, tol=1e-12, kink=True, target_vrel=5e-5)
h = x[1]-x[0]

def energy_parts(phi):
    phix = (np.roll(phi,-1)-np.roll(phi,1))/(2*h)
    phix[0]=(phi[1]-phi[0])/h; phix[-1]=(phi[-1]-phi[-2])/h
    V = 0.25*(1-phi**2)**2
    T = float(np.trapz(0.5*phix**2, x))
    Vint=float(np.trapz(V, x))
    return T, Vint, T+Vint, phix

# =============== Fig: kink profile + gradient ===============
T,V,E,phix = energy_parts(phi)
plt.figure(figsize=(6.4,5.2))
ax1 = plt.subplot(2,1,1)
ax1.plot(x, phi, lw=3, color="#d99800")
ax1.axhline(1.0, ls="--", color="gray", lw=1)
ax1.set_ylabel(r"$\Phi(x)$")
ax1.set_title("Kink profile and gradient")
ax1.grid(alpha=0.25)

ax2 = plt.subplot(2,1,2, sharex=ax1)
ax2.plot(x, np.abs(phix), lw=3, color="#d99800")
ax2.set_xlabel("x"); ax2.set_ylabel(r"$|\partial_x \Phi|$")
ax2.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(ROOT,"figures","kink_profile.png"), dpi=150)
plt.close()

# =============== Fig: c_s^2(αX) band ===============
# X := (∂x Φ)^2 in 1D; we plot analytic curve vs αX
aX = np.linspace(-0.25, 1.5, 400)
cs2 = (1.0 + aX)/(1.0 + 3.0*aX)
plt.figure(figsize=(6.4,4.2))
plt.plot(aX, cs2, lw=3, color="#d99800")
plt.title(r"Healthy band $(0<c_s^2<1)$ for $c_s^2=\frac{1+\alpha X}{1+3\alpha X}$")
plt.xlabel(r"$\alpha X$"); plt.ylabel(r"$c_s^2$")
plt.ylim(-0.5, 2.6)
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(ROOT,"figures","cs2_band.png"), dpi=150)
plt.close()

# =============== Fig: Newton residual vs iter + CPU time ===============
# Small Newton corrector on the EL=0 operator L = -d^2/dx^2 + V''(phi)
def Vpp(phi): return 3*phi**2 - 1.0
def thomas(a,b,c,d):
    n=len(d); ac=a.copy(); bc=b.copy(); cc=c.copy(); dc=d.copy()
    for i in range(1,n):
        m = ac[i-1]/bc[i-1]
        bc[i]-=m*cc[i-1]; dc[i]-=m*dc[i-1]
    x=np.zeros(n); x[-1]=dc[-1]/bc[-1]
    for i in range(n-2,-1,-1): x[i]=(dc[i]-cc[i]*x[i+1])/bc[i]
    return x

def newton_refine(phi, steps=12):
    res_hist=[]; t_hist=[]; t0=time.time()
    for k in range(steps+1):
        # residual r = -phi'' + dV/dphi
        phixx = (np.roll(phi,-1)-2*phi+np.roll(phi,1))/h**2
        phixx[0]=(phi[1]-2*phi[0]+phi[0])/h**2
        phixx[-1]=(phi[-1]-2*phi[-1]+phi[-2])/h**2
        dV = phi*(phi**2-1.0)
        r = -phixx + dV
        r[0]=0.0; r[-1]=0.0
        res = math.sqrt(float(np.trapz(r**2, x)))
        res_hist.append(max(res,1e-16)); t_hist.append(time.time()-t0)
        if k==steps: break
        # build tridiagonal for L δφ = -r, with L = -d2 + V''
        m=len(phi)-2; off = -1.0/h**2; main = 2.0/h**2 + Vpp(phi[1:-1])
        a=np.full(m-1, off); b=main.copy(); c=np.full(m-1, off)
        d = -r[1:-1].copy()
        d[0]  -= off*(-1.0 - phi[1])   # boundary shift: enforce δφ(0)=0
        d[-1] -= off*(+1.0 - phi[-2])  # δφ(R)=0
        dphi_inner = thomas(a,b,c,d)
        dphi = np.zeros_like(phi); dphi[1:-1]=dphi_inner
        phi += dphi
        phi[0]=-1.0; phi[-1]=+1.0
    return np.array(res_hist), np.array(t_hist)

resN, tN = newton_refine(phi.copy(), steps=12)
fig, ax1 = plt.subplots(figsize=(6.4,4.2))
ax1.semilogy(range(len(resN)), resN, marker="o", lw=2)
ax1.set_xlabel("Newton iteration"); ax1.set_ylabel(r"Dimensionless $L^2$ residual")
ax1.set_title("Dimensionless $L^2$-norm of the nonlinear residual vs. Newton iteration\nwith cumulative CPU time")
ax2 = ax1.twinx()
ax2.plot(range(len(tN)), tN, "--s", lw=2)
ax2.set_ylabel("Cumulative CPU time (s)")
ax1.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(ROOT,"figures","D1.png"), dpi=150)
plt.close()

# =============== Fig: Virial residual vs h (log–log) with ∝h^2 guide ===============
# Use absolute virial mismatch from your grid_refinement.csv
hvals = g["h"].values; Rvir = g["virial_mismatch_abs"].values
plt.figure(figsize=(6.4,3.8))
plt.loglog(hvals, Rvir, "-o")
# reference slope 2 through the finest point
hf, Rf = hvals.min(), Rvir[g["h"].idxmin()]
ref = (hvals/hf)**2 * Rf
plt.loglog(hvals, ref, "--", alpha=0.7)
plt.gca().invert_xaxis()
plt.xlabel("Grid spacing $h$"); plt.ylabel(r"Virial residual $|R_{\rm vir}(h)|$")
plt.title("Fig. D.2. Virial residual vs. grid spacing $h$ (log–log)")
plt.grid(alpha=0.25, which="both")
plt.tight_layout()
plt.savefig(os.path.join(ROOT,"figures","D2.png"), dpi=150)
plt.close()

# =============== Fig: Total energy vs outer radius R ===============
Rs = [20,24,28,32,36,40]
E_R=[]
for Rr in Rs:
    Nrr = int(round(Rr/finest_h))+1
    xr, phir, dfr, _ = gf_solve_1d("doublewell", R=Rr, N=Nrr, iters=40000, tol=1e-12, kink=True, target_vrel=5e-5)
    hr = xr[1]-xr[0]
    phix = (np.roll(phir,-1)-np.roll(phir,1))/(2*hr)
    phix[0]=(phir[1]-phir[0])/hr; phix[-1]=(phir[-1]-phir[-2])/hr
    V = 0.25*(1-phir**2)**2
    E_R.append(float(np.trapz(0.5*phix**2 + V, xr)))
Einf = E_R[-1] + (E_R[-1]-E_R[-2])  # simple visual guide
plt.figure(figsize=(6.4,4.2))
plt.plot(Rs, E_R, "-o", color="#d99800")
plt.axhline(Einf, ls="--", color="#d99800", alpha=0.7)
plt.xlabel(r"Outer boundary radius $R$"); plt.ylabel(r"Total energy $E(R)$")
plt.title("Fig. D.3. Total energy vs. outer boundary radius $R$")
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(ROOT,"figures","D3.png"), dpi=150)
plt.close()

# =============== Fig: Energy vs dilation E(λ) ===============
# Build φ_λ(x) = φ(λ x). For each λ, resample and integrate on the same x-window.
lams = np.linspace(0.5, 2.5, 9)
E_lam=[]
for lam in lams:
    xs = x/lam     # so that φ_λ(x) = φ(xs)
    phi_lam = np.interp(xs, x, phi, left=-1.0, right=+1.0)
    phix_l = (np.roll(phi_lam,-1)-np.roll(phi_lam,1))/(2*h)
    phix_l[0]=(phi_lam[1]-phi_lam[0])/h; phix_l[-1]=(phi_lam[-1]-phi_lam[-2])/h
    V = 0.25*(1-phi_lam**2)**2
    E_lam.append(float(np.trapz(0.5*phix_l**2 + V, x)))
plt.figure(figsize=(6.4,4.0))
plt.plot(lams, E_lam, "-", color="#d99800", lw=3)
plt.xlabel(r"$\lambda$"); plt.ylabel(r"$E(\lambda)$"); plt.title(r"Energy vs dilation $E(\lambda)$")
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(ROOT,"figures","E_lambda_curve.png"), dpi=150)
plt.close()

print("[OK] original-style figures refreshed.")
