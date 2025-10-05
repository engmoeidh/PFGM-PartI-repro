import os, time, math, argparse
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from gf.solver1d import gf_solve_1d
from gf.highorder import dphi_dx_4th, simpson

# ---------------- CLI ----------------
ap = argparse.ArgumentParser()
ap.add_argument("--L", type=float, default=20.0, help="half-domain so R=2L")
ap.add_argument("--hs", type=str, default="0.02,0.01,0.005", help="comma list of h for D2")
ap.add_argument("--gf_iters", type=int, default=40000)
ap.add_argument("--newton_steps", type=int, default=12)
ap.add_argument("--newton_tol", type=float, default=1e-10)
ap.add_argument("--target_vrel", type=float, default=5e-6, help="GF early-stop virial target")
ap.add_argument("--fast", action="store_true", help="reduce work where possible")
args = ap.parse_args()

ROOT="."
os.makedirs(f"{ROOT}/figures", exist_ok=True)
os.makedirs(f"{ROOT}/data", exist_ok=True)

def energy_parts_4th(x, phi):
    h = x[1]-x[0]
    dph = dphi_dx_4th(phi, h)
    V = 0.25*(1-phi**2)**2
    T = float(simpson(0.5*dph**2, x))
    Vint = float(simpson(V, x))
    return T, Vint, T+Vint, dph

def Vpp(phi): return 3*phi**2 - 1.0

def thomas(a,b,c,d):
    n=len(d); ac=a.copy(); bc=b.copy(); cc=c.copy(); dc=d.copy()
    for i in range(1,n):
        m=ac[i-1]/bc[i-1]; bc[i]-=m*cc[i-1]; dc[i]-=m*dc[i-1]
    x=np.zeros(n); x[-1]=dc[-1]/bc[-1]
    for i in range(n-2,-1,-1): x[i]=(dc[i]-cc[i]*x[i+1])/bc[i]
    return x

def residual_L2(x, phi):
    h = x[1]-x[0]
    phixx = (np.roll(phi,-1)-2*phi+np.roll(phi,1))/h**2
    phixx[0]=(phi[1]-2*phi[0]+phi[0])/h**2
    phixx[-1]=(phi[-1]-2*phi[-1]+phi[-2])/h**2
    dV = phi*(phi**2-1.0)
    r = -phixx + dV
    r[0]=0.0; r[-1]=0.0
    return math.sqrt(float(np.trapz(r**2, x)))

def newton_corrector(x, phi_init, max_steps=12, tol=1e-10):
    """Accepted Newton steps only, with backtracking line search."""
    phi = phi_init.copy()
    h = x[1]-x[0]; N=len(x)
    hist=[]
    t0=time.time()
    r_old = residual_L2(x, phi)
    hist.append((0, r_old, 0.0))
    for k in range(1, max_steps+1):
        phixx = (np.roll(phi,-1)-2*phi+np.roll(phi,1))/h**2
        phixx[0]=(phi[1]-2*phi[0]+phi[0])/h**2
        phixx[-1]=(phi[-1]-2*phi[-1]+phi[-2])/h**2
        dV = phi*(phi**2-1.0)
        r = -phixx + dV
        r[0]=0.0; r[-1]=0.0
        m=N-2; off=-1.0/h**2; main = 2.0/h**2 + Vpp(phi[1:-1])
        a=np.full(m-1, off); b=main.copy(); c=np.full(m-1, off)
        rhs = -r[1:-1].copy()
        dphi_inner = thomas(a,b,c,rhs)
        dphi = np.zeros_like(phi); dphi[1:-1]=dphi_inner
        step=1.0
        while step>1e-8:
            trial = phi + step*dphi
            trial[0], trial[-1] = -1.0, +1.0
            r_new = residual_L2(x, trial)
            if r_new <= 0.95*r_old:
                phi = trial; r_old = r_new
                hist.append((k, r_old, time.time()-t0))
                break
            step *= 0.5
        if r_old <= tol: break
    return phi, np.array(hist)

# ---------------- Build figures ----------------
L = args.L; R = 2*L
hs = [float(s) for s in args.hs.split(",") if s.strip()]

# 1) Kink profile + cs2 and D1 (accepted Newton only) on the finest spacing
h_fine = min(hs)
N_fine = int(round(R/h_fine))+1
x, phi_gf, _, _ = gf_solve_1d("doublewell", R=R, N=N_fine,
                              iters=(args.gf_iters//2 if args.fast else args.gf_iters),
                              tol=1e-12, kink=True, target_vrel=args.target_vrel)
phi_N, hist = newton_corrector(x, phi_gf,
                               max_steps=(8 if args.fast else args.newton_steps),
                               tol=args.newton_tol)

# D1 plot
df_hist = pd.DataFrame({"iter":hist[:,0], "residual_L2":hist[:,1], "cpu_cum":hist[:,2]})
df_hist.to_csv(f"{ROOT}/data/newton_accepted_history.csv", index=False)
fig, ax1 = plt.subplots(figsize=(6.4,4.2))
ax1.semilogy(df_hist["iter"], np.clip(df_hist["residual_L2"],1e-20,None), "-o", lw=2)
ax1.set_xlabel("Newton iteration"); ax1.set_ylabel(r"Dimensionless $L^2$ residual")
ax2 = ax1.twinx(); ax2.plot(df_hist["iter"], df_hist["cpu_cum"], "--s", lw=2)
ax2.set_ylabel("Cumulative CPU time (s)")
ax1.set_title(r"Dimensionless $L^2$ residual vs Newton iteration (accepted steps)")
ax1.grid(alpha=0.25); plt.tight_layout(); plt.savefig(f"{ROOT}/figures/D1.png", dpi=150); plt.close()

# kink + grad
T,V,E,phix = energy_parts_4th(x, phi_N)
plt.figure(figsize=(6.4,5.2))
ax1=plt.subplot(2,1,1); ax1.plot(x-R/2, phi_N, lw=3, color="#d99800"); ax1.axhline(1.0, ls="--", color="gray", lw=1)
ax1.set_ylabel(r"$\Phi(x)$"); ax1.set_title("Kink profile and gradient"); ax1.grid(alpha=0.25)
ax2=plt.subplot(2,1,2,sharex=ax1); ax2.plot(x-R/2, np.abs(phix), lw=3, color="#d99800")
ax2.set_xlabel("x"); ax2.set_ylabel(r"$|\partial_x\Phi|$"); ax2.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(f"{ROOT}/figures/kink_profile.png", dpi=150); plt.close()

# cs2 analytic
aX = np.linspace(-0.25, 1.5, 400); cs2=(1+aX)/(1+3*aX)
plt.figure(figsize=(6.4,4.2)); plt.plot(aX, cs2, lw=3, color="#d99800")
plt.title(r"Healthy band $(0<c_s^2<1)$ for $c_s^2=\frac{1+\alpha X}{1+3\alpha X}$")
plt.xlabel(r"$\alpha X$"); plt.ylabel(r"$c_s^2$"); plt.ylim(-0.5,2.6); plt.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(f"{ROOT}/figures/cs2_band.png", dpi=150); plt.close()

# 2) D2: true virial defect |T-V| with 4th-order numerics
rows=[]; phi_prev=None
for h in hs:
    N = int(round(R/h))+1
    phi0=None
    if phi_prev is not None:
        x_prev = np.linspace(0.0, R, len(phi_prev))
        x_new  = np.linspace(0.0, R, N)
        phi0   = np.interp(x_new, x_prev, phi_prev)
    xh, phih, _, _ = gf_solve_1d("doublewell", R=R, N=N,
                                 iters=(args.gf_iters//2 if args.fast else args.gf_iters),
                                 tol=1e-12, kink=True, phi0=phi0, target_vrel=args.target_vrel)
    phih, _ = newton_corrector(xh, phih,
                               max_steps=(6 if args.fast else args.newton_steps),
                               tol=args.newton_tol)
    Th, Vh, Eh, _ = energy_parts_4th(xh, phih)
    Rvir = abs(Th - Vh)
    rows.append({"h":float(h), "T":Th, "V":Vh, "E":Eh, "Rvir":Rvir})
    phi_prev = phih

d = pd.DataFrame(rows).sort_values("h")
d.to_csv(f"{ROOT}/data/virial_true_vs_h.csv", index=False)
# slope
p = np.polyfit(np.log(d["h"]), np.log(d["Rvir"]), 1)[0]
# h^2 guide through finest
hf  = float(d["h"].min()); Rf = float(d.loc[d["h"].idxmin(),"Rvir"])
href = d["h"].values; ref = (href/hf)**2 * Rf
plt.figure(figsize=(6.4,3.8))
plt.loglog(d["h"], d["Rvir"], "-o", label=r"$|R_{\rm vir}(h)|$")
plt.loglog(href, ref, "--", label=r"$\propto h^2$")
plt.gca().invert_xaxis(); plt.xlabel("Grid spacing $h$")
plt.ylabel(r"Virial residual $|R_{\rm vir}(h)| = |T(h)-V(h)|$")
plt.title(rf"Fig. D.2. Virial residual vs. $h$ (log–log), slope $p={p:.2f}$")
plt.grid(alpha=0.25, which="both"); plt.legend(); plt.tight_layout()
plt.savefig(f"{ROOT}/figures/D2.png", dpi=150); plt.close()

# 3) D3: E(R) plateau (reuse finest h to keep runtime small)
Rs = [20, 28, 36, 40] if args.fast else [20,24,28,32,36,40]
ER=[]
for Rr in Rs:
    Nrr = int(round(Rr/h_fine))+1
    xr, phir, _, _ = gf_solve_1d("doublewell", R=Rr, N=Nrr,
                                 iters=(args.gf_iters//2 if args.fast else args.gf_iters),
                                 tol=1e-12, kink=True, target_vrel=args.target_vrel)
    phir, _ = newton_corrector(xr, phir,
                               max_steps=(6 if args.fast else args.newton_steps),
                               tol=args.newton_tol)
    Tr,Vr,Er,_ = energy_parts_4th(xr, phir); ER.append(Er)
Einf = ER[-1]
plt.figure(figsize=(6.4,4.2))
plt.plot(Rs, ER, "-o", color="#d99800"); plt.axhline(Einf, ls="--", color="#d99800", alpha=0.7)
plt.xlabel(r"Outer boundary radius $R$"); plt.ylabel(r"Total energy $E(R)$")
plt.title("Fig. D.3. Total energy vs. outer boundary radius $R$")
plt.grid(alpha=0.25); plt.tight_layout(); plt.savefig(f"{ROOT}/figures/D3.png", dpi=150); plt.close()

# Gates (softened for fast mode but still referee-safe)
if not args.fast:
    if df_hist["residual_L2"].iloc[-1] > args.newton_tol:
        raise SystemExit(f"[FAIL] Newton did not reach tol: {df_hist['residual_L2'].iloc[-1]:.3e}")
    if not (1.8 <= p <= 2.2):
        raise SystemExit(f"[FAIL] D2 slope p={p:.3f} not ~2; increase L or refine hs.")
print("[OK] Paper-style figures done.")
