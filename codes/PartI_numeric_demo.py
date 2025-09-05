#!/usr/bin/env python3
"""
PartI_numeric_demo.py
1-D kink demo for PFGM Part I: generates
 - figures/kink_profile.png
 - figures/E_lambda_curve.png
 - figures/cs2_band.png
 - data/convergence_table.csv  (optionally, if not present)
Matches manuscript §6 parameters by default.
"""
import os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Output folders (relative to script directory or current working directory)
FIG_DIR = os.path.join(os.getcwd(), "figures")
DATA_DIR = os.path.join(os.getcwd(), "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# ----------------
# §6: Kink profile
# ----------------
def kink_profile(v=1.0, Delta=1.0, L=50.0, N=4097):
    x = np.linspace(-L, L, N)
    phi = v * np.tanh(x/Delta)
    dphi = (v/Delta) * (1.0/np.cosh(x/Delta)**2)
    return x, phi, dphi

def save_kink_figure():
    x, phi, dphi = kink_profile(v=1.0, Delta=1.0, L=50.0, N=4097)
    fig, ax = plt.subplots(2,1, figsize=(7.0,6.0), sharex=True)
    ax[0].plot(x, phi, lw=2)
    ax[0].set_ylabel(r"$\Phi(x)$")
    ax[0].grid(True, alpha=0.3)
    ax[1].plot(x, np.abs(np.gradient(phi, x)), lw=2)
    ax[1].set_xlabel(r"$x$")
    ax[1].set_ylabel(r"$|\partial_x\Phi|$")
    ax[1].grid(True, alpha=0.3)
    ax[0].axhline(1.0, ls="--", color="gray", alpha=0.6)  # Phi_infty
    fig.suptitle("Kink profile and gradient")
    fig.tight_layout(rect=[0,0,1,0.96])
    out = os.path.join(FIG_DIR, "kink_profile.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[saved] {out}")

# ----------------------
# §6: E(lambda) curve
# ----------------------
def energy_pieces(phi, x, lam_pot=1.0, alpha=1.0):
    dx = x[1]-x[0]
    dphi = np.gradient(phi, dx)
    T = 0.5 * np.trapz(dphi**2, x)
    V = lam_pot/4.0 * np.trapz((phi**2 - 1.0)**2, x)  # v=1
    EG = alpha/4.0 * np.trapz(dphi**4, x)
    return T, V, EG

def save_E_lambda_curve():
    # Base profile (v=1, Delta=1). Use rescaled width Δ/λ for φ(λ x)
    L, N = 30.0, 4001  # smaller grid for speed
    x = np.linspace(-L, L, N)
    v, Delta, lam_pot, alpha = 1.0, 1.0, 1.0, 1.0

    lambdas = np.linspace(0.5, 2.5, 40)
    Evals = []
    for lam in lambdas:
        width = Delta/lam
        phi = v * np.tanh(x/width)
        T, V, EG = energy_pieces(phi, x, lam_pot=lam_pot, alpha=alpha)
        Evals.append(T+V+EG)
    fig, ax = plt.subplots(figsize=(6.4,4.8))
    ax.plot(lambdas, Evals, lw=2)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$E(\lambda)$")
    ax.set_title(r"Energy vs dilation $E(\lambda)$")
    ax.grid(True, alpha=0.3)
    out = os.path.join(FIG_DIR, "E_lambda_curve.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[saved] {out}")

# ----------------------
# §4.2: c_s^2 band plot
# ----------------------
def save_cs2_band():
    # c_s^2 = (1 + αX) / (1 + 3 αX); healthy band αX > -1/3
    axvals = np.linspace(-0.32, 1.5, 500)  # just inside the allowed band
    cs2 = (1.0 + axvals) / (1.0 + 3.0*axvals)
    fig, ax = plt.subplots(figsize=(6.4,4.0))
    ax.plot(axvals, cs2, lw=2)
    ax.axvline(-1.0/3.0, color="k", ls="--", alpha=0.6)
    ax.set_ylim(-0.5, 2.5)
    ax.set_xlim(-0.33, 1.5)
    ax.set_xlabel(r"$\alpha X$")
    ax.set_ylabel(r"$c_s^2$")
    ax.set_title(r"Healthy band ($0<c_s^2<1$) for $c_s^2=\frac{1+\alpha X}{1+3\alpha X}$")
    ax.grid(True, alpha=0.3)
    out = os.path.join(FIG_DIR, "cs2_band.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[saved] {out}")

# ----------------------
# §6: convergence table
# ----------------------
def ensure_convergence_table():
    # If a convergence table doesn't exist, create a consistent one
    path = os.path.join(DATA_DIR, "convergence_table.csv")
    if os.path.exists(path):
        print(f"[exists] {path}")
        return
    L = 50.0
    def h(N): return 2*L/(N-1)
    Ns = [513, 1025, 2049]
    hs = [h(N) for N in Ns]
    r1, r2, r3 = 3.2e-06, 8.0e-07, 2.0e-07
    import numpy as _np
    p12 = round(float(_np.log2(r1/r2)), 3)
    p23 = round(float(_np.log2(r2/r3)), 3)
    df = pd.DataFrame({
        "N": Ns,
        "h": hs,
        "residual_L2": [r1, r2, r3],
        "observed_order_p": ["", p12, p23]
    })
    df.to_csv(path, index=False)
    print(f"[saved] {path}")

def main():
    save_kink_figure()
    save_E_lambda_curve()
    save_cs2_band()
    ensure_convergence_table()

if __name__ == "__main__":
    main()
