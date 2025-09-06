#!/usr/bin/env python3
"""
PartI_grid_solver.py
Radial finite-difference solver for the stationary PFGM equation (Appendix F),
with diagnostics to regenerate Figs D1–D3 and a grid refinement CSV.

Equation (dimensionless form, App. F):
  (1/r^2) d/dr [ r^2 φ' + α r^2 |φ'|^2 φ' ] = V'(φ)
We choose V(φ) = 1/2 φ^2  (so V'(φ)=φ), α>0.

Boundary conditions:
  φ'(0)=0 (regularity), φ(R)=0 (vacuum at outer boundary).

Outputs:
 - figures/D1.png (Newton residual vs iteration)
 - figures/D2.png (virial mismatch vs h)
 - figures/D3.png (total energy vs outer radius R)
 - data/grid_refinement.csv
 - data/virial_diagnostics.csv
"""
import os, math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags

FIG_DIR = os.path.join(os.getcwd(), "figures")
DATA_DIR = os.path.join(os.getcwd(), "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

def make_grid(R=20.0, N=1601):
    r = np.linspace(0.0, R, N)
    h = r[1]-r[0]
    return r, h

def residual(phi, r, h, alpha=0.3):
    # Compute residual F(phi) = (1/r^2) d/dr [ r^2 φ' + α r^2 |φ'|^2 φ' ] - φ
    N = len(r)
    # centered differences for φ'
    phi_p = np.zeros_like(phi)
    phi_p[1:-1] = (phi[2:] - phi[:-2])/(2*h)
    phi_p[0] = 0.0  # regularity φ'(0)=0
    phi_p[-1] = (phi[-1]-phi[-2])/h  # one-sided; will be adjusted by BC φ(R)=0

    flux = phi_p + alpha*(np.abs(phi_p)**2)*phi_p
    # radial divergence in conservative form
    # (1/r^2) d/dr ( r^2 * flux )
    r2 = r**2
    # midpoint r_{i+1/2}, r_{i-1/2}
    rph = (r[1:]+r[:-1])/2.0
    rph2 = rph**2

    # flux at midpoints
    flux_m = 0.5*(flux[1:]+flux[:-1])

    div = np.zeros_like(phi)
    # interior points: i=1..N-2
    div[1:-1] = ( rph2[1:]*flux_m[1:] - rph2[:-1]*flux_m[:-1] )/(h * r2[1:-1])
    # origin: use l'Hospital/regularity -> approximate with second point
    div[0] = div[1]
    # outer boundary: enforce Dirichlet φ(R)=0 via penalty in residual
    div[-1] = div[-2]
    # residual with V'(φ)=φ
    F = div - phi
    # enforce φ(R)=0 strongly in residual
    F[-1] = phi[-1]
    # enforce φ'(0)=0 weakly: already done via derivative definition
    return F


def jacobian(phi, r, h, alpha=0.3):
    """
    Build a strictly tridiagonal Jacobian with sizes:
      - main diagonal (offset 0): length N
      - lower diagonal (offset -1): length N-1  (rows 1..N-1, cols 0..N-2)
      - upper diagonal (offset +1): length N-1  (rows 0..N-2, cols 1..N-1)
    Enforce φ(R)=0 via a Dirichlet row at i=N-1.
    """
    N = len(r)

    # φ' (centered in the interior, one sided at the outer boundary)
    phi_p = np.zeros_like(phi)
    phi_p[1:-1] = (phi[2:] - phi[:-2])/(2*h)
    phi_p[0] = 0.0
    phi_p[-1] = (phi[-1] - phi[-2]) / h

    # d(flux)/d(φ'), where flux = φ' + α|φ'|^2 φ'
    dflux_dphip = 1.0 + 3.0*alpha*(phi_p**2)

    r2   = r**2
    rph  = 0.5*(r[1:] + r[:-1])
    rph2 = rph**2

    # allocate full-length bands
    a = np.zeros(N)  # lower band (we'll trim to length N-1 later)
    b = np.zeros(N)  # main band
    c = np.zeros(N)  # upper band (we'll trim to length N-1 later)

    # interior points: i = 1..N-2
    for i in range(1, N-1):
        # effective 'conductivities' at i±1/2
        k_p = dflux_dphip[i+1]     # at i+1/2
        k_m = dflux_dphip[i-1]     # at i-1/2

        Ap = rph2[i]   * k_p / (h * r2[i])
        Am = rph2[i-1] * k_m / (h * r2[i])

        # divergence stencil → (Am+Ap)/h on main, -Am/h on lower, -Ap/h on upper
        a[i] = -Am / h                  # contributes to offset -1
        c[i] = -Ap / h                  # contributes to offset +1
        b[i] = (Am + Ap) / h - 1.0      # -dV'/dφ with V'(φ)=φ → 1

    # origin row (i=0): natural reflecting BC φ'(0)=0 (simple regularization)
    b[0] = b[1]
    c[0] = c[1]
    a[0] = 0.0

    # Dirichlet row at i=N-1: φ(R)=0
    a[-1] = 0.0
    b[-1] = 1.0
    c[-1] = 0.0

    # *** Critical: trim off-diagonals to length N-1 ***
    a_trim = a[1:]      # rows 1..N-1, cols 0..N-2 → offset -1
    c_trim = c[:-1]     # rows 0..N-2, cols 1..N-1 → offset +1

    # Build sparse tridiagonal with explicit shape
    J = diags(
        diagonals=[a_trim, b, c_trim],
        offsets=[-1, 0, 1],
        shape=(N, N),
        format="csc"
    )
    return J


def energy_components(phi, r, h, alpha=0.3):
    # T = 1/2 ∫ |∇φ|^2 d^3x = 2π ∫ r^2 (φ')^2 dr * 2 (angular) => 4π ∫ r^2 (φ')^2 dr
    # For normalization, we report dimensionless T,V,E_G ignoring 4π factor consistently.
    N = len(r)
    phi_p = np.zeros_like(phi)
    phi_p[1:] = np.diff(phi)/h
    r2 = r**2
    T = 0.5*np.trapz(r2*(phi_p**2), r)          # up to constant factor
    V = 0.5*np.trapz(r2*(phi**2), r)            # V = 1/2 φ^2
    EG = 0.25*alpha*np.trapz(r2*(phi_p**4), r)  # quartic gradient
    return T, V, EG

def newton_solve(R=20.0, N=1601, alpha=0.3, tol=1e-10, maxit=30):
    r, h = make_grid(R, N)
    # initial guess: localized bump
    phi = np.exp(-(r/2.0)**2)
    res_hist = []

    for k in range(maxit):
        F = residual(phi, r, h, alpha=alpha)
        nrm = float(np.linalg.norm(F, 2))
        res_hist.append(nrm)
        if nrm < tol:
            break
        J = jacobian(phi, r, h, alpha=alpha)
        delta = spsolve(J, -F)
        # line search
        lam = 1.0
        for _ in range(10):
            phit = phi + lam*delta
            if np.linalg.norm(residual(phit, r, h, alpha=alpha), 2) < nrm:
                phi = phit
                break
            lam *= 0.5
    return r, h, phi, res_hist

def fig_D1(res_hist):
    fig, ax = plt.subplots(figsize=(6.4,4.8))
    ax.semilogy(res_hist, lw=2)
    ax.set_xlabel("Newton iteration")
    ax.set_ylabel(r"$\|F\|_2$")
    ax.set_title("Newton–Raphson convergence")
    ax.grid(True, which="both", alpha=0.3)
    out = os.path.join(FIG_DIR, "D1.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[saved] {out}")

def fig_D2(grid_table):
    fig, ax = plt.subplots(figsize=(6.4,4.8))
    ax.loglog(grid_table["h"], grid_table["virial_mismatch_abs"], 'o-', lw=2)
    ax.set_xlabel("grid spacing h")
    ax.set_ylabel(r"$|E_G - (T+3V)|$")
    ax.set_title("Virial mismatch vs grid spacing")
    ax.grid(True, which="both", alpha=0.3)
    out = os.path.join(FIG_DIR, "D2.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[saved] {out}")

def fig_D3(R_list, E_list):
    fig, ax = plt.subplots(figsize=(6.4,4.8))
    ax.plot(R_list, E_list, 'o-', lw=2)
    ax.set_xlabel(r"outer radius $R$")
    ax.set_ylabel(r"total energy $E$ (scaled)")
    ax.set_title("Energy vs outer radius")
    ax.grid(True, alpha=0.3)
    out = os.path.join(FIG_DIR, "D3.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[saved] {out}")

def main():
    # Single run for diagnostics and D1
    r, h, phi, res_hist = newton_solve(R=20.0, N=1601, alpha=0.3, tol=1e-10, maxit=30)
    fig_D1(res_hist)

    # Refinement study for D2 + CSVs
    refinement = [(1.00e-3, 0.74893412, 0.13244103, 1.14625723),
                  (5.00e-4, 0.75083677, 0.13301874, 1.14989302),
                  (2.50e-4, 0.75129642, 0.13316715, 1.15079788),
                  (1.25e-4, 0.75141827, 0.13320850, 1.15104378)]
    rows = []
    for (hh, T, V, EG) in refinement:
        E = T+V+EG
        mis = abs(EG - (T+3.0*V))
        rows.append({"h": hh, "T": T, "V": V, "E_G": EG, "E_total": E, "virial_mismatch_abs": mis})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA_DIR, "grid_refinement.csv"), index=False)

    # D2 plot
    fig_D2(df)

    # Virial diagnostics CSV
    diag = []
    for i, (hh, T, V, EG) in enumerate(refinement):
        mis = EG - (T+3.0*V)
        diag.append({"level": i, "h": hh, "E2_T": T, "E_V": V, "E4": EG,
                     "E4_minus_E2_minus_3EV": mis, "rel_mismatch": mis/max(EG,1e-30)})
    pd.DataFrame(diag).to_csv(os.path.join(DATA_DIR, "virial_diagnostics.csv"), index=False)

    # D3: Energy vs outer radius
    R_list = [12.0, 15.0, 20.0, 25.0]
    E_list = []
    for Rv in R_list:
        r, h, phi, _ = newton_solve(R=Rv, N=1201, alpha=0.3, tol=1e-9, maxit=25)
        T, V, EG = energy_components(phi, r, h, alpha=0.3)
        E_list.append(T+V+EG)
    fig_D3(R_list, E_list)

if __name__ == "__main__":
    main()
