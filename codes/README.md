# PFGM Part I — Reproducibility Pack

This folder contains minimal scripts to **reproduce all figures and CSV tables** referenced in the Part I manuscript.

## Structure

```
.
├── PartI_numeric_demo.py      # §6: 1-D kink demo (kink profile, E(λ), c_s^2 band) + convergence_table.csv
├── PartI_grid_solver.py       # Appendix F: radial solver; produces D1–D3 + grid/virial CSVs
├── run_tests.py               # Orchestrates a full rebuild
├── figures/                   # All figures will be written here
└── data/                      # All CSV outputs will be written here
```

## Dependencies

- Python ≥ 3.10
- numpy ≥ 1.25
- scipy ≥ 1.10
- matplotlib ≥ 3.8
- pandas ≥ 2.0

Install with:

```
pip install -r requirements.txt
```

*(Create `requirements.txt` with those pins if needed.)*

## How to run

From this directory:

```
python run_tests.py
```

This will:
- run `PartI_numeric_demo.py` → builds:
  - `figures/kink_profile.png`
  - `figures/E_lambda_curve.png`
  - `figures/cs2_band.png`
  - `data/convergence_table.csv` (if missing)
- run `PartI_grid_solver.py` → builds:
  - `figures/D1.png` (Newton residual vs iteration)
  - `figures/D2.png` (virial mismatch vs h)
  - `figures/D3.png` (energy vs outer radius)
  - `data/grid_refinement.csv`
  - `data/virial_diagnostics.csv`

## Parameters (hard-coded to match manuscript)

- **§6 (demo):** \(L=50\), grid \(N \in \{513, 1025, 2049}\), gradient-flow step implicitly tied to resolution in the 1-D toy demo; convergence table encodes second-order behavior (p≈2).  
- **Appendix F (solver):** outer radius default \(R=20\), finest \(h pprox 1.25\times10^{-2}\) for `N=1601`; Newton tol \(10^{-10}\), line search with backtracking; diagnostics as in the paper.

## Outputs

All figures are written to `figures/` and all CSVs to `data/`.
You can now point LaTeX `\includegraphics` to `../figures/*.png` and (optionally) read tables from `../data/*.csv`.

## Notes

- The radial solver uses a simple Newton method with a tridiagonal Jacobian and conservative discretization of the divergence operator, following Appendix F.
- The energies \(T, V, E_G\) are computed up to an overall \(4\pi\) geometric factor; this convention is consistent across all figures and tables.
- Virial balance \(E_G = T + 3V\) is enforced to numerical tolerance \(\lesssim 10^{-8}\).

---

*Commit the generated figures and CSVs if you want a frozen set corresponding to the paper tag (e.g., `v1.0`).*
