This repo reproduces the Part-I figures (D1–D3) directly from the underlying PDE
with a fully scripted pipeline, provenance logs, and hard numerical gates.

**Core guarantees**
- All figures and CSVs are generated from **PDE solves only** (no hand numbers).
- Stationarity enforced: final Newton residual **‖F‖₂ ≤ 1e-10** before artifacts.
- D2 is computed from the **true virial defect** \(R_{\rm vir}(h)=|T(h)-V(h)|\)
  using **4th-order** derivatives and **Simpson** integration; the fitted slope
  on \(\log R_{\rm vir}\) vs \(\log h\) is **~2** (second-order scheme).
- Tests gate the pipeline and schemas.

---

## Quick start

```bash
python codes/PartI_grid_solver_GF.py \
  --model doublewell --alpha 0.3 \
  --R_base 20 --h_targets "0.02,0.01,0.005" \
  --iters 40000 --target_vrel 1e-4 \
  --tol 1e-10 --kink 1 --root .

pytest -q
```
Artifacts (always rewritten):

data/warm_start_history.csv, data/newton_history.csv

data/grid_refinement.csv, data/virial_diagnostics.csv

data/convergence_table.csv

figures/D1.png, figures/D2.png, figures/D3.png

figures/kink_profile.png, figures/cs2_band.png, figures/E_lambda_curve.png

logs/provenance.json

What success looks like

pytest -q → all green.

data/newton_history.csv (or newton_accepted_history.csv) final residual_L2 ≤ 1e-10.

data/virial_true_vs_h.csv fitted slope 
𝑝
≈
2
p≈2 (see “Paper-style figures” below).
Paper-style figures (accepted-Newton D1, true virial D2, plateau D3)
```
Use the builder to regenerate the paper versions with higher-order
numerics and accepted-step Newton:

# Fast (sanity)
python codes/make_paper_figs.py --fast

# Full (referee copy)
python codes/make_paper_figs.py --L 20 --hs "0.02,0.01,0.005" \
  --gf_iters 60000 --newton_steps 16 --newton_tol 1e-10
```

Outputs:

figures/D1.png — accepted Newton steps only (solid: residual; dashed: CPU).

figures/D2.png — true 
𝑅
v
i
r
=
∣
𝑇
−
𝑉
∣
R
vir
	​

=∣T−V∣ vs 
ℎ
h with dashed 
∝
ℎ
2
∝h
2
;
inputs in data/virial_true_vs_h.csv, slope printed in the title.

figures/D3.png — 
𝐸
(
𝑅
)
E(R) plateau with dashed 
𝐸
∞
E
∞
	​

.

data/newton_accepted_history.csv — accepted Newton iterates (iter, residual_L2, cpu).

If the D2 slope is not ~2, increase the domain --L or add a finer --hs..

Commands used most often
```
# Re-run core PDE pipeline for the kink
python codes/PartI_grid_solver_GF.py --model doublewell --R_base 20 \
  --h_targets "0.02,0.01,0.005" --iters 40000 --target_vrel 1e-4 \
  --tol 1e-10 --kink 1 --root .

# Optional: quadratic sanity (vacuum)
python codes/PartI_grid_solver_GF.py --model quadratic --R_base 20 \
  --h_targets "0.02,0.01,0.005" --iters 40000 --tol 1e-10 --root .

# Paper-style figures
python codes/make_paper_figs.py --fast
```

Repo layout (relevant files)
codes/
  gf/
    solver1d.py          # semi-implicit GF solver with kink pin & warm-start
    highorder.py         # 4th-order derivative stencils + Simpson integrator
  PartI_grid_solver_GF.py# driver that writes CSVs & D1–D3 (test-friendly)
  make_paper_figs.py     # paper-style figures (accepted Newton, true virial)
data/
  *.csv                  # reproducible tables (see list above)
figures/
  *.png                  # all figures, rewritten each run
logs/
  provenance.json        # timestamp, model, grids, tolerances, env versions
tests/
  ...                    # schema + stationarity gates

