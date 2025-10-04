import numpy as np
import pandas as pd

def finite(series):
    a = series.to_numpy()
    return np.isfinite(a).all()

def test_convergence_rate_from_table():
    df = pd.read_csv("data/convergence_table.csv")  # N,h,residual_L2,observed_order_p
    for col in ("N","h","residual_L2","observed_order_p"):
        assert col in df.columns, f"Missing column {col} in convergence_table.csv"
    assert finite(df["h"]) and finite(df["residual_L2"]) and finite(df["observed_order_p"]), "NaNs in convergence table"
    # Accept either direct p or reconstructed p from residuals as >= ~2 (allow small slack)
    p_med_reported = float(np.nanmedian(df["observed_order_p"]))
    # Reconstruction (optional cross-check): p ~ log(e_i/e_{i+1}) / log(h_i/h_{i+1})
    df2 = df.dropna(subset=["h","residual_L2"]).sort_values("h", ascending=False)
    if len(df2) >= 3:
        num = np.log(df2["residual_L2"].values[:-1]/df2["residual_L2"].values[1:])
        den = np.log(df2["h"].values[:-1]/df2["h"].values[1:])
        p_rec = num/den
        p_med_rec = float(np.nanmedian(p_rec))
        assert p_med_rec > 1.8, f"Reconstructed convergence order too low: median={p_med_rec:.2f}"
    assert p_med_reported > 1.8, f"Reported convergence order too low: median={p_med_reported:.2f}"

def test_grid_refinement_virial_tight():
    g = pd.read_csv("data/grid_refinement.csv")  # h,T,V,E_G,E_total,virial_mismatch_abs
    for col in ("h","virial_mismatch_abs"):
        assert col in g.columns, "grid_refinement.csv must have h, virial_mismatch_abs"
    assert finite(g["h"]) and finite(g["virial_mismatch_abs"]), "NaNs in grid refinement"
    med = float(np.nanmedian(g["virial_mismatch_abs"]))
    assert med < 1e-3, f"Virial mismatch too large: median={med:.2e}"

def test_virial_diagnostics_rel_mismatch_tight():
    v = pd.read_csv("data/virial_diagnostics.csv")  # level,h,E2_T,E_V,E4,...,rel_mismatch
    assert all(c in v.columns for c in ("level","h","rel_mismatch")), "virial_diagnostics.csv must have level,h,rel_mismatch"
    assert finite(v["rel_mismatch"]), "NaNs in virial diagnostics"
    med = float(np.nanmedian(v["rel_mismatch"]))
    assert med < 1e-3, f"Virial relative mismatch too large: median={med:.2e}"
