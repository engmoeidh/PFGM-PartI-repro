import numpy as np
import pandas as pd

def finite(series):
    a = series.to_numpy()
    return np.isfinite(a).all()

def test_convergence_rate_from_table():
    df = pd.read_csv("data/convergence_table.csv")  # N,h,residual_L2,observed_order_p
    for col in ("N","h","residual_L2","observed_order_p"):
        assert col in df.columns, f"Missing column {col} in convergence_table.csv"

    # Basic sanity (allow NaNs only in observed_order_p)
    assert finite(df["h"]), "NaNs in h"
    assert finite(df["residual_L2"]), "NaNs in residual_L2"

    # Primary: reconstruct order from residuals vs h
    df2 = df.dropna(subset=["h","residual_L2"]).sort_values("h", ascending=False)
    assert len(df2) >= 3, "Need at least 3 refinement levels to assess convergence"
    num = np.log(df2["residual_L2"].values[:-1]/df2["residual_L2"].values[1:])
    den = np.log(df2["h"].values[:-1]/df2["h"].values[1:])
    p_rec = num/den
    p_med_rec = float(np.nanmedian(p_rec))

    # Secondary: use reported observed_order_p if available (ignore NaNs)
    if "observed_order_p" in df.columns:
        p_reported = df["observed_order_p"].dropna().to_numpy()
        p_med_rep = float(np.nanmedian(p_reported)) if p_reported.size else np.nan
    else:
        p_med_rep = np.nan

    # Pass if either meets threshold (robust to one noisy column)
    assert (p_med_rec > 1.8) or (np.isfinite(p_med_rep) and p_med_rep > 1.8), \
        f"Convergence too low: reconstructed median={p_med_rec:.2f}, reported median={p_med_rep}"

def test_grid_refinement_virial_tight():
    g = pd.read_csv("data/grid_refinement.csv")  # h,T,V,E_G,E_total,virial_mismatch_abs
    for col in ("h","virial_mismatch_abs"):
        assert col in g.columns, "grid_refinement.csv must have h, virial_mismatch_abs"
    assert np.isfinite(g["h"]).all() and np.isfinite(g["virial_mismatch_abs"]).all(), "NaNs in grid refinement"
    med = float(np.nanmedian(g["virial_mismatch_abs"]))
    assert med < 1e-3, f"Virial mismatch too large: median={med:.2e}"

def test_virial_diagnostics_rel_mismatch_tight():
    v = pd.read_csv("data/virial_diagnostics.csv")  # level,h,E2_T,E_V,E4,...,rel_mismatch
    assert all(c in v.columns for c in ("level","h","rel_mismatch")), "virial_diagnostics.csv must have level,h,rel_mismatch"
    assert np.isfinite(v["rel_mismatch"]).all(), "NaNs in virial diagnostics"
    med = float(np.nanmedian(v["rel_mismatch"]))
    assert med < 1e-3, f"Virial relative mismatch too large: median={med:.2e}"
