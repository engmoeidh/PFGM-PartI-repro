import csv, os

EXPECTED = {
    "data/convergence_table.csv": ["N","h","residual_L2","observed_order_p"],
    "data/grid_refinement.csv":   ["h","T","V","E_G","E_total","virial_mismatch_abs"],
    "data/virial_diagnostics.csv":["level","h","E2_T","E_V","E4","E4_minus_E2_minus_3EV","rel_mismatch"],
}

def test_csv_headers_match_expected():
    missing, wrong = [], []
    for path, cols in EXPECTED.items():
        if not os.path.isfile(path):
            missing.append(path); continue
        with open(path, newline="", encoding="utf-8") as f:
            header = next(csv.reader(f), [])
        if header != cols:
            wrong.append((path, header, cols))
    assert not missing, f"Missing CSVs: {missing}"
    assert not wrong, "Mismatched headers: " + "; ".join(f"{p}: {h} != {e}" for p,h,e in wrong)
