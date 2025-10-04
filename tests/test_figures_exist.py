import os

def test_required_figures_exist_and_nonempty():
    req = [
        "figures/D1.png",
        "figures/D2.png",
        "figures/D3.png",
        "figures/kink_profile.png",
        "figures/E_lambda_curve.png",
        "figures/cs2_band.png",
    ]
    missing = [p for p in req if not os.path.isfile(p)]
    zero = [p for p in req if os.path.isfile(p) and os.path.getsize(p)==0]
    assert not missing, f"Missing figures: {missing}"
    assert not zero, f"Zero-size figures: {zero}"
