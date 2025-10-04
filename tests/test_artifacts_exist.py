import os

REQUIRED = [
    "data/convergence_table.csv",
    "data/grid_refinement.csv",
    "data/virial_diagnostics.csv",
    "paper/main.tex",
]

def test_required_files_exist_and_nonempty():
    missing = [p for p in REQUIRED if not os.path.isfile(p)]
    assert not missing, f"Missing required files: {missing}"
    small = [p for p in REQUIRED if os.path.getsize(p) == 0]
    assert not small, f"Empty files detected: {small}"
