#!/usr/bin/env python3
"""
run_tests.py
End-to-end generator:
 - runs the 1-D demo and grid solver
 - ensures all CSVs are (re)built
 - confirms observed order ~ 2 for the §6 convergence table
"""
import os, sys, subprocess
import pandas as pd

<<<<<<< HEAD
# paths
THIS = os.path.abspath(__file__)              # .../PFGM-PartI-repro/codes/run_tests.py
ROOT = os.path.dirname(os.path.dirname(THIS)) # .../PFGM-PartI-repro
FIG_DIR = os.path.join(ROOT, "figures")
DATA_DIR = os.path.join(ROOT, "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

def run_script(rel_path_from_root):
    path = os.path.join(ROOT, rel_path_from_root)
    print(f"[run] {rel_path_from_root}")
    # use the current interpreter so venv is respected
    subprocess.check_call([sys.executable, path], cwd=ROOT)
=======
ROOT = os.path.dirname(os.path.abspath(__file__))  # .../PFGM-PartI-repro/codes
REPO = os.path.dirname(ROOT)                       # .../PFGM-PartI-repro
FIG_DIR = os.path.join(REPO, "figures")
DATA_DIR = os.path.join(REPO, "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

def run_script(rel_path):
    path = os.path.join(REPO, rel_path)
    print(f"[run] {rel_path}")
    try:
        subprocess.check_call([sys.executable, path], cwd=REPO)
    except subprocess.CalledProcessError as e:
        print(f"[error] {rel_path} failed with exit code {e.returncode}")
        raise
>>>>>>> 4a3bb852e248917ee69643173318afd6bb76893f

def check_convergence_table():
    path = os.path.join(DATA_DIR, "convergence_table.csv")
    if not os.path.exists(path):
<<<<<<< HEAD
        print("[warn] convergence_table.csv not found (it should be created by PartI_numeric_demo.py)")
=======
        print("[warn] convergence_table.csv not found; it should have been created by PartI_numeric_demo.py")
>>>>>>> 4a3bb852e248917ee69643173318afd6bb76893f
        return
    df = pd.read_csv(path)
    vals = df["observed_order_p"].values
    if len(vals) >= 3 and str(vals[1]).strip() != "" and str(vals[2]).strip() != "":
        p12 = float(vals[1]); p23 = float(vals[2])
        print(f"[info] observed orders: p12={p12}, p23={p23}")
        assert 1.8 <= p12 <= 2.2 and 1.8 <= p23 <= 2.2, "Observed orders deviate from 2"
    else:
<<<<<<< HEAD
        print("[warn] observed_order_p not populated; leaving as-is")
=======
        print("[warn] observed_order_p not populated as expected; continuing")
>>>>>>> 4a3bb852e248917ee69643173318afd6bb76893f

def main():
    run_script(os.path.join("codes", "PartI_numeric_demo.py"))
    run_script(os.path.join("codes", "PartI_grid_solver.py"))
    check_convergence_table()
    print("[ok] all artifacts generated")

if __name__ == "__main__":
    main()
