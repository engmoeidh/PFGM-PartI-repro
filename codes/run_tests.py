#!/usr/bin/env python3
"""
run_tests.py
End-to-end generator:
 - runs the 1-D demo and grid solver
 - ensures all CSVs are (re)built
 - confirms observed order ~ 2 for the §6 convergence table
"""
import os, math
import numpy as np
import pandas as pd
import subprocess

FIG_DIR = os.path.join(os.getcwd(), "figures")
DATA_DIR = os.path.join(os.getcwd(), "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

def run_script(path):
    print(f"[run] {path}")
    subprocess.check_call(["python", path])

def check_convergence_table():
    path = os.path.join(DATA_DIR, "convergence_table.csv")
    if not os.path.exists(path):
        print("[warn] convergence_table.csv not found; will be created by PartI_numeric_demo.py")
    else:
        df = pd.read_csv(path)
        # simple check: observed p ~ 2
        vals = df["observed_order_p"].values
        if str(vals[1]).strip() != "" and str(vals[2]).strip() != "":
            p12 = float(vals[1]); p23 = float(vals[2])
            print(f"[info] observed orders: p12={p12}, p23={p23}")
            assert 1.8 <= p12 <= 2.2 and 1.8 <= p23 <= 2.2, "Observed orders deviate from 2"

def main():
    run_script("PartI_numeric_demo.py")
    run_script("PartI_grid_solver.py")
    check_convergence_table()
    print("[ok] all artifacts generated")

if __name__ == "__main__":
    main()
