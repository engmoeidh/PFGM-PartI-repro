#!/usr/bin/env bash
set -Eeuo pipefail
mkdir -p logs
exec > >(tee "logs/run_all_ci.log") 2>&1
set -x
cd "$(dirname "$0")/.."

# venv
if [ -z "${VIRTUAL_ENV-}" ]; then
  if [ -f ".venv/Scripts/activate" ]; then
    source .venv/Scripts/activate
  else
    python -m venv .venv
    source .venv/Scripts/activate
    python -m pip install --upgrade pip
    pip install -r requirements.txt
  fi
fi

# analysis
[ -f codes/PartI_grid_solver.py ]  && python codes/PartI_grid_solver.py
[ -f codes/PartI_numeric_demo.py ] && python codes/PartI_numeric_demo.py

# tests
if command -v pytest >/dev/null 2>&1; then
  pytest -q
fi

# checksums
: > logs/data_sha256.txt
for f in data/convergence_table.csv data/grid_refinement.csv data/virial_diagnostics.csv; do
  if [ -f "$f" ]; then sha256sum "$f" >> logs/data_sha256.txt; else echo "MISSING $f" >> logs/data_sha256.txt; fi
done

# TeX bootstrap & build
initexmf --set-config-value [MPM]AutoInstall=1 || true
initexmf --update-fndb || true
if ! kpsewhich microtype.sty >/dev/null 2>&1; then
  ( mpm --install=microtype || mpm --admin --install=microtype || "/c/Users/HomePC/AppData/Local/Programs/MiKTeX/miktex/bin/x64/mpm.exe" --install=microtype ) || true
  initexmf --update-fndb || true
fi
cd paper
latexmk -C
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main.tex
