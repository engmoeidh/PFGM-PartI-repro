#!/usr/bin/env bash
set -Eeuo pipefail
trap 'ec=$?; echo ""; echo "✖ Error (exit $ec). See logs/run_all.log"; echo ""; read -n1 -r -p "Press any key to close..."; exit $ec' ERR

mkdir -p logs
exec > >(tee "logs/run_all.log") 2>&1
set -x

# --- ALWAYS start from the repo root, robustly ---
if git rev-parse --show-toplevel >/dev/null 2>&1; then
  cd "$(git rev-parse --show-toplevel)"
else
  # fallback to script directory/..
  cd "$(dirname "$0")/.."
fi

# 0) Python venv (Windows Git Bash)
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

# 1) ANALYSIS
[ -f codes/PartI_grid_solver.py ]  && python codes/PartI_grid_solver.py  || true
[ -f codes/PartI_numeric_demo.py ] && python codes/PartI_numeric_demo.py || true

# 2) TESTS
if command -v pytest >/dev/null 2>&1; then
  pytest -q || echo "pytest failed; continuing to build PDF so referees can still read."
fi

# 3) CHECKSUMS
mkdir -p data
: > logs/data_sha256.txt
for f in data/convergence_table.csv data/grid_refinement.csv data/virial_diagnostics.csv; do
  if [ -f "$f" ]; then sha256sum "$f" >> logs/data_sha256.txt; else echo "MISSING $f" >> logs/data_sha256.txt; fi
done

# 4) BUILD PAPER (resilient)
if [ ! -d paper ]; then
  echo "✖ No 'paper/' directory found in $(pwd)."; read -n1 -r -p "Press any key to close..."; exit 1
fi

# Try to ensure microtype; if not, compile a temp file without it.
need_microtype=0
kpsewhich microtype.sty >/dev/null 2>&1 || need_microtype=1
if [ $need_microtype -eq 1 ]; then
  initexmf --set-config-value [MPM]AutoInstall=1 || true
  initexmf --update-fndb || true
  ( mpm --install=microtype || mpm --admin --install=microtype || "/c/Users/HomePC/AppData/Local/Programs/MiKTeX/miktex/bin/x64/mpm.exe" --install=microtype ) || true
  initexmf --update-fndb || true
fi

cd paper
latexmk -C || true
if kpsewhich microtype.sty >/dev/null 2>&1; then
  latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main.tex
else
  echo "[TeX] microtype still missing — compiling a temporary file without microtype."
  cp main.tex main__tmp_no_microtype.tex
  sed -i "s/^[[:space:]]*\\\\usepackage[{]microtype[}]/% DISABLED microtype/gI" main__tmp_no_microtype.tex
  latexmk -C || true
  latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main__tmp_no_microtype.tex
  [ -f main__tmp_no_microtype.pdf ] && cp -f main__tmp_no_microtype.pdf main.pdf
fi
cd ..

echo ""
echo "✔ Pipeline finished."
echo "   PDF:                 paper/main.pdf"
echo "   Log (full details):  logs/run_all.log"
echo "   Data checksums:      logs/data_sha256.txt"
echo ""
read -n1 -r -p "Press any key to close..."
