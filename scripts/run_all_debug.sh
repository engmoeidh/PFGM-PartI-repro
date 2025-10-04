#!/usr/bin/env bash
# End-to-end reproducible run for PFGM Part I on Windows Git Bash.
# - Boots/uses .venv
# - Runs your analysis scripts
# - Runs pytest smoke tests
# - Hashes data CSVs
# - Builds paper PDF with MiKTeX auto-install bootstrap
# - NEVER auto-closes: shows errors and waits for a keypress

set -Eeuo pipefail
trap 'ec=$?; echo ""; echo "? Error (exit $ec). See logs/run_all.log"; echo ""; read -n1 -r -p "Press any key to close..."; exit $ec' ERR

mkdir -p logs
exec > >(tee "logs/run_all.log") 2>&1
set -x

# Always operate from repo root
cd "$(dirname "$0")/.."

# 0) Python venv (Windows Git Bash)
if [ -z "${VIRTUAL_ENV-}" ]; then
  if [ -f ".venv/Scripts/activate" ]; then
    # shellcheck disable=SC1091
    source .venv/Scripts/activate
  else
    python -m venv .venv
    # shellcheck disable=SC1091
    source .venv/Scripts/activate
    python -m pip install --upgrade pip
    pip install -r requirements.txt
  fi
fi

# 1) ANALYSIS  call your actual scripts if present
[ -f codes/PartI_grid_solver.py ]  && python codes/PartI_grid_solver.py  || true
[ -f codes/PartI_numeric_demo.py ] && python codes/PartI_numeric_demo.py || true

# 2) TESTS [G make sure at least the smoke test runs
if command -v pytest >/dev/null 2>&1; then
  pytest -q || echo "pytest failed; continuing to build PDF so referees can still read."
fi

# 3) DATA INTEGRITY [G compute SHA256 for the three key CSVs (if present)
mkdir -p data
: > logs/data_sha256.txt
for f in data/convergence_table.csv data/grid_refinement.csv data/virial_diagnostics.csv; do
  if [ -f "$f" ]; then
    sha256sum "$f" >> logs/data_sha256.txt
  else
    echo "MISSING $f" >> logs/data_sha256.txt
  fi
done

# 4) TEX BOOTSTRAP — ensure MiKTeX can install missing packages (e.g., microtype)
# Enable on-the-fly installs (no-op if already enabled)
initexmf --set-config-value [MPM]AutoInstall=1 || true
initexmf --update-fndb || true

# If microtype isn’t present, try several installers
if ! kpsewhich microtype.sty >/dev/null 2>&1; then
  ( mpm --install=microtype || mpm --admin --install=microtype || "/c/Users/HomePC/AppData/Local/Programs/MiKTeX/miktex/bin/x64/mpm.exe" --install=microtype ) || true
  initexmf --update-fndb || true
fi

# 5) BUILD PAPER — latexmk (MiKTeX)
cd paper
latexmk -C
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main.tex
cd ..

echo ""
echo "✔ Pipeline finished."
echo "   PDF:                 paper/main.pdf"
echo "   Log (full details):  logs/run_all.log"
echo "   Data checksums:      logs/data_sha256.txt"
echo ""
read -n1 -r -p "Press any key to close..."
