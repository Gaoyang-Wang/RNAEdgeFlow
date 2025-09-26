#!/usr/bin/env bash
set -euo pipefail

# RNAEdgeFlow installation helper:
# 1) create conda/mamba env  2) activate  3) ensure executables
echo "[INFO] Creating conda env: rnaedgeflow"
if command -v mamba >/dev/null 2>&1; then
  mamba env create -f env/environment.yml
else
  conda env create -f env/environment.yml
fi

# Activate
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rnaedgeflow

# Make scripts executable
chmod +x rnaedgeflow || true
chmod +x scripts/*.sh || true

echo "[INFO] Installation complete."
echo "[INFO] Activate with: conda activate rnaedgeflow"
echo "[INFO] Try: ./rnaedgeflow --help"
