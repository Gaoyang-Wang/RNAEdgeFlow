#!/usr/bin/env bash
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if ! grep -q "${REPO_DIR}" "${HOME}/.bashrc" 2>/dev/null; then
  echo "export PATH=\"${REPO_DIR}:\$PATH\"" >> "${HOME}/.bashrc"
  echo "[INFO] Added ${REPO_DIR} to PATH. Reopen shell or 'source ~/.bashrc'."
else
  echo "[INFO] PATH already contains ${REPO_DIR}."
fi
