#!/usr/bin/env bash
# Stage-5 GPU environment setup (see docs/gpu_setup.md).
#
# Run AFTER the NVIDIA driver is installed and the box has rebooted (Phase 1), from the repo root:
#     cd Galaxy-Collision-with-CUDA-FFT && bash scripts/setup_gpu_env.sh
#
# Installs Python 3.12 (if missing), a .venv, the project (+dev extras), then VERIFIES that Taichi
# resolves the CUDA backend (not a silent CPU fallback) and runs the baseline suite. Idempotent-ish:
# re-running recreates the venv from scratch.
set -euo pipefail

cd "$(dirname "$0")/.."   # repo root, regardless of where it's invoked from

echo "==> 0. Driver sanity (Phase 1 must be done + rebooted)"
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "ERROR: nvidia-smi not found. Install the NVIDIA driver and reboot first" >&2
  echo "       (see docs/gpu_setup.md, Phase 1)." >&2
  exit 1
fi
nvidia-smi || { echo "ERROR: nvidia-smi failed — driver not loaded (Secure Boot? reboot?)." >&2; exit 1; }

echo "==> 1. System build deps + Python 3.12 (repo pins >=3.11,<3.13; 22.04 ships 3.10)"
sudo apt-get update
sudo apt-get install -y git build-essential curl ca-certificates
if ! command -v python3.12 >/dev/null 2>&1; then
  echo "    python3.12 missing — adding deadsnakes PPA"
  sudo apt-get install -y software-properties-common
  sudo add-apt-repository -y ppa:deadsnakes/ppa
  sudo apt-get update
  sudo apt-get install -y python3.12 python3.12-venv python3.12-dev
else
  sudo apt-get install -y python3.12-venv || true
fi

echo "==> 2. Virtualenv + project (editable, with dev extras: pytest, ruff)"
rm -rf .venv
python3.12 -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"

echo "==> 3. Verify Taichi resolves the GPU (not a silent CPU fallback)"
python - <<'PY'
from galaxy_collision.config import SimConfig
from galaxy_collision.sim import init_backend
info = init_backend(SimConfig(backend="cuda"))
print(info)
assert not info["fell_back"], (
    "Taichi fell back to CPU — the CUDA backend is not available. Check nvidia-smi and that "
    "`python -c \"import taichi as ti; ti.init(arch=ti.cuda)\"` does not raise a cuInit error."
)
print("OK: Taichi CUDA backend is live on:", info["device"])
PY

echo "==> 4. Baseline suite"
python -m pytest -q

cat <<'DONE'

GPU environment ready.
  - Activate later with:   source .venv/bin/activate
  - Optional FFT oracle on GPU:   pip install cupy-cuda12x   (8 GB: oracle at reduced N only)
  - Next: install Claude Code (docs/gpu_setup.md, Phase 4), then start Stage-5 work.
DONE
