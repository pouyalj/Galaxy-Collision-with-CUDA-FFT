# Stage-5 GPU box setup (NVIDIA / CUDA)

Bring up a **fresh Ubuntu box with an NVIDIA GPU** for the Stage-5 CUDA work (D7 / §7
prerequisite). Target verified for: **RTX 3070 (8 GB), fresh Ubuntu 22.04 / 24.04.**

What the box needs and *doesn't*:

- **NVIDIA driver, CUDA-12-capable (≥535).** This is the only hard system requirement. Taichi's
  CUDA backend JITs its own PTX and links the driver's `libcuda` — you do **not** need the full
  CUDA toolkit / `nvcc` for the production multigrid path. (`nvidia-smi` will report
  `CUDA Version: 12.x`; that is the *driver's* capability, not an installed toolkit.)
- **Python 3.11 or 3.12** — the repo pins `>=3.11,<3.13` (Taichi wheels lag newer CPython).
  Ubuntu 24.04 ships 3.12; **22.04 ships 3.10, which is too old** — the script below pulls 3.12
  from deadsnakes when needed.
- **`cupy-cuda12x` is optional** — only for the zero-padded FFT *oracle* spot-check. Its wheel
  bundles the CUDA libs, so it needs the driver, not a toolkit. The production solver is multigrid.

> **8 GB VRAM note.** 100M particles is ~3.8 GB of state (§5.2), so the 3070 clears the §7
> 100M-particle gate — approach it via the benchmark ladder (1M → 10M → 30M → 100M). Two caveats:
> (1) the **FFT oracle at 256³** builds a 512³ buffer (plus a transient coordinate-grid spike,
> AGENT.md §11) that will not fit 8 GB — run the oracle as a reduced-N or CPU spot-check, never the
> production path; (2) wiring Taichi's `device_memory_fraction` / `device_memory_GB` into
> `init_backend` is likely the first Stage-5 code task so the 100M run controls its VRAM pool.

---

## Phase 1 — NVIDIA driver (manual; needs a reboot)

```bash
sudo apt update && sudo apt -y full-upgrade
sudo apt install -y ubuntu-drivers-common
ubuntu-drivers devices            # shows the recommended driver for the 3070
sudo ubuntu-drivers install       # installs it (e.g. nvidia-driver-550)
# explicit alternative: sudo apt install -y nvidia-driver-550
sudo reboot
```

After the reboot, confirm the driver is live:

```bash
nvidia-smi   # must list "NVIDIA GeForce RTX 3070", Driver Version >= 535, CUDA Version 12.x
```

If `nvidia-smi` reports *no devices* or the module isn't loaded, **Secure Boot** is most likely
blocking the unsigned kernel module: either enroll the MOK key when prompted during install, or
disable Secure Boot in BIOS, then reboot.

## Phase 2 + 3 — Python env, repo, and the CUDA verification gate (scripted)

With the driver up, clone the repo and run the env script — it installs Python 3.12 (if missing),
creates the venv, installs the project, **verifies Taichi actually resolves the GPU**, and runs the
baseline suite:

```bash
git clone https://github.com/pouyalj/Galaxy-Collision-with-CUDA-FFT.git
cd Galaxy-Collision-with-CUDA-FFT
bash scripts/setup_gpu_env.sh
```

Or do it by hand:

```bash
# Python 3.12 (skip the PPA line on 24.04, which already has 3.12)
sudo add-apt-repository -y ppa:deadsnakes/ppa && sudo apt update     # 22.04 only
sudo apt install -y git build-essential python3.12 python3.12-venv python3.12-dev

python3.12 -m venv .venv && source .venv/bin/activate
pip install -U pip
pip install -e ".[dev]"            # taichi, numpy, h5py, matplotlib, pyyaml + pytest, ruff
# optional FFT-oracle on GPU:  pip install cupy-cuda12x
```

**Verification gate** (the script runs these; here they are to run by hand):

```bash
# 1. Taichi resolves CUDA, not a silent CPU fallback. init_backend WARNs loudly if it fell back
#    (the safety property in tests/test_backend.py), so a run never *looks* like it used a GPU it
#    didn't.
python -c "from galaxy_collision.config import SimConfig; \
from galaxy_collision.sim import init_backend; \
i=init_backend(SimConfig(backend='cuda')); print(i); assert not i['fell_back']"

# 2. The smoke sim on CUDA, and the full baseline suite.
hello-sim --backend cuda
pytest -q
```

Expected: the `init_backend` line reports a device naming the **RTX 3070** with `fell_back: False`,
and the suite is green on the new box. If it prints a `WARNING … fell back to CPU` (resolved arch
`x64`/`arm64`), Taichi can't see CUDA — re-check `nvidia-smi` and that `taichi` imported the GPU
runtime (`python -c "import taichi as ti; ti.init(arch=ti.cuda)"` should not raise a `cuInit` error).

## Phase 4 — Claude Code on the box

So the staged review workflow continues here. The current Linux install is the native installer
(no Node needed):

```bash
curl -fsSL https://claude.ai/install.sh | bash
# open a new shell (or source your profile), then:
claude doctor                              # verify the install
cd ~/Galaxy-Collision-with-CUDA-FFT
claude                                     # first run → OAuth login in the browser
```

(Fallback if you prefer npm: install Node 18+ and `npm i -g @anthropic-ai/claude-code`.)

---

## Quick troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `nvidia-smi`: command not found / no devices | Driver not installed or not loaded — redo Phase 1; check Secure Boot (MOK enrollment). |
| `init_backend` WARNs *fell back to CPU* | Taichi can't reach CUDA: driver missing, or `taichi` import sees no GPU. Confirm `nvidia-smi`, then `ti.init(arch=ti.cuda)`. |
| `cuInit` error from `ti.init(arch=ti.cuda)` | Driver/userspace mismatch — reboot after the driver install; ensure `libcuda.so` is present (`ldconfig -p | grep libcuda`). |
| OOM at 30–100M particles | 8 GB is tight at 100M. Set `device_memory_fraction`/`device_memory_GB` in `ti.init` (Stage-5 wiring), run headless (no X) to free VRAM, and climb the 1M→10M→30M→100M ladder. |
| `pip install -e .` rejects the Python | You're on 3.10 (Ubuntu 22.04 default). Use `python3.12` explicitly (deadsnakes). |

## Where this leaves Stage 5

Once the gate is green, the first code tasks (already logged as deferrals) are: device-resident
state and the `device_memory` knob (RV6), multigrid perf tuning (RV5), and the benchmark matrix
(steps/sec at 1M/10M/30M/100M + per-stage profile) that is the §7 exit gate.
