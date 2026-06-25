# Galaxy Collision

A cross-platform, research-grade **Particle-Mesh (PM) gravitational N-body** simulation of a
Milky Way × Andromeda galaxy collision. Two spiral galaxies are placed apart, pushed toward each
other, and evolved through the collision on a 256³ density grid — scaling to **100M particles** on
a single GPU.

[![CI](https://github.com/pouyalj/Galaxy-Collision-with-CUDA-FFT/actions/workflows/ci.yml/badge.svg)](https://github.com/pouyalj/Galaxy-Collision-with-CUDA-FFT/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg)

Originally a 2020 University of Toronto Scarborough PHYD57 project
([paper PDF](https://static1.squarespace.com/static/5d5dcd310b9b0100013bcbe1/t/5f05fb53b5c64717ff32882e/1594227543673/PHYD57_NBody-Project_No1-2020.pdf)),
it is being rebuilt as **one portable source** that runs on CPU, NVIDIA CUDA, and Apple-Silicon
GPU from the same kernels via [Taichi](https://www.taichi-lang.org/).

## Highlights

- **One source, three backends** — the same `@ti.kernel` code compiles to CPU, CUDA, and (in
  progress) Apple Metal; pick at runtime with `--backend`.
- **Research-grade physics** — Cloud-In-Cell deposit/gather, a symplectic kick-drift-kick
  leapfrog, force softening, central black holes carrying real grid mass, and **open
  (isolated) boundary conditions** that fix the periodic-wrap artifact of a plain FFT.
- **Pluggable Poisson solver** — an open-BC geometric **multigrid** (the portable default) plus a
  zero-padded **FFT oracle** (NumPy, or CuPy/cuFFT on NVIDIA) as a spectrally-exact ground truth.
- **Fast at scale** — a 100M-particle, 400-Myr collision runs in **~3.5 minutes** on an 8 GB
  RTX 3070 (see [Performance](#performance)).
- **Validated** — a stable Plummer sphere, a two-body Kepler orbit, multigrid-vs-oracle agreement,
  conservation diagnostics, and CPU↔CUDA determinism, all in CI.

> **Status:** Stages 0–5 complete (correct CPU reference → paper reproduction → CUDA scale-up).
> Apple/Metal (Stage 6) and the realtime viewer/movie output (Stage 7) are next. See
> [`AGENT.md`](AGENT.md) for the full architecture and staged roadmap.

## Installation

Requires **Python 3.11 or 3.12** (Taichi has no 3.13 wheels yet).

```bash
git clone https://github.com/pouyalj/Galaxy-Collision-with-CUDA-FFT.git
cd Galaxy-Collision-with-CUDA-FFT
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"          # core + pytest/ruff
```

Taichi selects the GPU at runtime — **no CUDA toolkit install is required** for the simulation
itself (Taichi ships its own JIT). For an NVIDIA box, `scripts/setup_gpu_env.sh` provisions the
venv and verifies the CUDA backend is live; see [`docs/gpu_setup.md`](docs/gpu_setup.md).

Optional extras:

```bash
pip install -e ".[cuda]"         # CuPy/cuFFT GPU path for the FFT validation oracle
```

## Quickstart

```bash
hello-sim                                   # no-op smoke run — confirms the toolchain + backend
galaxy-sim --config configs/plummer.yaml    # a real PM N-body run (Plummer sphere)
galaxy-sim --config configs/paper_4v.yaml --backend cuda   # the 4v collision on the GPU
paper-repro --config configs/paper_4v.yaml  # collision → paper figures (heavy; see the config)
galaxy-bench --n 10000000 --backend cuda    # benchmark the force chain at 10M particles
```

Every run prints the **device it actually ran on** and warns loudly if a requested GPU backend
silently fell back to CPU.

## Usage

| Command | What it does |
|---|---|
| `hello-sim [--config C] [--backend B]` | Minimal end-to-end smoke run (toolchain + backend check). |
| `galaxy-sim --config C [--backend B]` | Run the full PM N-body pipeline from a YAML config; writes HDF5/npz snapshots + diagnostics. |
| `paper-repro --config C [--backend B]` | Run a two-galaxy collision and render the paper's density-projection + Sun-like-tracer figures. |
| `galaxy-bench --n N [--grid G] [--solver S] [--backend B]` | Throughput + per-stage profile of the force chain. |

Runs are driven by **YAML configs** (validated against a schema). Ships with `configs/`:
`plummer.yaml` (equilibrium test), `paper_4v.yaml` / `paper_2v.yaml` (the two collision speeds),
and `smoke.yaml`. Key knobs: `backend` (`cpu`/`cuda`/`metal`), `n_particles` (or `particle_mass`),
`grid_size`, `dt`, `steps`, `ic_preset`, `solver` (`multigrid`/`fft`), `output_cadence`, `seed`.

## Performance

On an 8 GB NVIDIA RTX 3070 (256³ grid, multigrid solver):

| N particles | throughput | peak VRAM |
|---|---|---|
| 1M | 33 steps/s | 0.8 GB |
| 10M | 25 steps/s | 1.3 GB |
| 30M | 16 steps/s | 2.5 GB |
| 100M | 7 steps/s | 6.5 GB |

A 100M-particle, 800-step (400-Myr) collision completes in ~3.5 minutes. Full benchmark, the
per-stage breakdown, and the tuning that got there (adaptive warm-start multigrid cycling) are in
[`docs/performance.md`](docs/performance.md).

## How it works

Per timestep, all state stays resident on the device:

```
kick(½dt) → drift(dt) → deposit (CIC) → Poisson solve (∇²Φ = 4πGρ) → g = −∇Φ → gather → kick(½dt)
```

Mass is binned onto the grid by Cloud-In-Cell, Poisson's equation is solved with open boundary
conditions (multigrid, or the FFT oracle), the potential is differentiated into forces, and
particles are advanced by a symplectic leapfrog. The Poisson solve sits behind a small interface
so the solver is swappable. Units are (kpc, Myr, M☉) with a single unit-tested gravitational
constant. The full design, decision record, and physics rationale live in [`AGENT.md`](AGENT.md).

## Repository layout

```
src/galaxy_collision/
  config.py        # SimConfig schema + YAML load/dump
  units.py         # (kpc, Myr, M☉) unit system + derived G
  data.py          # Structure-of-Arrays particle + grid state
  ic.py            # two-galaxy + Plummer initial conditions
  deposit.py       # CIC deposit, −∇Φ, force gather (Taichi kernels)
  integrator.py    # KDK leapfrog + direct softened force
  solver/          # PoissonSolver interface: multigrid (default) + fft_oracle
  diagnostics*.py  # energy/momentum/tracer (host fp64 + device reductions)
  io.py            # HDF5/npz snapshots
  bench.py         # per-stage profiler + benchmark (galaxy-bench)
  sim.py           # orchestration + CLI entry points
  viz/             # paper-reproduction figures (realtime viewer = Stage 7)
configs/  ·  tests/  ·  docs/  ·  legacy/   # YAML configs · pytest suite · docs · the 2020 CUDA source
```

## Development

```bash
ruff check .                    # lint (what CI runs)
pytest                          # full test suite
GALAXY_TEST_ARCH=cuda pytest    # re-run the suite on the GPU backend
```

CI (GitHub Actions) runs ruff + pytest + a CPU smoke run on every push. The test suite covers the
units/G, IC generation, the solvers (convergence + analytic anchor + oracle agreement), the
integrator (Kepler), a Plummer-equilibrium end-to-end gate, conservation diagnostics, and
CPU↔CUDA determinism. Contributor setup is in [`docs/development.md`](docs/development.md).

## Documentation

| Doc | Contents |
|---|---|
| [`AGENT.md`](AGENT.md) | The project bible — architecture, decision record, staged plan, review log. |
| [`docs/development.md`](docs/development.md) | Contributor quickstart. |
| [`docs/gpu_setup.md`](docs/gpu_setup.md) | Provisioning a CUDA workstation. |
| [`docs/performance.md`](docs/performance.md) | Benchmark matrix + per-stage profile (Stage 5). |
| [`docs/paper_reproduction.md`](docs/paper_reproduction.md) | Reproducing the 2020 paper's figures. |
| [`docs/stage5_plan.md`](docs/stage5_plan.md) | The CUDA scale-up plan. |

## Roadmap

Stages 0–5 are done: scaffold & CI → units & data model → initial conditions → correct CPU
reference → validation & FFT oracle & paper reproduction → **CUDA scale-up to 100M**. Up next:
**Stage 6** (Apple/Metal as a first-class GPU target) and **Stage 7** (realtime Taichi GGUI viewer
+ batch→movie output). See [`AGENT.md` §7](AGENT.md) for per-stage detail and exit criteria.

## License & origin

Apache-2.0 — see [`LICENSE`](LICENSE). Originally *"Galaxy Collisions With CUDA and FFT"* (Adams,
Lajevardi, LeBlanc, Movaseghi; PHYD57, 2020); the original CUDA source is preserved under
[`legacy/`](legacy/).
