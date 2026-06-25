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

- **One source, three backends** — the same `@ti.kernel` code compiles to CPU, NVIDIA CUDA, and
  Apple Metal; pick at runtime with `--backend`. (Metal has no hardware fp64, so the device
  reductions use a Kahan-compensated fp32 path there — D6.)
- **Research-grade physics** — Cloud-In-Cell deposit/gather, a symplectic kick-drift-kick
  leapfrog, force softening, central black holes carrying real grid mass, and **open
  (isolated) boundary conditions** that fix the periodic-wrap artifact of a plain FFT.
- **Pluggable Poisson solver** — an open-BC geometric **multigrid** (the portable default) plus a
  zero-padded **FFT oracle** (NumPy, or CuPy/cuFFT on NVIDIA) as a spectrally-exact ground truth.
- **Fast at scale** — a 100M-particle, 400-Myr collision runs in **~3.5 minutes** on an 8 GB
  RTX 3070 (see [Performance](#performance)).
- **Validated** — a stable Plummer sphere, a two-body Kepler orbit, multigrid-vs-oracle agreement,
  conservation diagnostics, and cross-backend determinism (CPU↔CUDA + CPU↔Metal), all in CI.

> **Status:** Stages 0–6 complete (correct CPU reference → paper reproduction → CUDA scale-up →
> Apple/Metal as a first-class backend). The realtime viewer/movie output (Stage 7) is next. See
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

Installing the package (`pip install -e .`) puts four commands on your `PATH`. What each one does:

| Command | What it does |
|---|---|
| `hello-sim` | Minimal no-op run that just confirms the toolchain and your chosen backend work. Start here. |
| `galaxy-sim` | Runs the full PM N-body pipeline from a YAML config; writes HDF5/npz snapshots + diagnostics. |
| `paper-repro` | Runs a two-galaxy collision and renders the paper's density-projection + Sun-like-tracer figures. |
| `galaxy-bench` | Times the force chain (throughput + per-stage profile). No physics output — a speed test. |

A run is configured two ways: most settings live in a **YAML config file** (passed with `--config`),
and a few can be overridden on the command line with **flags**. Both are spelled out below.

### Command-line flags

**`hello-sim`, `galaxy-sim`, `paper-repro`** take the same two core flags:

- **`--config PATH`** — the YAML file describing the run (see [Config file settings](#config-file-settings-yaml)).
  - `galaxy-sim` and `paper-repro`: **required**.
  - `hello-sim`: optional, defaults to `configs/smoke.yaml`.
- **`--backend NAME`** — which processor to run on, overriding whatever the config says. Options:
  - `cpu` — runs anywhere, no GPU needed (slowest; fine for small tests).
  - `cuda` — an NVIDIA GPU.
  - `metal` — an Apple-Silicon GPU (M-series Macs).
  - *Default:* whatever the config file's `backend` is. Every run prints the device it actually used
    and warns loudly if a requested GPU silently fell back to CPU.

**`paper-repro`** adds two figure-related flags:

- **`--tracer-radius FLOAT`** — galactocentric radius (in kpc) of the "Sun-like" tracer particle whose
  path is plotted. *Default:* `8.32` (the Sun's distance from the Galactic center).
- **`--out PATH`** — directory to write the figures into. *Default:* `figures`.

**`galaxy-bench`** is configured entirely by flags (no YAML needed):

- **`--n INT`** — number of particles. *Default:* `1000000` (1M).
- **`--grid INT`** — grid cells per side; the box is this many kpc across. *Default:* `256` (a 256³ grid).
- **`--solver NAME`** — Poisson solver. Options: `multigrid` (the portable production default) or `fft`
  (the spectrally-exact validation oracle). *Default:* `multigrid`.
- **`--backend NAME`** — `cpu` · `cuda` · `metal` (as above). *Default:* `cuda` — pass `--backend metal`
  on a Mac or `--backend cpu` with no GPU.
- **`--preset NAME`** — which initial condition to benchmark: `two_galaxy_4v` · `two_galaxy_2v` ·
  `plummer` · `hello`. *Default:* `two_galaxy_4v`.
- **`--dt FLOAT`** — timestep in Myr (millions of years). *Default:* `0.5`.
- **`--warmup INT`** — untimed steps run first (to exclude one-time GPU compilation). *Default:* `10`.
- **`--measure INT`** — timed steps the throughput number is averaged over. *Default:* `10`.

### Config file settings (YAML)

`galaxy-sim`/`paper-repro` read these from the `--config` file; ready-made configs ship in `configs/`
(`smoke.yaml` — the no-op test · `plummer.yaml` — a single equilibrium sphere · `paper_4v.yaml` /
`paper_2v.yaml` — the two collision speeds). Every setting, its options, and its default:

- **`backend`** — processor to run on: `cpu` · `cuda` · `metal`. *Default:* `cpu`.
- **`ic_preset`** — the initial condition (what's in the box at t=0):
  - `hello` — empty no-op (smoke test only).
  - `plummer` — a single, settled spherical star cluster; used to check the sim stays stable.
  - `two_galaxy_4v` — the Milky Way × Andromeda collision at the **higher** approach speed (the pair is
    less bound and spreads into a large diffuse remnant).
  - `two_galaxy_2v` — the same collision at the **lower** approach speed (the pair stays bound and compact).
  - *Default:* `hello`.
- **`solver`** — how gravity (Poisson's equation) is solved: `multigrid` (portable, runs on every
  backend — the production default) or `fft` (spectrally-exact, used to validate `multigrid`).
  *Default:* `multigrid`.
- **Particle count — set *one* of these** (or neither, to use the preset's own default):
  - **`n_particles`** — the number of particles directly (e.g. `10000000`).
  - **`particle_mass`** — solar masses (M☉) per particle; the count is derived as
    (total galaxy mass) ÷ this. Setting both is an error.
- **`grid_size`** — density-grid cells per side; the simulation box is this many kpc across, in 1 kpc³
  cells. Larger = finer resolution but more memory and compute. *Default:* `256`.
- **`dt`** — timestep in Myr. *Default:* `0.01` (the collision configs use `0.5`).
- **`steps`** — how many timesteps to integrate. *Default:* `1`.
- **`separation`** — *(two-galaxy presets only)* initial center-to-center distance between the two
  galaxies, in kpc. *Default:* `90.0`.
- **`impact_parameter`** — *(two-galaxy presets only)* sideways offset between their paths, in kpc;
  `0` is a head-on collision. *Default:* `0.0`.
- **`output_cadence`** — write a snapshot every N steps; `0` disables snapshot output entirely.
  *Default:* `0`.
- **`output_dir`** — folder snapshots are written to. *Default:* `outputs`.
- **`seed`** — random-number seed; the same seed reproduces the same run exactly. *Default:* `0`.
- **`name`** — a label for the run (used in output filenames/logs). *Default:* `hello-sim`.

## Performance

256³ grid, multigrid solver, on the two GPU backends (NVIDIA RTX 3070, 8 GB · Apple M5 Pro, 64 GB):

| N particles | CUDA (RTX 3070) | Metal (M5 Pro) |
|---|---|---|
| 1M | 33 steps/s | 19 steps/s |
| 10M | 25 steps/s | 5.6 steps/s |
| 30M | 16 steps/s | 2.1 steps/s |
| 100M | 7 steps/s (6.5 GB VRAM) | 0.65 steps/s (fits in unified memory) |

A 100M-particle, 800-step (400-Myr) collision completes in ~3.5 minutes on the 3070. Metal is
competitive where the grid solve dominates (1M) but the CIC **deposit** (global write-scatter) is its
throughput ceiling at scale — ~10× CUDA's deposit at 100M, an accepted, characterized limitation.
Full benchmark matrix, per-stage breakdown, the CUDA solve tuning, and the Metal deposit spike are in
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
GALAXY_TEST_ARCH=cuda pytest    # re-run the suite on a GPU backend (or =metal on Apple Silicon)
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
| [`docs/performance.md`](docs/performance.md) | Benchmark matrix + per-stage profile (CUDA, Stage 5; Metal, Stage 6). |
| [`docs/paper_reproduction.md`](docs/paper_reproduction.md) | Reproducing the 2020 paper's figures. |
| [`docs/stage5_plan.md`](docs/stage5_plan.md) | The CUDA scale-up plan (Stage 5). |
| [`docs/stage6_plan.md`](docs/stage6_plan.md) | The Apple/Metal plan (Stage 6). |

## Roadmap

Stages 0–6 are done: scaffold & CI → units & data model → initial conditions → correct CPU
reference → validation & FFT oracle & paper reproduction → **CUDA scale-up to 100M** → **Apple/Metal
as a first-class backend** (Kahan-fp32 reductions, 3-way parity, benchmark; the CIC deposit is the
Metal throughput ceiling at scale — see [Performance](#performance)). Up next: **Stage 7** (realtime
Taichi GGUI viewer + batch→movie output) and **Stage 8** (research campaigns). See
[`AGENT.md` §7](AGENT.md) for per-stage detail and exit criteria.

## License & origin

Apache-2.0 — see [`LICENSE`](LICENSE). Originally *"Galaxy Collisions With CUDA and FFT"* (Adams,
Lajevardi, LeBlanc, Movaseghi; PHYD57, 2020); the original CUDA source is preserved under
[`legacy/`](legacy/).
