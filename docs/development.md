# Development

The modernized simulator is a Python package (`galaxy_collision`) using
[Taichi](https://www.taichi-lang.org/) kernels that compile to CPU / CUDA / Apple
Metal from one source. See [`AGENT.md`](../AGENT.md) for the full architecture and
the staged implementation plan; this file is just the quickstart.

## Setup

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

The optional `cuda` extra (`pip install -e ".[dev,cuda]"`) pulls in CuPy for the
Stage-4 zero-padded FFT validation oracle — NVIDIA only; not needed otherwise.

## Run the smoke simulation

```bash
hello-sim                      # uses configs/smoke.yaml on the CPU backend
hello-sim --backend cuda       # if an NVIDIA GPU + CUDA runtime are present
hello-sim --config path.yaml   # any YAML matching the SimConfig schema
```

`hello-sim` is a Stage-0 no-op step loop — it proves the package, config loader,
and Taichi backend all work. Real physics arrives at Stage 3.

## Checks (what CI runs)

```bash
ruff check .   # lint
pytest         # config + smoke tests (the Taichi test self-skips without Taichi)
hello-sim --backend cpu
```

### Re-running the suite on a GPU backend

The physics/solver tests pin Taichi to CPU for portable CI. To re-pass the same suite on a GPU
(AGENT.md §7 requires this for Stages 5/6), set `GALAXY_TEST_ARCH` — a `conftest.py` fixture
redirects the per-test `ti.init` onto that backend:

```bash
GALAXY_TEST_ARCH=cuda pytest -q     # or: metal
```

Tests whose purpose is a specific backend's reporting/fallback semantics, or that compare two
backends, opt out with `@pytest.mark.fixed_arch`.

## Layout

```
src/galaxy_collision/
  config.py        run-config schema + YAML loader (Stage 0)
  sim.py           orchestration + hello-sim CLI (Stage 0)
  units.py         (kpc, Myr, M_sun) system + derived G (Stage 1)
  data.py          SoA particle/grid fields + memory estimator (Stage 1)
  ic.py            two-galaxy initial conditions (Stage 2)
  solver/          Poisson solvers — multigrid + FFT oracle (Stage 3+)
  viz/             paper_repro static figures (Stage 4/4B); realtime GGUI + movie (Stage 7)
configs/           YAML run configs (smoke.yaml is the baseline)
tests/             pytest suite
legacy/            the original 2020 CUDA .cu files, preserved for reference
```
