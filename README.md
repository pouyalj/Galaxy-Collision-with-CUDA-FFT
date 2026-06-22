# Galaxy Collision

A Particle-Mesh (PM) gravitational N-body simulation of a Milky Way × Andromeda
galaxy collision. Two spiral galaxies are started apart, pushed toward each other,
and evolved through the collision on a 256³ density grid.

Originally a 2020 University of Toronto Scarborough PHYD57 project
([paper PDF](https://static1.squarespace.com/static/5d5dcd310b9b0100013bcbe1/t/5f05fb53b5c64717ff32882e/1594227543673/PHYD57_NBody-Project_No1-2020.pdf)),
now being rebuilt as **one portable source** that runs on CPU, NVIDIA CUDA, and
Apple-Silicon GPU via [Taichi](https://www.taichi-lang.org/), with research-grade
physics at 10–100M particles.

> **Status: Stages 0–3 done; Stage 4 in progress.** On top of the scaffold, config loader,
> CI, units, SoA data model, and two-galaxy/Plummer initial conditions, there is a **correct
> CPU simulation**: CIC deposit/gather, both Poisson solvers (open-boundary multigrid plus a
> zero-padded FFT oracle), a KDK leapfrog with Plummer softening, fp64 conservation
> diagnostics, and HDF5/npz snapshots — driven by `galaxy-sim`. It is validated by a stable
> Plummer sphere, a two-body Kepler orbit, and multigrid-vs-oracle agreement. **Stage 4**
> (validation hardening + paper reproduction) is underway — Checkpoints **4A** (principled
> multigrid tolerances, disks launched in equilibrium with the grid's own measured rotation
> curve, a big-run energy metric, a collision smoke test) and **4B** (a Sun-like tracer plus
> the density-projection-sequence and tracer-path figures via the `paper-repro` CLI) are in;
> the production 4v/2v runs are Checkpoint 4C. See [`AGENT.md`](AGENT.md) for the full
> architecture and staged plan, and [`docs/development.md`](docs/development.md) to get started.

## Quickstart

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
hello-sim                                  # no-op smoke simulation on CPU
galaxy-sim --config configs/plummer.yaml   # real PM N-body run (a Plummer sphere)
paper-repro --config configs/paper_4v.yaml # two-galaxy collision -> paper figures (heavy; see config)
pytest                                     # run the test suite
```

The original 2020 CUDA source is preserved under [`legacy/`](legacy/).

## License

Apache-2.0 — see [`LICENSE`](LICENSE).

## Things to revisit

A code review of Stages 0–1 flagged a few small, non-urgent items — **most are now fixed**. Full
detail and status are in [`AGENT.md` §11](AGENT.md).

Resolved (2026-06-16):

- ✅ The plain-math tests no longer skip when the GPU library (Taichi) is absent — they always run.
- ✅ The docs' memory note now matches the code: each particle stores its own mass (needed for the
  galaxies' central black holes), so 100M particles is ~3.2 GB, not 2.8 GB.
- ✅ Install requires **Python 3.11 or 3.12** (the GPU library has no 3.13 wheels yet); this is now
  enforced by the package metadata.
- ✅ Input checks are consistent (counts must be positive everywhere).

Resolved (2026-06-20, Stage 3):

- ✅ **Grid-memory accounting updated** for the real solver buffers: the estimate now counts ρ, Φ,
  the acceleration field, and the multigrid hierarchy (~0.57 GB at 256³), not just ρ/Φ.

Still open (deferred to Stage 5, performance):

- The open-boundary multigrid is *correct* but not yet *fast* (~0.5–0.7 convergence/cycle) — to be
  tuned when chasing GPU performance. The boundary-condition moments are also computed on the host
  for now; moving them onto the device is part of the same Stage-5 performance work.
