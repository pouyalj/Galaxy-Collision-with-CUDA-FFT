# Galaxy Collision

A Particle-Mesh (PM) gravitational N-body simulation of a Milky Way × Andromeda
galaxy collision. Two spiral galaxies are started apart, pushed toward each other,
and evolved through the collision on a 256³ density grid.

Originally a 2020 University of Toronto Scarborough PHYD57 project
([paper PDF](https://static1.squarespace.com/static/5d5dcd310b9b0100013bcbe1/t/5f05fb53b5c64717ff32882e/1594227543673/PHYD57_NBody-Project_No1-2020.pdf)),
now being rebuilt as **one portable source** that runs on CPU, NVIDIA CUDA, and
Apple-Silicon GPU via [Taichi](https://www.taichi-lang.org/), with research-grade
physics at 10–100M particles.

> **Status: Stages 0–2 done.** Scaffold, config loader, CI, and a trivial
> `hello-sim` are in place, plus a tested (kpc, Myr, M☉) unit system, the SoA
> particle/grid data model, and a reproducible two-galaxy initial-condition
> generator (disk + bulge + central black hole, 4v/2v approach). Next is Stage 3
> (the correct CPU simulation). See [`AGENT.md`](AGENT.md) for the full architecture
> and staged plan, and [`docs/development.md`](docs/development.md) to get started.

## Quickstart

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
hello-sim          # runs the no-op smoke simulation on CPU
pytest             # run the test suite
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

Still open:

- **Grid-memory accounting will need updating in later stages** as the Poisson solver adds its own
  grid buffers (the multigrid hierarchy in Stage 3, the FFT buffer in Stage 4).
