# Galaxy Collision

A Particle-Mesh (PM) gravitational N-body simulation of a Milky Way × Andromeda
galaxy collision. Two spiral galaxies are started apart, pushed toward each other,
and evolved through the collision on a 256³ density grid.

Originally a 2020 University of Toronto Scarborough PHYD57 project
([paper PDF](https://static1.squarespace.com/static/5d5dcd310b9b0100013bcbe1/t/5f05fb53b5c64717ff32882e/1594227543673/PHYD57_NBody-Project_No1-2020.pdf)),
now being rebuilt as **one portable source** that runs on CPU, NVIDIA CUDA, and
Apple-Silicon GPU via [Taichi](https://www.taichi-lang.org/), with research-grade
physics at 10–100M particles.

> **Status: Stages 0–1 done.** Scaffold, config loader, CI, and a trivial
> `hello-sim` are in place, plus a tested (kpc, Myr, M☉) unit system and the SoA
> particle/grid data model. Next is Stage 2 (initial conditions); real physics
> arrives at Stage 3. See [`AGENT.md`](AGENT.md) for the full architecture and
> staged plan, and [`docs/development.md`](docs/development.md) to get started.

## Quickstart

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
hello-sim          # runs the no-op smoke simulation on CPU
pytest             # run the test suite
```

The original 2020 CUDA source is preserved under [`legacy/`](legacy/).

## License

MIT — see [`LICENSE`](LICENSE).

## Things to revisit

A code review of Stages 0–1 flagged a few small, non-urgent items. **Nothing here breaks the current
build** — they're little improvements and decisions to make later. Full detail is in
[`AGENT.md` §11](AGENT.md).

- **Some tests are skipped when the GPU library (Taichi) isn't installed — including a few that don't
  actually need it.** A handful of plain-math tests share a file with the GPU tests and get skipped
  together, so they don't run on a machine without Taichi. Worth rearranging so they always run.
- **Each particle now stores its own mass.** That's needed for the galaxies' central black holes, but
  it makes a particle a bit bigger than the original plan assumed — so the memory note in the docs
  should be bumped up to match (about 3.2 GB for 100M particles, not 2.8 GB).
- **Use Python 3.11 or 3.12 for now.** The GPU library doesn't support the newest Python (3.13) yet,
  so installing on 3.13 can fail.
- **A couple of tiny cleanups.** Make the input checks consistent, and remember to update the memory
  estimate when more grid data is added in later stages.
