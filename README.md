# Galaxy Collision

A Particle-Mesh (PM) gravitational N-body simulation of a Milky Way × Andromeda
galaxy collision. Two spiral galaxies are started apart, pushed toward each other,
and evolved through the collision on a 256³ density grid.

Originally a 2020 University of Toronto Scarborough PHYD57 project
([paper PDF](https://static1.squarespace.com/static/5d5dcd310b9b0100013bcbe1/t/5f05fb53b5c64717ff32882e/1594227543673/PHYD57_NBody-Project_No1-2020.pdf)),
now being rebuilt as **one portable source** that runs on CPU, NVIDIA CUDA, and
Apple-Silicon GPU via [Taichi](https://www.taichi-lang.org/), with research-grade
physics at 10–100M particles.

> **Status: Stages 0–5 done.** On top of the scaffold, config loader, CI, units, SoA data
> model, and two-galaxy/Plummer initial conditions, there is a **correct CPU simulation**: CIC
> deposit/gather, both Poisson solvers (open-boundary multigrid plus a zero-padded FFT oracle),
> a KDK leapfrog with Plummer softening, fp64 conservation diagnostics, and HDF5/npz snapshots
> — driven by `galaxy-sim`. It is validated by a stable Plummer sphere, a two-body Kepler orbit,
> and multigrid-vs-oracle agreement. **Stage 4** hardened that engine (principled solver
> tolerances; disks launched in equilibrium with the grid's own measured rotation curve) and
> **qualitatively reproduced the paper's headline figures** — the 4v & 2v Milky-Way × Andromeda
> collisions at 256³ / 10M particles plus the Sun-like tracer path, via the `paper-repro` CLI
> ([`docs/paper_reproduction.md`](docs/paper_reproduction.md); the gross outcome and 4v↔2v contrast
> match, with significant cold-disk numerical heating at this resolution — a warm disk is a later
> refinement). **Stage 5 (CUDA scale-up) is done** on an RTX 3070: the force chain runs
> device-resident on the GPU, the multigrid solve was tuned for a **6.5–12.5× throughput gain**
> (adaptive warm-start cycling + a thread-local reduction fix), the FFT oracle gained a CuPy/cuFFT
> GPU path, and a **100M-particle collision** ran on the 8 GB card (400 Myr in ~3.5 min, 6.5 GB
> peak) — see [`docs/performance.md`](docs/performance.md). Suite green on both CPU and CUDA.
> Next: Stage 6 (Apple/Metal). See [`AGENT.md`](AGENT.md) for the architecture and staged plan,
> [`docs/stage5_plan.md`](docs/stage5_plan.md) for the Stage-5 plan, and
> [`docs/development.md`](docs/development.md) to get started.

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

Resolved (2026-06-25, Stage 5):

- ✅ **The multigrid is now fast.** Rather than tune the per-cycle rate, Stage 5 exploits that the
  time loop warm-starts every solve — one V-cycle already hits the floor — so the solver
  early-terminates adaptively (cold start keeps the full cap). Combined with a thread-local
  reduction fix, the solve is ~13× faster and the whole step **6.5–12.5×** faster
  ([`docs/performance.md`](docs/performance.md)).
- ✅ **The boundary-condition moments now run on the device** (reduced on-GPU each solve; only 10
  scalars cross to the host), along with the conservation diagnostics — no per-step full-grid
  transfers.
