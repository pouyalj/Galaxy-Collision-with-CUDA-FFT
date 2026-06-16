# Galaxy Collision

A Particle-Mesh (PM) gravitational N-body simulation of a Milky Way × Andromeda
galaxy collision. Two spiral galaxies are started apart, pushed toward each other,
and evolved through the collision on a 256³ density grid.

Originally a 2020 University of Toronto Scarborough PHYD57 project
([paper PDF](https://static1.squarespace.com/static/5d5dcd310b9b0100013bcbe1/t/5f05fb53b5c64717ff32882e/1594227543673/PHYD57_NBody-Project_No1-2020.pdf)),
now being rebuilt as **one portable source** that runs on CPU, NVIDIA CUDA, and
Apple-Silicon GPU via [Taichi](https://www.taichi-lang.org/), with research-grade
physics at 10–100M particles.

> **Status: Stage 0 (scaffold).** The package, config loader, CI, and a trivial
> `hello-sim` smoke run are in place. Real physics arrives at Stage 3. See
> [`AGENT.md`](AGENT.md) for the full architecture and staged plan, and
> [`docs/development.md`](docs/development.md) to get started.

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
