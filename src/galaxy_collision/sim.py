"""Simulation orchestration and CLI entry points.

Stage 0 ships only ``hello-sim``: a trivial, no-op Taichi step loop that proves the
package imports, a config loads, the chosen backend initializes, and a kernel runs.
It is the green baseline every later stage must keep passing — not real physics.

The real ``kick -> drift -> deposit -> solve -> grad -> kick`` pipeline (AGENT.md
§5.1) replaces the no-op step starting at Stage 3.
"""

import argparse
from pathlib import Path
from typing import Any

from galaxy_collision.config import SimConfig, load_config

# Map our backend strings to Taichi archs. Imported lazily inside functions so that
# config-only code (and CI lint) does not require the Taichi runtime.
_ARCH_NAMES = {"cpu": "cpu", "cuda": "cuda", "metal": "metal"}


def run_hello_sim(config: SimConfig) -> dict[str, Any]:
    """Run the no-op smoke simulation and return a small summary dict.

    Initializes Taichi on the configured backend, allocates a tiny SoA position
    field, and advances it ``config.steps`` times with a kernel that does nothing
    physical. Verifies the whole toolchain end-to-end without committing to any
    physics decisions.
    """
    import taichi as ti

    arch = getattr(ti, _ARCH_NAMES[config.backend])
    ti.init(arch=arch, random_seed=config.seed, offline_cache=False)

    # A handful of placeholder particles — enough to exercise a kernel, cheap enough
    # to run anywhere (CI included).
    n = min(config.resolve_n_particles() or 16, 1024)
    pos = ti.Vector.field(3, dtype=ti.f32, shape=n)

    @ti.kernel
    def seed_positions():
        for i in pos:
            pos[i] = ti.Vector([ti.f32(i), 0.0, 0.0])

    @ti.kernel
    def noop_step(dt: ti.f32):
        # Intentionally a no-op drift: add zero velocity * dt. Real integrator lands
        # in Stage 3. The multiply-by-dt keeps dt in the kernel signature so the
        # smoke test exercises a parameterized launch.
        for i in pos:
            pos[i] += ti.Vector([0.0, 0.0, 0.0]) * dt

    seed_positions()
    for _ in range(config.steps):
        noop_step(config.dt)
    ti.sync()

    return {
        "name": config.name,
        "backend": config.backend,
        "arch": _ARCH_NAMES[config.backend],
        "n_particles": n,
        "steps": config.steps,
        "dt": config.dt,
        "status": "ok",
    }


def hello_sim_cli(argv: list[str] | None = None) -> int:
    """CLI entry point: ``hello-sim [--config PATH] [--backend NAME]``."""
    parser = argparse.ArgumentParser(
        prog="hello-sim",
        description="Run the Stage-0 no-op smoke simulation.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "configs" / "smoke.yaml",
        help="Path to a YAML run config (default: configs/smoke.yaml).",
    )
    parser.add_argument(
        "--backend",
        choices=tuple(_ARCH_NAMES),
        default=None,
        help="Override the backend from the config.",
    )
    args = parser.parse_args(argv)

    config = load_config(args.config)
    if args.backend is not None:
        config.backend = args.backend
        config.validate()

    result = run_hello_sim(config)
    print(
        f"hello-sim ok: name={result['name']} backend={result['backend']} "
        f"n={result['n_particles']} steps={result['steps']} dt={result['dt']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(hello_sim_cli())
