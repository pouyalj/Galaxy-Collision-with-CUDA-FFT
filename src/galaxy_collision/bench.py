"""Per-stage profiler + throughput benchmark for the PM force chain (Stage 5 / 5B).

AGENT.md §5.6: the Stage-5 exit gate wants steps/sec at N ∈ {1M, 10M, 30M, 100M} on CUDA
*plus* a per-stage breakdown (deposit / solve / grad / gather / integrate). This module is
that instrument. It rebuilds the same device-resident force chain ``run_simulation`` runs —
``deposit → solve → potential_to_accel → gather → KDK`` — but strips snapshots / diagnostics /
history so the numbers are pure kernel cost.

Two measurement passes, deliberately separate:

* **Throughput** — ``measure`` steps timed wall-clock with a *single* ``ti.sync()`` at the end,
  so async kernel launches overlap exactly as in a real run. This is the steps/sec number.
* **Per-stage** — the same step decomposed, each stage wrapped in ``ti.sync()`` so its device
  time is attributed correctly. The syncs serialize the pipeline, so this pass is a bit slower
  than throughput — read the per-stage table for *relative* cost, not absolute steps/sec.

A few warmup steps run first so Taichi's one-time JIT (heavy for the 256³ multigrid kernels)
and the solver warm-start are excluded from both measurements.
"""

from __future__ import annotations

import argparse
import time
from collections.abc import Callable
from typing import Any

from galaxy_collision.config import SimConfig


def _gpu_mem_used_mib() -> float | None:
    """Best-effort current GPU memory use (MiB) via nvidia-smi; None if unavailable."""
    import subprocess

    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10, check=True,
        )
        return float(out.stdout.strip().splitlines()[0])
    except Exception:
        return None


def _timed(ti, fn: Callable[[], None]) -> float:
    """Run ``fn`` and return its wall time, bracketed by ``ti.sync()`` for true device time."""
    ti.sync()
    t0 = time.perf_counter()
    fn()
    ti.sync()
    return time.perf_counter() - t0


# The stage order of one KDK step. Both half-kicks + the drift are attributed to "integrate".
_STAGES = ("integrate", "deposit", "solve", "grad", "gather")


def run_benchmark(
    config: SimConfig,
    *,
    warmup: int = 3,
    measure: int = 10,
    deposit_fn: Callable[..., None] | None = None,
) -> dict[str, Any]:
    """Benchmark the force chain for ``config``; return throughput + per-stage breakdown.

    ``deposit_fn`` overrides the deposit step (to compare deposition variants); it is called
    as ``deposit_fn(parts, rho, dx)`` like :func:`deposit.deposit_density` (the default).
    """
    import taichi as ti

    from galaxy_collision import ic as ic_mod
    from galaxy_collision import units
    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.deposit import deposit_density, gather_acceleration, potential_to_accel
    from galaxy_collision.integrator import drift, kick
    from galaxy_collision.sim import DX, _build_ic, init_backend
    from galaxy_collision.solver import make_solver

    backend = init_backend(config)
    deposit_fn = deposit_fn or deposit_density

    icr = _build_ic(config)
    n, gs = icr.n, config.grid_size
    parts = ParticleState(n)
    ic_mod.load_into_particle_state(icr, parts)
    grid = GridState(gs)
    acc_x = ti.field(ti.f32, shape=n)
    acc_y = ti.field(ti.f32, shape=n)
    acc_z = ti.field(ti.f32, shape=n)
    solver = make_solver(config.solver, gs, dx=DX, grav_constant=units.G)
    half = 0.5 * config.dt
    warm = {"on": False}

    # Match run_simulation: launch the disks in equilibrium so the warm-start quality (and thus
    # the adaptive solver's per-step cycle count) is representative of a real run, not an
    # artificially hard transient from cold/wrong velocities.
    if config.ic_preset in ("two_galaxy_4v", "two_galaxy_2v"):
        ic_mod.equilibrate_disk_velocities(
            icr, dx=DX, solver=solver, rho=grid.rho, phi=grid.phi,
            ax=grid.ax, ay=grid.ay, az=grid.az, parts=parts,
        )

    def do_deposit() -> None:
        deposit_fn(parts, grid.rho, DX)

    def do_solve() -> None:
        solver.solve(grid.rho, grid.phi, warm_start=warm["on"])
        warm["on"] = True

    def do_grad() -> None:
        potential_to_accel(grid.phi, grid.ax, grid.ay, grid.az, DX)

    def do_gather() -> None:
        gather_acceleration(parts, grid.ax, grid.ay, grid.az, acc_x, acc_y, acc_z, DX)

    def do_kick() -> None:
        kick(parts.vel_x, parts.vel_y, parts.vel_z, acc_x, acc_y, acc_z, n, half)

    def do_drift() -> None:
        drift(parts.pos_x, parts.pos_y, parts.pos_z,
              parts.vel_x, parts.vel_y, parts.vel_z, n, config.dt)

    def step() -> None:
        do_kick()
        do_drift()
        do_deposit()
        do_solve()
        do_grad()
        do_gather()
        do_kick()

    # Prime the cached accel, then warm up (excludes JIT + solver cold-start).
    do_deposit()
    do_solve()
    do_grad()
    do_gather()
    for _ in range(warmup):
        step()
    ti.sync()
    mem_mib = _gpu_mem_used_mib() if config.backend == "cuda" else None

    # Throughput: one sync at the end so launches overlap like a real run. Record the adaptive
    # solver's per-step cycle count (if the solver exposes it) to explain the solve cost.
    cycles: list[int] = []
    t0 = time.perf_counter()
    for _ in range(measure):
        step()
        if hasattr(solver, "last_cycles"):
            cycles.append(solver.last_cycles)
    ti.sync()
    wall = time.perf_counter() - t0
    mean_cycles = sum(cycles) / len(cycles) if cycles else None

    # Per-stage: each stage individually synced (serializes — relative cost, not throughput).
    per_stage = {s: 0.0 for s in _STAGES}
    for _ in range(measure):
        per_stage["integrate"] += _timed(ti, do_kick)
        per_stage["integrate"] += _timed(ti, do_drift)
        per_stage["deposit"] += _timed(ti, do_deposit)
        per_stage["solve"] += _timed(ti, do_solve)
        per_stage["grad"] += _timed(ti, do_grad)
        per_stage["gather"] += _timed(ti, do_gather)
        per_stage["integrate"] += _timed(ti, do_kick)
    per_stage_ms = {s: per_stage[s] / measure * 1e3 for s in _STAGES}

    return {
        "device": backend["device"],
        "backend_resolved": backend["resolved"],
        "n_particles": n,
        "grid_size": gs,
        "solver": config.solver,
        "steps_per_sec": measure / wall,
        "ms_per_step": wall / measure * 1e3,
        "per_stage_ms": per_stage_ms,
        "gpu_mem_used_mib": mem_mib,
        "mean_solver_cycles": mean_cycles,
    }


def _format_report(r: dict[str, Any]) -> str:
    lines = [
        f"device      : {r['device']}  (arch={r['backend_resolved']})",
        f"N particles : {r['n_particles']:,}",
        f"grid        : {r['grid_size']}^3   solver={r['solver']}",
        f"throughput  : {r['steps_per_sec']:.2f} steps/s   ({r['ms_per_step']:.1f} ms/step)",
    ]
    if r["gpu_mem_used_mib"] is not None:
        lines.append(f"GPU mem     : {r['gpu_mem_used_mib']:.0f} MiB used")
    if r["mean_solver_cycles"] is not None:
        lines.append(f"solver      : {r['mean_solver_cycles']:.1f} V-cycles/step (mean, adaptive)")
    lines.append("per-stage (ms/step, individually synced):")
    total = sum(r["per_stage_ms"].values())
    for s in _STAGES:
        ms = r["per_stage_ms"][s]
        pct = 100.0 * ms / total if total > 0 else 0.0
        lines.append(f"    {s:<10s} {ms:8.2f} ms   {pct:5.1f}%")
    lines.append(f"    {'(sum)':<10s} {total:8.2f} ms")
    return "\n".join(lines)


def bench_cli(argv: list[str] | None = None) -> int:
    """CLI: ``galaxy-bench --n 10000000 [--grid 256] [--solver multigrid] [--backend cuda]``."""
    p = argparse.ArgumentParser(prog="galaxy-bench", description="Benchmark the PM force chain.")
    p.add_argument("--n", type=int, default=1_000_000, help="particle count")
    p.add_argument("--grid", type=int, default=256, help="grid size per side")
    p.add_argument("--solver", default="multigrid", help="multigrid | fft")
    p.add_argument("--backend", default="cuda", help="cpu | cuda | metal")
    p.add_argument("--preset", default="two_galaxy_4v", help="IC preset")
    p.add_argument("--dt", type=float, default=0.5)
    p.add_argument("--warmup", type=int, default=10)
    p.add_argument("--measure", type=int, default=10)
    args = p.parse_args(argv)

    config = SimConfig(
        name="bench", backend=args.backend, ic_preset=args.preset, n_particles=args.n,
        dt=args.dt, steps=0, solver=args.solver, grid_size=args.grid, output_cadence=0, seed=0,
    )
    r = run_benchmark(config, warmup=args.warmup, measure=args.measure)
    print(_format_report(r))
    return 0


if __name__ == "__main__":
    raise SystemExit(bench_cli())
