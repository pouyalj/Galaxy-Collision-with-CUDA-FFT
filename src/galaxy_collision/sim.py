"""Simulation orchestration and CLI entry points.

``hello-sim`` (Stage 0) is a trivial, no-op Taichi step loop that proves the package
imports, a config loads, the backend initializes, and a kernel runs — the green smoke
baseline every later stage must keep passing.

``run_simulation`` (Stage 3) is the real thing: the device-resident
**kick → drift → deposit → solve → grad → gather → kick** KDK pipeline (AGENT.md §5.1).
The main loop owns acceleration *priming* — it evaluates the force chain once before the
loop and then trusts :func:`~galaxy_collision.integrator.kdk_step` to recompute it exactly
once per step (the single point where the cached accel and the positions are reconciled,
per that function's contract). Conservation diagnostics are recorded along the way and
snapshots written at ``output_cadence``.
"""

import argparse
import platform
import subprocess
from pathlib import Path
from typing import Any

from galaxy_collision.config import SimConfig, load_config

# Map our backend strings to Taichi archs. Imported lazily inside functions so that
# config-only code (and CI lint) does not require the Taichi runtime.
_ARCH_NAMES = {"cpu": "cpu", "cuda": "cuda", "metal": "metal"}

# Cell size: the box is ``grid_size`` kpc with 1 kpc^3 cells (AGENT.md §3.4), node-centered.
DX = 1.0

# Taichi reports the resolved CPU arch as its native name, not "cpu".
_CPU_ARCHS = frozenset({"arm64", "x64", "cpu"})


def _backend_honored(requested: str, resolved: str) -> bool:
    """True if Taichi initialized on the requested backend (no silent fallback).

    Taichi falls back to CPU when a GPU backend is unavailable (e.g. ``cuda`` on a Mac),
    reporting the native CPU arch (``arm64``/``x64``). A ``cpu`` request is honored by any
    CPU arch; ``cuda``/``metal`` must resolve to themselves.
    """
    if requested == "cpu":
        return resolved in _CPU_ARCHS
    return requested == resolved


def _hardware_description(resolved_arch: str) -> str:
    """Human-readable description of the device a resolved Taichi arch runs on."""
    cpu = platform.processor() or platform.machine() or "CPU"
    if platform.system() == "Darwin":
        try:
            out = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True, text=True, timeout=2,
            )
            cpu = out.stdout.strip() or cpu
        except (OSError, subprocess.SubprocessError):
            pass
    if resolved_arch == "cuda":
        try:
            out = subprocess.run(
                ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                capture_output=True, text=True, timeout=5,
            )
            names = [ln.strip() for ln in out.stdout.splitlines() if ln.strip()]
            if names:
                extra = f" (+{len(names) - 1} more)" if len(names) > 1 else ""
                return f"NVIDIA {names[0]}{extra}"
        except (OSError, subprocess.SubprocessError):
            pass
        return "NVIDIA CUDA GPU"
    if resolved_arch == "metal":
        return f"{cpu} — Apple GPU (Metal)"
    if resolved_arch == "vulkan":
        return f"{cpu} — Vulkan GPU"
    return f"{cpu} (CPU)"


def init_backend(config: SimConfig) -> dict[str, Any]:
    """Initialize Taichi on the configured backend and report the detected hardware.

    Prints one line announcing the device the simulation will actually run on. If the
    requested backend is unavailable, Taichi silently falls back (e.g. ``cuda`` → CPU on a
    Mac); we detect that and print a clear WARNING instead, so a run never *looks* like it
    used a GPU it didn't. Returns ``{requested, resolved, device, fell_back}``.
    """
    import taichi as ti

    arch = getattr(ti, _ARCH_NAMES[config.backend])
    ti.init(arch=arch, random_seed=config.seed, offline_cache=False)
    try:
        resolved = ti.lang.impl.current_cfg().arch.name
    except Exception:  # pragma: no cover - defensive against Taichi internals moving
        resolved = "unknown"
    info = {
        "requested": config.backend,
        "resolved": resolved,
        "device": _hardware_description(resolved),
        "fell_back": not _backend_honored(config.backend, resolved),
    }
    if info["fell_back"]:
        print(
            f"WARNING: backend '{info['requested']}' is unavailable here — Taichi fell back to "
            f"'{resolved}'. Running on {info['device']}."
        )
    else:
        print(f"galaxy-collision: running on {info['device']}  "
              f"(backend={info['requested']}, arch={resolved})")
    return info


def run_hello_sim(config: SimConfig) -> dict[str, Any]:
    """Run the no-op smoke simulation and return a small summary dict.

    Initializes Taichi on the configured backend, allocates a tiny SoA position
    field, and advances it ``config.steps`` times with a kernel that does nothing
    physical. Verifies the whole toolchain end-to-end without committing to any
    physics decisions.
    """
    backend = init_backend(config)
    import taichi as ti

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
        "backend_resolved": backend["resolved"],
        "device": backend["device"],
        "arch": _ARCH_NAMES[config.backend],
        "n_particles": n,
        "steps": config.steps,
        "dt": config.dt,
        "status": "ok",
    }


def _build_ic(config: SimConfig, plummer_model=None):
    """Dispatch to the right IC builder for the configured preset (Stage 3 presets only)."""
    from galaxy_collision import ic as ic_mod

    if config.ic_preset == "plummer":
        return ic_mod.build_plummer_ic(config, model=plummer_model)
    if config.ic_preset in ("two_galaxy_4v", "two_galaxy_2v"):
        return ic_mod.build_ic(config)
    raise NotImplementedError(
        f"ic_preset {config.ic_preset!r} has no Stage-3 simulation builder "
        f"(available: plummer, two_galaxy_4v, two_galaxy_2v)"
    )


def _check_box_fit(icr, grid_size: int, tol: float = 0.10) -> None:
    """Fail loudly if too much mass falls outside the [0, grid_size−1] box.

    The CIC deposit silently drops out-of-box weight (open BCs, legacy bug #2 guard), so a
    badly-placed IC (e.g. a two-galaxy run on too small a grid — centers can land outside the
    box, dropping ~half the mass) would otherwise corrupt the density, forces, and energy
    while still reporting 'ok'. Better to refuse than lie. The 10% threshold catches gross
    misplacement while tolerating the intrinsic fat tail of a Plummer sphere (a few % of its
    mass legitimately sits beyond any modest box).
    """
    import numpy as np

    outside = np.any((icr.pos < 0.0) | (icr.pos > grid_size - 1), axis=1)
    frac = float(icr.mass[outside].sum() / icr.mass.sum())
    if frac > tol:
        raise ValueError(
            f"{frac:.1%} of the IC mass lies outside the {grid_size}^3 box — the initial "
            f"conditions do not fit (check grid_size vs galaxy placement/separation)."
        )


def run_simulation(
    config: SimConfig,
    *,
    history_cadence: int | None = None,
    write_snapshots: bool = True,
    solver_kwargs: dict[str, Any] | None = None,
    plummer_model=None,
    direct_pe_softening: float | None = None,
) -> dict[str, Any]:
    """Run the full PM N-body pipeline and return a summary (incl. a diagnostics history).

    ``history_cadence`` controls how often conservation diagnostics are recorded (default:
    ~20 samples across the run). Snapshots are written to ``config.output_dir`` every
    ``config.output_cadence`` steps (plus a final one) when ``write_snapshots`` and that
    cadence are set. ``solver_kwargs`` is forwarded to the Poisson solver (e.g. ``n_cycles``
    for multigrid). ``plummer_model`` overrides the default :class:`PlummerModel` for the
    ``plummer`` preset (lets callers pick a denser, faster-evolving sphere for tests).

    ``direct_pe_softening`` (kpc), when set, also records ``energy_direct`` in the history
    using the exact O(N²) softened pair-sum PE. This is the *clean* conservation metric — it
    omits the grid PE's fluctuating CIC self-energy term — but is O(N²), so enable it only
    for small N (tests/validation), never for production-scale runs. Leave ``None`` otherwise.
    """
    import numpy as np
    import taichi as ti

    from galaxy_collision import diagnostics, units
    from galaxy_collision import ic as ic_mod
    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.deposit import deposit_density, gather_acceleration, potential_to_accel
    from galaxy_collision.integrator import kdk_step
    from galaxy_collision.io import snapshot_from_states, write_snapshot
    from galaxy_collision.solver import make_solver

    backend = init_backend(config)

    icr = _build_ic(config, plummer_model=plummer_model)
    n, gsize = icr.n, config.grid_size
    _check_box_fit(icr, gsize)
    parts = ParticleState(n)
    ic_mod.load_into_particle_state(icr, parts)
    grid = GridState(gsize)

    acc_x = ti.field(ti.f32, shape=n)
    acc_y = ti.field(ti.f32, shape=n)
    acc_z = ti.field(ti.f32, shape=n)
    ax_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    ay_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    az_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))

    solver = make_solver(
        config.solver, gsize, dx=DX, grav_constant=units.G, **(solver_kwargs or {})
    )

    # The force chain. After the first (cold) call, warm-start the iterative solver from
    # the previous step's potential — Φ moves little per step, so a few cycles suffice.
    warm = {"on": False}

    def accel_fn():
        deposit_density(parts, grid.rho, DX)
        solver.solve(grid.rho, grid.phi, warm_start=warm["on"])
        potential_to_accel(grid.phi, ax_g, ay_g, az_g, DX)
        gather_acceleration(parts, ax_g, ay_g, az_g, acc_x, acc_y, acc_z, DX)
        warm["on"] = True

    accel_fn()  # prime: the main loop owns the single accel/positions reconciliation point

    mass_np = parts.mass.to_numpy().astype(np.float64)

    def measure(step: int) -> dict[str, Any]:
        pos = np.stack(
            [parts.pos_x.to_numpy(), parts.pos_y.to_numpy(), parts.pos_z.to_numpy()], axis=1
        ).astype(np.float64)
        vel = np.stack(
            [parts.vel_x.to_numpy(), parts.vel_y.to_numpy(), parts.vel_z.to_numpy()], axis=1
        ).astype(np.float64)
        ke = diagnostics.kinetic_energy(mass_np, vel)
        pe = diagnostics.potential_energy_grid(grid.rho.to_numpy(), grid.phi.to_numpy(), DX)
        rec = {
            "step": step,
            "time": step * config.dt,
            "energy": ke + pe,
            "kinetic": ke,
            "potential": pe,
            "momentum": diagnostics.linear_momentum(mass_np, vel),
            "ang_momentum": diagnostics.angular_momentum(mass_np, pos, vel),
            "half_mass_radius": diagnostics.lagrangian_radius(mass_np, pos, 0.5),
        }
        if direct_pe_softening is not None:
            pe_d = diagnostics.potential_energy_direct(
                mass_np, pos, softening=direct_pe_softening, grav=units.G
            )
            rec["potential_direct"] = pe_d
            rec["energy_direct"] = ke + pe_d
        return rec

    hist_cad = history_cadence or max(1, config.steps // 20)
    snap_dir = Path(config.output_dir)
    snapshots: list[str] = []
    # Scalar/vector diagnostics worth persisting alongside each snapshot.
    _diag_keys = ("energy", "kinetic", "potential", "momentum", "ang_momentum", "half_mass_radius")

    def _snapshot_due(step: int) -> bool:
        if not (write_snapshots and config.output_cadence > 0):
            return False
        # Every output_cadence steps, plus always the final state (even if off-cadence).
        return step % config.output_cadence == 0 or step == config.steps

    def write_snap(step: int, m: dict[str, Any]) -> None:
        diag = {k: np.asarray(m[k]) for k in _diag_keys}
        snap = snapshot_from_states(
            step, step * config.dt, parts, grid=grid, diagnostics=diag, dx=DX
        )
        snapshots.append(str(write_snapshot(snap, snap_dir / f"snapshot_{step:05d}.h5")))

    m0 = measure(0)
    history = [m0]
    if _snapshot_due(0):
        write_snap(0, m0)
    for step in range(1, config.steps + 1):
        kdk_step(parts, acc_x, acc_y, acc_z, accel_fn, config.dt)
        want_hist = step % hist_cad == 0 or step == config.steps
        want_snap = _snapshot_due(step)
        if want_hist or want_snap:
            m = measure(step)
            if want_hist:
                history.append(m)
            if want_snap:
                write_snap(step, m)
    ti.sync()

    e0 = history[0]["energy"]
    ef = history[-1]["energy"]
    drift = abs((ef - e0) / e0) if e0 != 0.0 else float("nan")
    drift_max = max(abs((h["energy"] - e0) / e0) for h in history) if e0 != 0.0 else float("nan")
    return {
        "name": config.name,
        "backend": config.backend,
        "backend_resolved": backend["resolved"],
        "device": backend["device"],
        "preset": config.ic_preset,
        "solver": config.solver,
        "n_particles": n,
        "grid_size": gsize,
        "steps": config.steps,
        "dt": config.dt,
        "energy_initial": e0,
        "energy_final": ef,
        "energy_drift": drift,  # endpoint (first vs last sample)
        "energy_drift_max": drift_max,  # worst sample over the run
        "history": history,
        "snapshots": snapshots,
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


def run_sim_cli(argv: list[str] | None = None) -> int:
    """CLI entry point: ``galaxy-sim --config PATH [--backend NAME]``."""
    parser = argparse.ArgumentParser(
        prog="galaxy-sim",
        description="Run the Stage-3 PM N-body simulation from a YAML config.",
    )
    parser.add_argument("--config", type=Path, required=True, help="Path to a YAML run config.")
    parser.add_argument(
        "--backend", choices=tuple(_ARCH_NAMES), default=None, help="Override the config backend."
    )
    args = parser.parse_args(argv)

    config = load_config(args.config)
    if args.backend is not None:
        config.backend = args.backend
        config.validate()

    result = run_simulation(config)
    print(
        f"galaxy-sim ok: preset={result['preset']} solver={result['solver']} "
        f"n={result['n_particles']} grid={result['grid_size']}^3 steps={result['steps']} "
        f"| max energy drift={result['energy_drift_max']:.3e}, "
        f"{len(result['snapshots'])} snapshots"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(hello_sim_cli())
