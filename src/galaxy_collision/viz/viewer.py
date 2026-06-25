"""Realtime 3D/2D particle viewer (Taichi GGUI) — Stage 7 / 7B, AGENT.md §5.7 (D25–D27).

Watch the collision evolve live on the GPU: the same device-resident force chain
``run_simulation`` runs (deposit → solve → grad → gather → KDK), rendered each frame by Taichi
**GGUI**. The renderer runs in-process on the *same* arch as the sim — confirmed on Apple
**Metal** (and CUDA / Vulkan), so on a Mac the whole thing (physics + render) stays on the GPU.

Two render modes, toggled live with **M**:

* **3D** — ``scene.particles`` over a fixed **subsample** (drawing 10–100M points is neither
  feasible nor legible), colored by galaxy id; orbit the camera with RMB-drag + WASD.
* **2D** — the full-N device density projection (``viz.project``, the same kernel the movie uses) →
  ``canvas.set_image``; cycle the projection plane (xy/xz/yz) with **P**.

Interactive controls: **SPACE** pause · **N** single-step (while paused) · **R** restart to t=0 ·
**[ / ]** slower/faster (steps per frame) · **- / =** smaller/larger points (3D) · **ESC/Q** quit.

One renderer, two run modes: **interactive** (``offscreen=False``, needs a display — *not*
CI-tested) and **offscreen** (``offscreen=True``, ``show_window=False`` — renders ``frames`` PNGs
headless; the CI smoke path). The physics always runs at the config's full N; only the *drawn* 3D
set is thinned.

This module deliberately omits ``from __future__ import annotations`` — PEP 563 would stringize
the ``ti.template()`` kernel annotations (see ``deposit.py``).
"""

import math
from pathlib import Path

import numpy as np
import taichi as ti

# 2D-mode display: fixed log dynamic range (decades below the per-frame peak) for the density
# colormap. A fixed window is stable frame-to-frame (a per-frame min would flicker) and standard
# for density displays; the movie/paper figures use their own LogNorm and are unaffected.
_LOG_DECADES = 4.0
_LUT_N = 256  # colormap lookup-table resolution

# Per-galaxy point colors (gid 0 = Milky Way, 1 = Andromeda); anything else → grey.
_COLOR_MW = (0.55, 0.72, 1.0)   # cool white-blue
_COLOR_AND = (1.0, 0.66, 0.28)  # warm orange
_COLOR_OTHER = (0.7, 0.7, 0.7)

# 2D-mode projection planes cycled by 'P': (in-plane axis a, b), projecting out the third.
_PROJ_PLANES = ((0, 1), (0, 2), (1, 2))  # xy, xz, yz
_PLANE_NAMES = {(0, 1): "xy", (0, 2): "xz", (1, 2): "yz"}


@ti.kernel
def _pack_positions(
    pos_x: ti.template(), pos_y: ti.template(), pos_z: ti.template(),
    stride: ti.i32, inv_box: ti.f32, m: ti.i32, out_pos: ti.template(),
):
    """Subsample (every ``stride``-th particle) and normalize into a centered unit cube.

    The box spans ``[0, grid_size]`` kpc per axis; ``inv_box = 1/grid_size`` maps it to ``[0,1]``,
    then we shift by −0.5 so the box is centered on the origin for a stable camera."""
    for i in range(m):
        p = i * stride
        out_pos[i] = ti.Vector([
            pos_x[p] * inv_box - 0.5,
            pos_y[p] * inv_box - 0.5,
            pos_z[p] * inv_box - 0.5,
        ])


@ti.kernel
def _max_field(img: ti.template(), out: ti.template()):
    """Device max-reduction of ``img`` into the 0-D field ``out`` (f32, Metal-safe; 0-D → TLS)."""
    out[None] = 0.0
    for i, j in img:
        ti.atomic_max(out[None], img[i, j])


@ti.kernel
def _colormap_log(
    img: ti.template(), lut: ti.template(), log_vmin: ti.f32, inv_logspan: ti.f32,
    out: ti.template(),
):
    """LogNorm a density image through a 256-entry colormap LUT — **fully on device**.

    ``out[i,j] = lut[ round(clamp01((ln img − log_vmin)·inv_logspan) · 255) ]``. Indexed ``[i,j]`` =
    ``[x,y]`` to match both ``project._project`` and GGUI ``set_image`` (no transpose). Only the
    scalar peak crosses to the host (to derive ``log_vmin``); the image never does."""
    for i, j in img:
        v = img[i, j]
        t = 0.0
        if v > 0.0:
            t = (ti.log(v) - log_vmin) * inv_logspan
        t = ti.min(1.0, ti.max(0.0, t))
        idx = ti.cast(t * (_LUT_N - 1) + 0.5, ti.i32)
        idx = ti.min(_LUT_N - 1, ti.max(0, idx))
        out[i, j] = ti.Vector([lut[idx, 0], lut[idx, 1], lut[idx, 2]])


@ti.kernel
def _max_grid3(field: ti.template(), out: ti.template()):
    """Device max-reduction of a 3D field into 0-D ``out`` (Metal-safe; 0-D → TLS)."""
    out[None] = 0.0
    for i, j, k in field:
        ti.atomic_max(out[None], field[i, j, k])


@ti.kernel
def _density_color(
    pos_x: ti.template(), pos_y: ti.template(), pos_z: ti.template(),
    rho: ti.template(), stride: ti.i32, m: ti.i32, inv_dx: ti.f32, gsize: ti.i32,
    lut: ti.template(), log_vmin: ti.f32, inv_logspan: ti.f32, out_col: ti.template(),
):
    """Color each drawn (subsampled) point by the **local 3D density** ``rho`` at its position,
    through the same inferno LogNorm LUT as the 2D view — turning the 3D point cloud into a
    *density cloud*. Nearest-grid-node gather (cheap; the grid is 1 kpc, fine for coloring)."""
    for i in range(m):
        p = i * stride
        ix = ti.min(gsize - 1, ti.max(0, ti.cast(pos_x[p] * inv_dx + 0.5, ti.i32)))
        iy = ti.min(gsize - 1, ti.max(0, ti.cast(pos_y[p] * inv_dx + 0.5, ti.i32)))
        iz = ti.min(gsize - 1, ti.max(0, ti.cast(pos_z[p] * inv_dx + 0.5, ti.i32)))
        v = rho[ix, iy, iz]
        t = 0.0
        if v > 0.0:
            t = (ti.log(v) - log_vmin) * inv_logspan
        t = ti.min(1.0, ti.max(0.0, t))
        idx = ti.min(_LUT_N - 1, ti.max(0, ti.cast(t * (_LUT_N - 1) + 0.5, ti.i32)))
        out_col[i] = ti.Vector([lut[idx, 0], lut[idx, 1], lut[idx, 2]])


def _inferno_lut() -> np.ndarray:
    """The 256×3 'inferno' colormap as an f32 LUT (one-time host build, uploaded to the device)."""
    import matplotlib

    return matplotlib.colormaps["inferno"](np.linspace(0.0, 1.0, _LUT_N))[:, :3].astype(np.float32)


def _setup_chain(config):
    """Build the device-resident force chain for ``config`` (mirrors ``bench.run_benchmark``)."""
    from galaxy_collision import ic as ic_mod
    from galaxy_collision import units
    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.sim import DX, _build_ic
    from galaxy_collision.solver import make_solver

    icr = _build_ic(config)
    n, gs = icr.n, config.grid_size
    parts = ParticleState(n)
    ic_mod.load_into_particle_state(icr, parts)
    grid = GridState(gs)
    acc = tuple(ti.field(ti.f32, shape=n) for _ in range(3))
    solver = make_solver(config.solver, gs, dx=DX, grav_constant=units.G)
    # Launch the disks in equilibrium so the opening frames aren't a cold transient (as bench does).
    if config.ic_preset in ("two_galaxy_4v", "two_galaxy_2v"):
        ic_mod.equilibrate_disk_velocities(
            icr, dx=DX, solver=solver, rho=grid.rho, phi=grid.phi,
            ax=grid.ax, ay=grid.ay, az=grid.az, parts=parts,
        )
    return icr, parts, grid, acc, solver, n, gs, DX


def run_viewer(
    config,
    *,
    max_points: int = 300_000,
    steps_per_frame: int = 1,
    point_radius: float = 0.0,
    bins: int = 512,
    offscreen: bool = False,
    frames: int = 0,
    out_dir: str | Path = "frames",
    resolution=(1024, 768),
    color_by: str = "galaxy",
):
    """Run the PM sim and render its particles with GGUI (3D points + 2D density projection).

    ``max_points`` caps the *drawn* 3D particle count (subsampled); the physics runs at full N.
    ``steps_per_frame`` advances the sim that many KDK steps between draws. ``bins`` is the 2D
    projection resolution. Interactive by default; ``offscreen=True`` renders ``frames`` PNGs into
    ``out_dir`` with no window (headless, 3D mode). Returns a small summary dict.
    """
    from galaxy_collision.deposit import deposit_density, gather_acceleration, potential_to_accel
    from galaxy_collision.integrator import kdk_step
    from galaxy_collision.sim import init_backend
    from galaxy_collision.viz.project import scatter_density

    backend = init_backend(config)
    icr, parts, grid, acc, solver, n, gs, dx = _setup_chain(config)
    acc_x, acc_y, acc_z = acc
    # Stash the t=0 state (host arrays) so 'R' can restart without reallocating device fields.
    pos0, vel0 = icr.pos.copy(), icr.vel.copy()

    # 3D render fields: ceil-division stride so `max_points` is a true cap (m ≤ max_points).
    stride = max(1, -(-n // max_points))
    m = len(range(0, n, stride))
    rpos = ti.Vector.field(3, ti.f32, shape=m)
    rcol = ti.Vector.field(3, ti.f32, shape=m)
    gid_sub = icr.gid[::stride][:m].astype(np.int32)
    col = np.empty((m, 3), dtype=np.float32)
    col[:] = _COLOR_OTHER
    col[gid_sub == 0] = _COLOR_MW
    col[gid_sub == 1] = _COLOR_AND
    rcol.from_numpy(col)
    inv_box = 1.0 / float(gs)
    radius = point_radius if point_radius > 0.0 else max(0.0012, 0.06 / m**0.5)

    # 2D mode fields (all device-resident): the projection target, a colormap LUT, a 0-D peak
    # buffer, and the RGB image fed to set_image. Only the scalar peak crosses to the host.
    proj_img = ti.field(ti.f32, shape=(bins, bins))
    rgb_img = ti.Vector.field(3, ti.f32, shape=(bins, bins))
    lut = ti.field(ti.f32, shape=(_LUT_N, 3))
    lut.from_numpy(_inferno_lut())
    maxbuf = ti.field(ti.f32, shape=())
    inv_logspan = 1.0 / (_LOG_DECADES * math.log(10.0))
    extent_full = (0.0, float(gs), 0.0, float(gs))

    def reload_ic():
        from galaxy_collision import ic as ic_mod
        icr.pos[:], icr.vel[:] = pos0, vel0
        ic_mod.load_into_particle_state(icr, parts)

    warm = {"on": False}

    def accel_fn():
        """Recompute a(x) into ``acc_*`` from current positions — the one force eval ``kdk_step``
        calls each step (and ``prime`` calls once up front). Warm-starts the solve after step 1."""
        deposit_density(parts, grid.rho, dx)
        solver.solve(grid.rho, grid.phi, warm_start=warm["on"])
        warm["on"] = True
        potential_to_accel(grid.phi, grid.ax, grid.ay, grid.az, dx)
        gather_acceleration(parts, grid.ax, grid.ay, grid.az, acc_x, acc_y, acc_z, dx)

    def prime():
        warm["on"] = False  # first solve is cold; subsequent ones warm-start
        accel_fn()

    def step_once():
        # Reuse the canonical KDK ordering (integrator.kdk_step) — no inline re-implementation.
        kdk_step(parts, acc_x, acc_y, acc_z, accel_fn, config.dt)

    prime()

    window = ti.ui.Window("Galaxy Collision", resolution, show_window=not offscreen, vsync=True)
    canvas = window.get_canvas()
    canvas.set_background_color((0.0, 0.0, 0.02))
    scene = window.get_scene()
    camera = ti.ui.Camera()
    camera.position(0.0, 0.0, 1.6)
    camera.lookat(0.0, 0.0, 0.0)
    camera.up(0.0, 1.0, 0.0)

    # Mutable view state driven by the keyboard.
    st = {"mode": "3d", "paused": False, "spf": max(1, steps_per_frame), "plane": 0,
          "radius": radius, "step": 0, "single": False,
          "color": color_by if color_by in ("galaxy", "density") else "galaxy"}

    def restart():
        reload_ic()
        prime()  # re-primes acc_* (and resets the warm-start) for the reloaded positions
        st["step"] = 0

    def draw_3d():
        _pack_positions(parts.pos_x, parts.pos_y, parts.pos_z, stride, inv_box, m, rpos)
        if st["color"] == "density":
            # Recolor the cloud by local density each frame (a 3D "density cloud"): grid max →
            # LogNorm window → per-point nearest-node gather through the inferno LUT.
            _max_grid3(grid.rho, maxbuf)
            vmax = float(maxbuf[None]) or 1.0
            log_vmin = math.log(vmax) - _LOG_DECADES * math.log(10.0)
            _density_color(parts.pos_x, parts.pos_y, parts.pos_z, grid.rho, stride, m,
                           1.0 / dx, gs, lut, log_vmin, inv_logspan, rcol)
        scene.set_camera(camera)
        scene.ambient_light((0.45, 0.45, 0.5))
        scene.point_light(pos=(2.0, 2.0, 2.0), color=(1.0, 1.0, 1.0))
        scene.particles(rpos, radius=st["radius"], per_vertex_color=rcol)
        canvas.scene(scene)

    def draw_2d():
        # Fully device-resident: scatter → device max → device LogNorm+colormap → set_image.
        # Only the scalar peak crosses to the host (to set the LogNorm window); the image does not.
        scatter_density(parts, proj_img, extent_full, axes=_PROJ_PLANES[st["plane"]])
        _max_field(proj_img, maxbuf)
        vmax = float(maxbuf[None])
        if vmax <= 0.0:
            vmax = 1.0
        log_vmin = math.log(vmax) - _LOG_DECADES * math.log(10.0)
        _colormap_log(proj_img, lut, log_vmin, inv_logspan, rgb_img)
        canvas.set_image(rgb_img)

    def advance():
        if st["paused"] and not st["single"]:
            return
        nsteps = 1 if st["single"] else st["spf"]
        for _ in range(nsteps):
            step_once()
            st["step"] += 1
        st["single"] = False

    drawn = 0
    if offscreen:
        from galaxy_collision.progress import ProgressBar

        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        nframes = frames if frames > 0 else 1
        bar = ProgressBar(nframes, label="rendering")
        for f in range(nframes):
            for _ in range(st["spf"]):
                step_once()
                st["step"] += 1
            draw_3d()
            window.save_image(str(out / f"frame_{f:05d}.png"))
            drawn += 1
            bar.update(f + 1)
        bar.finish()
        return _summary(backend, n, m, stride, drawn, offscreen)

    while window.running:
        for e in window.get_events(ti.ui.PRESS):
            k = e.key
            if k in (ti.ui.ESCAPE, "q"):
                window.running = False
            elif k == ti.ui.SPACE:
                st["paused"] = not st["paused"]
            elif k == "n":
                st["single"] = True
            elif k == "r":
                restart()
            elif k == "m":
                st["mode"] = "2d" if st["mode"] == "3d" else "3d"
            elif k == "p":
                st["plane"] = (st["plane"] + 1) % len(_PROJ_PLANES)
            elif k == "c":
                st["color"] = "density" if st["color"] == "galaxy" else "galaxy"
                if st["color"] == "galaxy":
                    rcol.from_numpy(col)  # restore the static galaxy-id colors
            elif k == "[":
                st["spf"] = max(1, st["spf"] - 1)
            elif k == "]":
                st["spf"] += 1
            elif k == "-":
                st["radius"] = max(0.0005, st["radius"] * 0.8)
            elif k == "=":
                st["radius"] *= 1.25
        if st["mode"] == "3d":
            camera.track_user_inputs(window, movement_speed=0.02, hold_key=ti.ui.RMB)

        advance()
        draw_3d() if st["mode"] == "3d" else draw_2d()

        gui = window.get_gui()
        with gui.sub_window("Galaxy Collision", 0.02, 0.02, 0.34, 0.22):
            gui.text(f"step {st['step']}   t = {st['step'] * config.dt:.1f} Myr")
            gui.text(f"mode {st['mode'].upper()}"
                     + (f"  plane {_PLANE_NAMES[_PROJ_PLANES[st['plane']]]}"
                        if st["mode"] == "2d" else f"  color {st['color']}  drawn {m:,}/{n:,}"))
            gui.text(f"speed {st['spf']} steps/frame" + ("   [PAUSED]" if st["paused"] else ""))
            gui.text("SPACE pause  N step  R restart  M 2D/3D  C color  P plane  [ ] speed  ESC")
        window.show()
        drawn += 1

    return _summary(backend, n, m, stride, drawn, offscreen)


def _summary(backend, n, m, stride, drawn, offscreen):
    return {
        "device": backend["device"],
        "backend_resolved": backend["resolved"],
        "n_particles": n,
        "drawn_points": m,
        "stride": stride,
        "frames_rendered": drawn,
        "offscreen": offscreen,
    }


def view_cli(argv: list[str] | None = None) -> int:
    """CLI: ``galaxy-view --config C [--backend B] [--max-points N] [--steps-per-frame K]
    [--bins B] [--offscreen --frames N --out DIR]``."""
    import argparse

    from galaxy_collision.config import load_config

    p = argparse.ArgumentParser(
        prog="galaxy-view",
        description="Realtime GGUI 3D/2D viewer for a PM N-body run (RMB-drag orbit; M = 2D).",
    )
    p.add_argument("--config", type=Path, required=True, help="Path to a YAML run config.")
    p.add_argument("--backend", default=None, help="Override the config backend (cpu|cuda|metal).")
    p.add_argument("--max-points", type=int, default=300_000,
                   help="Cap on drawn (subsampled) 3D particles; physics still runs at full N.")
    p.add_argument("--steps-per-frame", type=int, default=1,
                   help="KDK steps advanced between rendered frames.")
    p.add_argument("--bins", type=int, default=512, help="2D density-projection resolution.")
    p.add_argument("--radius", type=float, default=0.0, help="3D point radius (0 = auto).")
    p.add_argument("--color", choices=("galaxy", "density"), default="galaxy",
                   help="3D point color: galaxy id (default) or local density (density cloud; "
                        "toggle live with C).")
    p.add_argument("--offscreen", action="store_true",
                   help="Headless: render PNG frames with no window (CI / quick check).")
    p.add_argument("--frames", type=int, default=0, help="Frames to render in --offscreen mode.")
    p.add_argument("--out", type=Path, default=Path("frames"), help="Output dir for --offscreen.")
    args = p.parse_args(argv)

    config = load_config(args.config)
    if args.backend is not None:
        config.backend = args.backend
        config.validate()

    r = run_viewer(
        config, max_points=args.max_points, steps_per_frame=args.steps_per_frame,
        point_radius=args.radius, bins=args.bins, offscreen=args.offscreen,
        frames=args.frames, out_dir=args.out, color_by=args.color,
    )
    print(f"galaxy-view: {r['device']} | N={r['n_particles']:,} drawn={r['drawn_points']:,} "
          f"(stride {r['stride']}) | frames={r['frames_rendered']} offscreen={r['offscreen']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(view_cli())
