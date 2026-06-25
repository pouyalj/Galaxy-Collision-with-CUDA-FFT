"""Batch → movie: device-projected density frames → MP4/GIF (Stage 7 / 7A).

The frames come from ``run_simulation(frame_cadence=…)`` (device-projected surface-density images,
see :mod:`viz.project`), so this module never touches the multi-GB position snapshots. It renders
each frame with the **same** LogNorm + colormap + surface-density convention as the static paper
figures (:mod:`viz.paper_repro`), with one global color scale across the sequence so brightness is
comparable over time, then encodes via ``imageio`` (``imageio-ffmpeg`` bundles ffmpeg — no system
install needed). If the encoder is unavailable it degrades to writing PNG frames + the ffmpeg
command. This is the DISLIN replacement of the original 2020 code.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def _global_lognorm(frames: np.ndarray):
    """A single LogNorm spanning the positive surface-density values across all frames."""
    from matplotlib.colors import LogNorm

    pos = frames[frames > 0]
    vmin = float(pos.min()) if pos.size else 1.0
    vmax = float(frames.max()) if frames.max() > vmin else vmin * 10.0  # guard degenerate
    return LogNorm(vmin=vmin, vmax=vmax)


def _render_frame_rgb(img, *, extent, norm, cmap, title, figsize=(6.0, 5.4), dpi=100) -> np.ndarray:
    """Render one surface-density image to an (H, W, 3) uint8 array (matplotlib Agg, off-screen).

    Fixed figsize×dpi (no ``bbox_inches="tight"``) so every frame has identical pixel dimensions —
    a hard requirement for video encoders.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=figsize, dpi=dpi)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes((0.12, 0.12, 0.74, 0.80))
    masked = np.ma.masked_where(img <= 0, img)
    im = ax.imshow(masked, origin="lower", extent=extent, cmap=cmap, norm=norm, aspect="equal")
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_title(title)
    cax = fig.add_axes((0.88, 0.12, 0.03, 0.80))
    fig.colorbar(im, cax=cax, label=r"$\Sigma$ [M$_\odot$/kpc$^2$]")
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[:, :, :3].copy()


def save_panel(img, path, *, extent, time=None, cmap="inferno") -> Path:
    """Save a single surface-density panel as a PNG (e.g. the RV19 100M headline figure)."""
    import imageio.v2 as imageio

    title = "surface density" if time is None else f"t = {time:.0f} Myr"
    rgb = _render_frame_rgb(img, extent=extent, norm=_global_lognorm(np.asarray(img)[None]),
                            cmap=cmap, title=title)
    path = Path(path)
    imageio.imwrite(path, rgb)
    return path


def write_movie(frames, times, path, *, extent, fps=12, cmap="inferno") -> Path:
    """Encode ``frames`` (``(n, bins, bins)`` surface-density) to an MP4/GIF at ``path``.

    Falls back to a directory of PNG frames + the ffmpeg command if no encoder is installed.
    """
    frames = np.asarray(frames)
    times = np.asarray(times)
    norm = _global_lognorm(frames)
    title_at = (lambda k: f"t = {times[k]:.0f} Myr") if times.size else (lambda k: "")
    rgb = [_render_frame_rgb(frames[k], extent=extent, norm=norm, cmap=cmap, title=title_at(k))
           for k in range(len(frames))]

    path = Path(path)
    try:
        import imageio.v2 as imageio

        with imageio.get_writer(path, **_writer_kwargs(path, fps)) as w:
            for frame in rgb:
                w.append_data(frame)
        return path
    except Exception as exc:  # graceful degrade: dump frames + the ffmpeg command
        frame_dir = path.with_suffix("")
        frame_dir.mkdir(parents=True, exist_ok=True)
        import imageio.v2 as imageio

        for k, frame in enumerate(rgb):
            imageio.imwrite(frame_dir / f"frame_{k:04d}.png", frame)
        print(f"movie encode failed ({exc}); wrote PNG frames to {frame_dir}. Stitch with:\n"
              f"  ffmpeg -framerate {fps} -i {frame_dir}/frame_%04d.png "
              f"-pix_fmt yuv420p {path}")
        return frame_dir


def _writer_kwargs(path: Path, fps: int) -> dict[str, Any]:
    """imageio writer kwargs by container: mp4 takes fps (+ macro_block_size pads to even dims for
    libx264); the GIF (Pillow) writer wants per-frame duration in ms, not fps."""
    if path.suffix.lower() == ".gif":
        return {"duration": 1000.0 / fps}
    if path.suffix.lower() == ".mp4":
        return {"fps": fps, "macro_block_size": 16}
    return {"fps": fps}


def encode_image_files(image_paths, path, *, fps=12) -> Path:
    """Encode existing image files (e.g. the viewer's GGUI offscreen 3D PNGs) into an MP4/GIF —
    the encode half of the 3D (``--view particles``) movie path, where the frames are already
    rendered RGB images rather than density arrays."""
    import imageio.v2 as imageio

    path = Path(path)
    with imageio.get_writer(path, **_writer_kwargs(path, fps)) as w:
        for p in image_paths:
            w.append_data(imageio.imread(p))
    return path


_AXES = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}


def _movie_particles_3d(config, args) -> int:
    """3D movie: render the collision with the viewer's GGUI offscreen particle renderer (same
    look as ``galaxy-view``'s 3D mode) and encode the frames. Reuses the proven offscreen path —
    it renders correctly-oriented PNGs to a temp dir, which we then stitch into the movie."""
    import shutil
    import tempfile

    from galaxy_collision.viz.viewer import run_viewer

    cadence = args.frame_cadence or max(1, config.steps // 60)
    n_frames = max(1, config.steps // cadence)
    with tempfile.TemporaryDirectory() as td:
        run_viewer(
            config, offscreen=True, frames=n_frames, steps_per_frame=cadence,
            out_dir=td, max_points=args.max_points, point_radius=args.radius,
            color_by=args.color,
        )
        pngs = sorted(Path(td).glob("frame_*.png"))
        out = encode_image_files(pngs, args.out, fps=args.fps)
        if args.panel is not None and pngs:
            shutil.copyfile(pngs[-1], args.panel)
    print(f"galaxy-movie ok: {len(pngs)} 3D frames (≤{args.max_points:,} pts) -> {out}")
    if args.panel is not None:
        print(f"  panel -> {args.panel}")
    return 0


def galaxy_movie_cli(argv: list[str] | None = None) -> int:
    """CLI: ``galaxy-movie --config C [--view density|particles] [--backend B] [--n N] …``."""
    import argparse

    from galaxy_collision.config import load_config
    from galaxy_collision.sim import run_simulation

    p = argparse.ArgumentParser(
        prog="galaxy-movie",
        description="Run a sim and render a collision movie (2D density or 3D particles).")
    p.add_argument("--config", type=Path, required=True, help="YAML run config.")
    p.add_argument("--view", choices=("density", "particles"), default="density",
                   help="density = 2D surface-density projection (default); "
                        "particles = 3D point cloud, like galaxy-view's 3D mode.")
    p.add_argument("--backend", choices=("cpu", "cuda", "metal"), default=None)
    p.add_argument("--n", type=int, default=None,
                   help="particle count, overriding the config's n_particles (e.g. 100000000).")
    p.add_argument("--frame-cadence", type=int, default=None,
                   help="steps between frames (default: ~60 frames across the run).")
    p.add_argument("--bins", type=int, default=512, help="2D image resolution per side (density).")
    p.add_argument("--axes", choices=tuple(_AXES), default="xy", help="projection plane (density).")
    p.add_argument("--max-points", type=int, default=300_000,
                   help="3D: cap on drawn (subsampled) points (particles view).")
    p.add_argument("--radius", type=float, default=0.0,
                   help="3D: point radius (0 = auto; particles view).")
    p.add_argument("--color", choices=("galaxy", "density"), default="galaxy",
                   help="3D point color: galaxy id (default) or local density (3D density cloud).")
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--out", type=Path, default=Path("collision.mp4"))
    p.add_argument("--panel", type=Path, default=None, help="also save the final frame as a PNG.")
    args = p.parse_args(argv)

    config = load_config(args.config)
    if args.backend is not None:
        config.backend = args.backend
    if args.n is not None:
        # Override the config's particle count; clear particle_mass so the two don't both
        # resolve a count (SimConfig forbids setting both).
        config.n_particles = args.n
        config.particle_mass = None
    config.validate()

    if args.view == "particles":
        return _movie_particles_3d(config, args)

    cadence = args.frame_cadence or max(1, config.steps // 60)
    result = run_simulation(
        config, write_snapshots=False,
        frame_cadence=cadence, frame_bins=args.bins, frame_axes=_AXES[args.axes],
    )
    frames, times, extent = result["frames"], result["frame_times"], result["frame_extent"]
    out = write_movie(frames, times, args.out, extent=extent, fps=args.fps)
    print(f"galaxy-movie ok: {len(frames)} frames @ {args.bins}^2 ({args.axes}) -> {out}")
    if args.panel is not None:
        save_panel(frames[-1], args.panel, extent=extent, time=float(times[-1]))
        print(f"  panel -> {args.panel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(galaxy_movie_cli())
