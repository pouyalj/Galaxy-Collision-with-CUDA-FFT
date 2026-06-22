"""Static paper-reproduction figures + driver (AGENT.md §5.7, Stage 4 / 4B).

The DISLIN ``make_image`` replacement: **static** matplotlib figures of the collision —
projected-density panels over time and the Sun-like tracer's path. (The realtime Taichi GGUI
viewer and the frames→ffmpeg movie pipeline remain Stage 7; this module is only the paper's
signature still figures.)

Density is projected by 2D-histogramming the *particle* positions (mass-weighted) rather than
the coarse 1 kpc density grid, so the panels are sharper than the legacy grid projection. Uses
the headless **Agg** backend, so it renders in CI / over SSH with no display.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _plt():
    """Import pyplot on the non-interactive Agg backend (headless/CI-safe)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _hist2d(pos, mass, extent, bins, axes):
    """Mass-weighted 2D histogram of ``pos`` projected onto ``axes`` (surface density), oriented
    for ``imshow(origin="lower")``."""
    a0, a1 = axes
    h, _, _ = np.histogram2d(
        pos[:, a0],
        pos[:, a1],
        bins=bins,
        range=[[extent[0], extent[1]], [extent[2], extent[3]]],
        weights=mass,
    )
    return h.T


def projected_density(ax, pos, mass, *, extent, bins=300, axes=(0, 1), cmap="inferno", norm=None):
    """Draw a projected surface-density image of (``pos``, ``mass``) onto ``ax``. Returns the image
    (use a shared ``norm`` across panels to make them comparable). Empty cells are masked."""
    from matplotlib.colors import LogNorm

    h = _hist2d(pos, mass, extent, bins, axes)
    if norm is None:
        pos_vals = h[h > 0]
        vmin = pos_vals.min() if pos_vals.size else 1.0
        vmax = h.max() if h.max() > vmin else vmin * 10.0  # avoid a degenerate vmin==vmax LogNorm
        norm = LogNorm(vmin=vmin, vmax=vmax)
    masked = np.ma.masked_where(h <= 0, h)
    im = ax.imshow(masked, origin="lower", extent=extent, cmap=cmap, norm=norm, aspect="equal")
    return im


def density_sequence(snapshots, *, extent, bins=300, axes=(0, 1), cmap="inferno"):
    """Row of projected-density panels at increasing times — the collision sequence figure. A
    single shared log-norm (computed over all panels) makes brightness comparable across time."""
    from matplotlib.colors import LogNorm

    plt = _plt()
    hists = [_hist2d(s.pos, s.mass, extent, bins, axes) for s in snapshots]
    allpos = np.concatenate([h[h > 0] for h in hists if (h > 0).any()] or [np.array([1.0])])
    vmin = allpos.min()
    vmax = allpos.max() if allpos.max() > vmin else vmin * 10.0  # guard degenerate vmin==vmax
    norm = LogNorm(vmin=vmin, vmax=vmax)

    n = len(snapshots)
    fig, axs = plt.subplots(1, n, figsize=(3.2 * n, 3.6), squeeze=False)
    im = None
    for ax, h, snap in zip(axs[0], hists, snapshots, strict=True):
        masked = np.ma.masked_where(h <= 0, h)
        im = ax.imshow(masked, origin="lower", extent=extent, cmap=cmap, norm=norm, aspect="equal")
        ax.set_title(f"t = {snap.time:.0f} Myr")
        ax.set_xlabel("x [kpc]")
    axs[0][0].set_ylabel("y [kpc]")
    if im is not None:
        fig.colorbar(im, ax=axs[0], fraction=0.018, pad=0.01, label=r"$\Sigma$ [M$_\odot$/kpc$^2$]")
    return fig


def tracer_trajectory(
    tracer_path, *, times=None, background=None, extent=None, bins=300, axes=(0, 1), cmap="inferno"
):
    """Plot the tracer particle path(s) (``tracer_path`` shape ``(n_samples, n_tracers, 3)``),
    optionally over a faint final-state density (``background`` = a Snapshot, with ``extent``)."""
    from matplotlib.colors import LogNorm

    plt = _plt()
    a0, a1 = axes
    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    if background is not None and extent is not None:
        h = _hist2d(background.pos, background.mass, extent, bins, axes)
        masked = np.ma.masked_where(h <= 0, h)
        ax.imshow(
            masked,
            origin="lower",
            extent=extent,
            cmap=cmap,
            norm=LogNorm(),
            alpha=0.65,
            aspect="equal",
        )

    tp = np.asarray(tracer_path)
    for t in range(tp.shape[1]):
        ax.plot(tp[:, t, a0], tp[:, t, a1], "-", color="cyan", lw=1.5, alpha=0.9)
        ax.plot(
            tp[0, t, a0], tp[0, t, a1], "o", color="lime", ms=7, label="start" if t == 0 else None
        )
        ax.plot(
            tp[-1, t, a0], tp[-1, t, a1], "*", color="white", ms=14, label="end" if t == 0 else None
        )
    if extent is not None:
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    span = f" over {times[-1] - times[0]:.0f} Myr" if times is not None and len(times) > 1 else ""
    ax.set_title(f"Sun-like tracer path{span}")
    ax.legend(loc="upper right", framealpha=0.6)
    return fig


def run_paper_repro(
    config, *, tracer_radius=8.32, out_dir="figures", solver_kwargs=None, history_cadence=None
):
    """Build the two-galaxy IC (geometry from ``config.separation``/``config.impact_parameter``),
    select the Sun-like tracer, run the collision (writing snapshots), and emit the
    density-sequence + tracer-path figures. Returns ``{figures, result, tracer_index}``."""
    from galaxy_collision import ic as ic_mod
    from galaxy_collision import sim
    from galaxy_collision.io import read_snapshot

    icr = ic_mod.build_ic(
        config, separation=config.separation, impact_parameter=config.impact_parameter
    )
    tracer = ic_mod.select_tracer(icr, target_radius=tracer_radius, gid=0)
    res = sim.run_simulation(
        config,
        icr=icr,
        tracer_indices=[tracer],
        solver_kwargs=solver_kwargs,
        history_cadence=history_cadence,
    )

    plt = _plt()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    extent = (0.0, float(config.grid_size), 0.0, float(config.grid_size))
    figures: dict[str, str] = {}

    snaps = [read_snapshot(p) for p in res["snapshots"]]
    if snaps:
        fig = density_sequence(snaps, extent=extent)
        path = out / f"density_sequence_{config.ic_preset}.png"
        fig.savefig(path, dpi=130, bbox_inches="tight")
        plt.close(fig)  # release the figure (in-process batch callers would otherwise accumulate)
        figures["density_sequence"] = str(path)
    if res["tracer_path"] is not None:
        fig = tracer_trajectory(
            res["tracer_path"],
            times=res["tracer_times"],
            background=snaps[-1] if snaps else None,
            extent=extent,
        )
        path = out / f"tracer_path_{config.ic_preset}.png"
        fig.savefig(path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        figures["tracer_path"] = str(path)
    return {"figures": figures, "result": res, "tracer_index": tracer}


def paper_repro_cli(argv=None) -> int:
    """CLI: ``paper-repro --config PATH [--tracer-radius KPC] [--out DIR]``.

    Geometry (separation, impact parameter) lives in the YAML config so a run is fully
    reproducible from it; this CLI only adds the tracer radius and output location.
    """
    import argparse

    from galaxy_collision.config import load_config

    parser = argparse.ArgumentParser(
        prog="paper-repro",
        description="Run a two-galaxy collision and emit the paper-reproduction figures.",
    )
    parser.add_argument("--config", type=Path, required=True, help="Path to a YAML run config.")
    parser.add_argument("--backend", default=None, help="Override the config backend.")
    parser.add_argument(
        "--tracer-radius",
        type=float,
        default=8.32,
        help="Galactocentric radius of the Sun-like tracer (kpc).",
    )
    parser.add_argument("--out", type=Path, default=Path("figures"), help="Figure output dir.")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    if args.backend is not None:
        config.backend = args.backend
        config.validate()

    out = run_paper_repro(config, tracer_radius=args.tracer_radius, out_dir=str(args.out))
    print(
        f"paper-repro ok: preset={config.ic_preset} "
        f"tracer_idx={out['tracer_index']} figures={list(out['figures'])} -> {args.out}"
    )
    return 0


__all__ = [
    "projected_density",
    "density_sequence",
    "tracer_trajectory",
    "run_paper_repro",
    "paper_repro_cli",
]
