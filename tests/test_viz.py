"""Paper-reproduction figures + driver (Stage 4 / 4B).

The figure functions are Taichi-free — they take plain arrays / :class:`Snapshot` objects — so
they test without a GPU runtime; only the end-to-end driver smoke needs Taichi. All rendering
uses the headless Agg backend (see paper_repro._plt), so these run in CI with no display.
"""

from __future__ import annotations

import numpy as np
import pytest

mpl = pytest.importorskip("matplotlib", reason="matplotlib not installed")


def _two_blobs(n=4000, box=160.0, seed=0):
    """Two Gaussian blobs (fake galaxies) as (pos, mass, gid)."""
    rng = np.random.default_rng(seed)
    c0 = np.array([box * 0.35, box * 0.5, box * 0.5])
    c1 = np.array([box * 0.65, box * 0.5, box * 0.5])
    p0 = rng.normal(c0, 8.0, (n // 2, 3))
    p1 = rng.normal(c1, 8.0, (n // 2, 3))
    pos = np.concatenate([p0, p1]).clip(0, box)
    mass = np.full(n, 1.0e6)
    gid = np.concatenate([np.zeros(n // 2, np.int32), np.ones(n // 2, np.int32)])
    return pos, mass, gid


def _snapshot(step, time, pos, mass, gid):
    from galaxy_collision.io import Snapshot

    return Snapshot(step=step, time=time, pos=pos, vel=np.zeros_like(pos), mass=mass, gid=gid)


def test_projected_density_returns_masked_image():
    from galaxy_collision.viz.paper_repro import _plt, projected_density

    plt = _plt()
    pos, mass, _ = _two_blobs()
    fig, ax = plt.subplots()
    im = projected_density(ax, pos, mass, extent=(0, 160, 0, 160), bins=64)
    arr = im.get_array()
    assert arr.shape == (64, 64)
    assert arr.count() > 0  # some unmasked (non-empty) cells
    plt.close(fig)


def test_density_sequence_saves_figure(tmp_path):
    from galaxy_collision.viz.paper_repro import density_sequence

    snaps = [_snapshot(s, float(s), *_two_blobs(seed=s)) for s in (0, 5, 10)]
    fig = density_sequence(snaps, extent=(0, 160, 0, 160), bins=64)
    assert len(fig.axes) >= 3  # 3 panels (+ a colorbar axis)
    out = tmp_path / "seq.png"
    fig.savefig(out)
    assert out.is_file() and out.stat().st_size > 2000
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_tracer_trajectory_saves_figure(tmp_path):
    from galaxy_collision.viz.paper_repro import tracer_trajectory

    # A short curved path for one tracer: (n_samples, 1, 3).
    t = np.linspace(0, 1, 20)
    path = np.stack([56 + 4 * np.cos(6 * t), 80 + 4 * np.sin(6 * t), np.full_like(t, 80.0)], axis=1)
    path = path[:, None, :]
    pos, mass, gid = _two_blobs()
    bg = _snapshot(10, 50.0, pos, mass, gid)
    fig = tracer_trajectory(
        path, times=np.linspace(0, 50, 20), background=bg, extent=(0, 160, 0, 160), bins=64
    )
    out = tmp_path / "tracer.png"
    fig.savefig(out)
    assert out.is_file() and out.stat().st_size > 2000
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_run_paper_repro_smoke(tmp_path):
    """End-to-end driver: tuned IC -> tracer select -> run with snapshots -> both figures."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision.config import SimConfig
    from galaxy_collision.viz.paper_repro import run_paper_repro

    cfg = SimConfig(
        name="repro-smoke",
        backend="cpu",
        ic_preset="two_galaxy_4v",
        solver="multigrid",
        grid_size=160,
        n_particles=4000,
        dt=1.0,
        steps=12,
        output_cadence=4,
        output_dir=str(tmp_path / "out"),
        seed=0,
    )
    out = run_paper_repro(
        cfg,
        out_dir=str(tmp_path / "figs"),
        solver_kwargs={"n_cycles": 8},
        history_cadence=2,
    )
    assert set(out["figures"]) == {"density_sequence", "tracer_path"}
    for path in out["figures"].values():
        from pathlib import Path

        assert Path(path).is_file()
    res = out["result"]
    # One tracer, recorded at every history sample; it is a galaxy-0 disk particle.
    assert res["tracer_path"].shape[1:] == (1, 3)
    assert res["tracer_path"].shape[0] == len(res["history"])
    assert res["tracer_indices"] == [out["tracer_index"]]
