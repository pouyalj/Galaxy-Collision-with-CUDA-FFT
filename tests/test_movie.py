"""The batch→movie encoder (Stage 7 / 7A) renders frames to a playable file + a single panel.

Exercises the render + encode path on synthetic surface-density frames (no sim/GPU needed), so it
is fast and runs wherever the optional ``viz`` extra (imageio) is installed; it self-skips
otherwise. The projection that *produces* frames is covered by ``test_project.py``.
"""

from __future__ import annotations

import numpy as np
import pytest


def test_write_movie_and_panel(tmp_path):
    pytest.importorskip("imageio", reason="viz extra (imageio) not installed")
    pytest.importorskip("matplotlib")
    from galaxy_collision.viz.movie import save_panel, write_movie

    rng = np.random.default_rng(0)
    frames = np.abs(rng.normal(size=(5, 48, 48))) * 1e7  # surface densities (M_sun/kpc^2)
    times = np.array([0.0, 10.0, 20.0, 30.0, 40.0])
    extent = (0.0, 256.0, 0.0, 256.0)

    mp4 = write_movie(frames, times, tmp_path / "m.mp4", extent=extent, fps=8)
    assert mp4.exists() and mp4.stat().st_size > 0

    gif = write_movie(frames, times, tmp_path / "m.gif", extent=extent, fps=8)
    assert gif.exists() and gif.stat().st_size > 0

    panel = save_panel(frames[-1], tmp_path / "p.png", extent=extent, time=40.0)
    assert panel.exists() and panel.stat().st_size > 0
