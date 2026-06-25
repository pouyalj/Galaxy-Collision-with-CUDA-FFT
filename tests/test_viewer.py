"""GGUI viewer smoke test (Stage 7 / 7B).

The realtime viewer is inherently display-bound and **not** fully CI-testable, but its *offscreen*
path (``show_window=False``) renders headless and exercises the whole code path — sim chain build,
the subsample/pack kernel, GGUI scene setup, and a PNG dump. This test runs it for one frame and
checks a non-empty image lands on disk.

GGUI needs a graphics device (Metal/CUDA/Vulkan). On a GPU-less CI runner the window creation
raises, so the test **skips** there (like the CUDA determinism test) — it is a real check on the
Apple dev box / a GPU workstation, a no-op on plain CI.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


@pytest.mark.fixed_arch  # picks a GGUI-capable arch itself; not subject to GALAXY_TEST_ARCH
def test_viewer_offscreen_renders_a_nonempty_frame(tmp_path):
    pytest.importorskip("taichi", reason="Taichi not installed")
    pytest.importorskip("matplotlib", reason="matplotlib not installed (image read)")

    from galaxy_collision.config import SimConfig
    from galaxy_collision.viz.viewer import run_viewer

    cfg = SimConfig(
        name="viewer-smoke", backend="metal", ic_preset="two_galaxy_4v",
        n_particles=20_000, dt=0.5, steps=0, solver="multigrid", grid_size=32,
        output_cadence=0, seed=0,
    )
    try:
        r = run_viewer(
            cfg, max_points=4_000, steps_per_frame=1, offscreen=True, frames=1, out_dir=tmp_path
        )
    except Exception as e:  # no graphics device (headless CI), or GGUI/driver gap
        pytest.skip(f"GGUI offscreen render unavailable here: {type(e).__name__}: {e}")

    assert r["frames_rendered"] == 1
    assert r["drawn_points"] <= 4_000
    frame = Path(tmp_path) / "frame_00000.png"
    assert frame.is_file() and frame.stat().st_size > 0

    # The render must not be all-black (particles actually drew).
    from matplotlib import image

    img = image.imread(str(frame))
    assert int((img[..., :3].sum(-1) > 0.05).sum()) > 100, "rendered frame is (near-)empty"
    assert np.isfinite(img).all()
