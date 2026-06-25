"""The device 2D density projection (Stage 7 / 7A) must match the host histogram reference.

``viz.project.project_density`` is the on-device replacement for the per-frame host histogram, so
it must reproduce ``viz.paper_repro._hist2d`` (the surface-density convention the static figures
already use) — same binning, same units, same imshow orientation.
"""

from __future__ import annotations

import numpy as np
import pytest


def test_project_density_matches_hist2d_reference():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import ParticleState
    from galaxy_collision.viz.paper_repro import _hist2d
    from galaxy_collision.viz.project import project_density

    ti.init(arch=ti.cpu)
    n, bins = 50_000, 128
    extent = (0.0, 256.0, 0.0, 256.0)
    rng = np.random.default_rng(0)
    # Particles strictly inside the extent (so right-edge binning conventions don't differ) plus a
    # few outside (must be dropped by both, as histogram2d drops out-of-range).
    pos = rng.uniform(10.0, 246.0, (n, 3)).astype(np.float32)
    pos[:200] = rng.uniform(300.0, 400.0, (200, 3))  # outside the box → dropped
    mass = rng.uniform(500.0, 1500.0, n).astype(np.float32)

    parts = ParticleState(n)
    parts.pos_x.from_numpy(pos[:, 0])
    parts.pos_y.from_numpy(pos[:, 1])
    parts.pos_z.from_numpy(pos[:, 2])
    parts.mass.from_numpy(mass)

    img = ti.field(ti.f32, shape=(bins, bins))
    for axes in [(0, 1), (0, 2), (1, 2)]:
        dev = project_density(parts, img, extent, axes=axes)
        ref = _hist2d(pos.astype(np.float64), mass.astype(np.float64), extent, bins, axes)
        assert dev.shape == ref.shape
        # f32 device accumulation vs f64 numpy — agree to f32 round-off.
        np.testing.assert_allclose(dev, ref, rtol=1e-4, atol=ref.max() * 1e-5)
