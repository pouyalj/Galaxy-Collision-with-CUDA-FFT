"""Device-side 2D density projection (Stage 7 / 7A).

Rendering a density frame at 100M by pulling all positions to the host is ~2.4 GB/frame — a
non-starter for a movie or a live 2D view. This module bins particle mass into a small ``(bins,
bins)`` image **on the device** (a deposit-style atomic scatter, one axis projected out); only
that image (~1 MB at 512²) crosses to the host. The same kernel feeds both the batch→movie path
and the realtime viewer's 2D mode (Stage 7B).

The output matches ``viz.paper_repro._hist2d`` exactly — a mass-weighted 2D histogram divided by
the bin area (surface density, M_sun/kpc²), transposed for ``imshow(origin="lower")`` — so movie
frames and the static paper figures share one convention. Binning is **nearest-bin (NGP), not
CIC** — it must match ``np.histogram2d``, and a LogNorm image doesn't benefit from CIC smoothing.
Mass outside the extent is dropped (open boundaries, like the deposit); a particle exactly on the
*upper* edge lands in ``bin == bins`` and is likewise dropped, where ``np.histogram2d`` keeps it in
the last bin — a measure-zero difference for real particle distributions.

f32 accumulation keeps it Metal-safe (no fp64) and is plenty for a LogNorm image. This module
deliberately omits ``from __future__ import annotations`` — PEP 563 would stringize the
``ti.template()`` kernel annotations (see ``deposit.py``).
"""

import numpy as np
import taichi as ti


@ti.kernel
def _project(
    pos_a: ti.template(),
    pos_b: ti.template(),
    mass: ti.template(),
    n: ti.i32,
    lo_a: ti.f32,
    lo_b: ti.f32,
    binw_a: ti.f32,
    binw_b: ti.f32,
    bins: ti.i32,
    img: ti.template(),
):
    """Scatter ``mass`` into ``img[ia, ib]`` (mass per bin) over the in-plane axes ``a``, ``b``."""
    for idx in ti.grouped(img):
        img[idx] = 0.0
    for p in range(n):
        ia = ti.cast(ti.floor((pos_a[p] - lo_a) / binw_a), ti.i32)
        ib = ti.cast(ti.floor((pos_b[p] - lo_b) / binw_b), ti.i32)
        if 0 <= ia < bins and 0 <= ib < bins:
            ti.atomic_add(img[ia, ib], mass[p])


# Which ParticleState position field backs each axis index.
_AXIS = {0: "pos_x", 1: "pos_y", 2: "pos_z"}


def scatter_density(parts, img, extent, axes=(0, 1)) -> tuple[float, float]:
    """Device-only: scatter ``parts``' mass into ``img`` (``(bins, bins)`` mass-per-bin) and
    return the in-plane bin widths ``(binw_a, binw_b)``. **No host copy** — ``img`` stays on the
    device for a GPU consumer (the viewer's 2D colormap, Stage 7B). ``project_density`` wraps this
    and adds the host transfer + surface-density conversion.

    Note the field is indexed ``img[ia, ib]`` = ``[a, b]`` natively (no transpose), which is exactly
    the ``[x, y]`` layout GGUI ``canvas.set_image`` expects — so a device consumer needs no flip."""
    a, b = axes
    bins = img.shape[0]
    lo_a, hi_a, lo_b, hi_b = (float(v) for v in extent)
    binw_a = (hi_a - lo_a) / bins
    binw_b = (hi_b - lo_b) / bins
    _project(
        getattr(parts, _AXIS[a]), getattr(parts, _AXIS[b]), parts.mass, parts.n,
        lo_a, lo_b, binw_a, binw_b, bins, img,
    )
    return binw_a, binw_b


def project_density(parts, img, extent, axes=(0, 1)) -> np.ndarray:
    """Fill ``img`` (a preallocated square ``(bins, bins)`` Taichi f32 field) with the device
    projection of ``parts`` and return the **surface density** as a NumPy array.

    ``extent = (a_lo, a_hi, b_lo, b_hi)`` are the in-plane ranges (kpc) for axes ``a, b`` = ``axes``
    (e.g. ``(0, 1)`` = the x–y plane, projecting out z). The returned array is oriented for
    ``imshow(origin="lower", extent=extent)`` — identical to ``paper_repro._hist2d``. ``img`` is
    caller-owned so it is allocated once (Taichi forbids new fields after the first kernel launch).
    """
    binw_a, binw_b = scatter_density(parts, img, extent, axes)
    return img.to_numpy().T / (binw_a * binw_b)  # mass/bin → surface density, imshow-oriented


__all__ = ["project_density", "scatter_density"]
