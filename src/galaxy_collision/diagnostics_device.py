"""Device-resident conservation diagnostics (Stage 5 / RV6b).

The Stage-3 ``diagnostics`` module computes energy/momentum/half-mass on the **host** in
fp64, fed by ``field.to_numpy()`` at the output cadence. That is correct and clear, but at
10–100M particles each cadence sample copies the full position/velocity/grid fields to the
host (~2.5 GB at 100M). This module provides the same scalars as **device reductions** so a
history-only sample copies nothing larger than a small histogram — the host ``diagnostics``
functions stay as the numeric reference these kernels are unit-tested against
(``tests/test_diagnostics.py``).

Accumulators are **f64** on CPU/CUDA. **Metal has no hardware fp64** (Stage 6 / RV15, D22),
so there the same scalars are reduced in a **Kahan-compensated fp32** path
(``*_kahan`` kernels, selected via :func:`galaxy_collision.backend.supports_fp64`): partials
accumulated in fp32 with a running compensation term, finished in host fp64. The radial
histogram stays fp32 either way — each particle mass (~1e4 M☉) is far above the ULP of a
bin total (~1e8), so the scatter sum is accurate without compensation.

This module deliberately omits ``from __future__ import annotations`` — PEP 563 would
stringize the ``ti.template()`` kernel annotations (see ``deposit.py``).
"""

import numpy as np
import taichi as ti

from galaxy_collision.backend import supports_fp64

# Lane count for the grid-stride Kahan particle reduction (Metal). Each lane Kahan-sums
# ~N/lanes particles sequentially; partials (11 × lanes fp32) finish in host fp64. Grid
# reductions instead parallelize over the (i,j) column (no fixed lane count needed).
_KAHAN_LANES = 8192

# Particle scalar-reduction buffer layout (all f64):
#   0      : M               (total mass)
#   1      : KE              (½ Σ m |v|²)
#   2,3,4  : P               (Σ m v — linear momentum)
#   5,6,7  : L               (Σ m r×v — angular momentum about the ORIGIN)
#   8,9,10 : Σ m r           (for the center of mass)
_PSCALARS = 11


@ti.kernel
def _reduce_particle_scalars(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    vel_x: ti.template(),
    vel_y: ti.template(),
    vel_z: ti.template(),
    mass: ti.template(),
    n: ti.i32,
    out: ti.template(),
):
    """One pass: total mass, KE, linear momentum, angular momentum (about origin), Σm·r."""
    for c in range(_PSCALARS):
        out[c] = 0.0
    for i in range(n):
        m = ti.cast(mass[i], ti.f64)
        x = ti.cast(pos_x[i], ti.f64)
        y = ti.cast(pos_y[i], ti.f64)
        z = ti.cast(pos_z[i], ti.f64)
        vx = ti.cast(vel_x[i], ti.f64)
        vy = ti.cast(vel_y[i], ti.f64)
        vz = ti.cast(vel_z[i], ti.f64)
        out[0] += m
        out[1] += 0.5 * m * (vx * vx + vy * vy + vz * vz)
        out[2] += m * vx
        out[3] += m * vy
        out[4] += m * vz
        out[5] += m * (y * vz - z * vy)
        out[6] += m * (z * vx - x * vz)
        out[7] += m * (x * vy - y * vx)
        out[8] += m * x
        out[9] += m * y
        out[10] += m * z


@ti.kernel
def _reduce_grid_pe(rho: ti.template(), phi: ti.template(), n: ti.i32, out: ti.template()):
    """Raw Σ ρΦ as a single f64 reduction into ``out[0]``.

    The ½·dx³ factor of the grid PE ½ΣρΦ·dx³ is applied host-side in :meth:`sample`, so this
    kernel is purely the reduction (no trailing scalar op to reason about). f64 accumulator —
    see the module docstring; summation order is non-deterministic across launches/backends, so
    it is exact only to fp64 round-off (immaterial under the run tolerances)."""
    out[0] = 0.0
    for i, j, k in ti.ndrange(n, n, n):
        out[0] += ti.cast(rho[i, j, k], ti.f64) * ti.cast(phi[i, j, k], ti.f64)


@ti.kernel
def _reduce_particle_scalars_kahan(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    vel_x: ti.template(),
    vel_y: ti.template(),
    vel_z: ti.template(),
    mass: ti.template(),
    n: ti.i32,
    nlanes: ti.i32,
    part: ti.template(),
):
    """Kahan-fp32 counterpart of :func:`_reduce_particle_scalars` for Metal (RV15/D22).

    Grid-stride over the 1-D particle index (``nlanes`` lanes, division-free): each lane
    Kahan-sums its slice of particles into an 11-vector of fp32 partials at ``part[:, lane]``;
    the host finishes the ``nlanes→1`` sum in fp64. The 11 components match the f64 kernel's
    buffer layout (M, KE, P, L, Σm·r)."""
    for lane in range(nlanes):
        s = ti.Vector.zero(ti.f32, _PSCALARS)
        c = ti.Vector.zero(ti.f32, _PSCALARS)
        idx = lane
        while idx < n:
            m = mass[idx]
            x, y, z = pos_x[idx], pos_y[idx], pos_z[idx]
            vx, vy, vz = vel_x[idx], vel_y[idx], vel_z[idx]
            v = ti.Vector([
                m,
                0.5 * m * (vx * vx + vy * vy + vz * vz),
                m * vx, m * vy, m * vz,
                m * (y * vz - z * vy), m * (z * vx - x * vz), m * (x * vy - y * vx),
                m * x, m * y, m * z,
            ])
            yk = v - c
            t = s + yk
            c = (t - s) - yk
            s = t
            idx += nlanes
        for cc in ti.static(range(_PSCALARS)):
            part[cc, lane] = s[cc]


@ti.kernel
def _reduce_grid_pe_kahan(rho: ti.template(), phi: ti.template(), n: ti.i32, part: ti.template()):
    """Kahan-fp32 counterpart of :func:`_reduce_grid_pe` for Metal (RV15/D22).

    Parallelizes over the ``(i,j)`` column; each lane Kahan-sums ``ρΦ`` along ``k`` into
    ``part[i,j]``; the host finishes the ``n²→1`` sum in fp64. The ½·dx³ factor is applied
    host-side, as in the f64 path."""
    for i, j in ti.ndrange(n, n):
        s = 0.0
        c = 0.0
        for k in range(n):
            v = rho[i, j, k] * phi[i, j, k]
            yk = v - c
            t = s + yk
            c = (t - s) - yk
            s = t
        part[i, j] = s


@ti.kernel
def _radial_mass_histogram(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    mass: ti.template(),
    n: ti.i32,
    cx: ti.f32,
    cy: ti.f32,
    cz: ti.f32,
    inv_binw: ti.f32,
    nbins: ti.i32,
    hist: ti.template(),
):
    """Mass-weighted histogram of particle radius about the center of mass ``(cx,cy,cz)``.

    Bin ``b`` collects particles with ``floor(r/binw) == b``; the last bin is a catch-all for
    ``r`` beyond the range, so all mass is accounted for. The host turns the cumulative
    histogram into the half-mass radius by linear interpolation at the half-mass crossing.
    """
    for b in range(nbins):
        hist[b] = 0.0
    for i in range(n):
        rx = ti.cast(pos_x[i], ti.f64) - cx
        ry = ti.cast(pos_y[i], ti.f64) - cy
        rz = ti.cast(pos_z[i], ti.f64) - cz
        r = ti.sqrt(rx * rx + ry * ry + rz * rz)
        b = ti.cast(r * inv_binw, ti.i32)
        if b < 0:
            b = 0
        if b >= nbins:
            b = nbins - 1
        hist[b] += ti.cast(mass[i], ti.f64)


@ti.kernel
def _radial_mass_histogram_f32(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    mass: ti.template(),
    n: ti.i32,
    cx: ti.f32,
    cy: ti.f32,
    cz: ti.f32,
    inv_binw: ti.f32,
    nbins: ti.i32,
    hist: ti.template(),
):
    """fp32 counterpart of :func:`_radial_mass_histogram` for Metal (no hardware fp64).

    Identical binning; the radius and the mass scatter are computed in fp32. The scatter sum
    is accurate without compensation (each mass ≫ the ULP of its bin total — see the module
    docstring), and the half-mass radius is binning/interpolation-limited regardless of the
    accumulator dtype (test tolerance 5e-3)."""
    for b in range(nbins):
        hist[b] = 0.0
    for i in range(n):
        rx = pos_x[i] - cx
        ry = pos_y[i] - cy
        rz = pos_z[i] - cz
        r = ti.sqrt(rx * rx + ry * ry + rz * rz)
        b = ti.cast(r * inv_binw, ti.i32)
        if b < 0:
            b = 0
        if b >= nbins:
            b = nbins - 1
        hist[b] += mass[i]


class DeviceDiagnostics:
    """Reusable scratch buffers + reductions for device-side conservation diagnostics.

    Allocate once per run (the buffers are tiny) and call :meth:`sample` each history
    cadence. Returns a dict with the same keys/units as the host ``measure`` path.
    """

    def __init__(self, grid_size: int, *, nbins: int = 2048, r_max: float | None = None):
        self.grid_size = grid_size
        self.nbins = nbins
        # Cover the whole box (+ a margin for particles that drift outside) so every particle
        # lands in a real bin; the box diagonal is the natural cap.
        self.r_max = float(r_max) if r_max is not None else float(grid_size) * np.sqrt(3.0)
        # Bin width sets the half-mass resolution (linear interpolation gives sub-bin accuracy).
        # At 256³/2048 bins this is ~0.22 kpc — fine for galaxy/Plummer half-mass radii of tens
        # of kpc, but a sub-kpc-half-mass object would need more bins (or a smaller r_max). Unlike
        # the host `lagrangian_radius`, this trades a hair of exactness (binning) for residency.
        self.binw = self.r_max / nbins
        # f64 reductions on CPU/CUDA; Kahan-fp32 partials (finished in host fp64) on Metal,
        # which has no hardware fp64 (Stage 6 / RV15, D22).
        self._fp64 = supports_fp64()
        if self._fp64:
            self._pbuf = ti.field(ti.f64, shape=_PSCALARS)
            self._gbuf = ti.field(ti.f64, shape=1)
            self._hist = ti.field(ti.f64, shape=nbins)
        else:
            self._ppart = ti.field(ti.f32, shape=(_PSCALARS, _KAHAN_LANES))
            self._gpart = ti.field(ti.f32, shape=(grid_size, grid_size))
            self._hist = ti.field(ti.f32, shape=nbins)

    def _half_mass_radius(self, parts, cx: float, cy: float, cz: float, total: float) -> float:
        hist_kernel = _radial_mass_histogram if self._fp64 else _radial_mass_histogram_f32
        hist_kernel(
            parts.pos_x, parts.pos_y, parts.pos_z, parts.mass, parts.n,
            cx, cy, cz, 1.0 / self.binw, self.nbins, self._hist,
        )
        hist = self._hist.to_numpy().astype(np.float64)  # f32 on Metal → f64 for the cumsum
        cum = np.cumsum(hist)
        target = 0.5 * total
        b = int(np.searchsorted(cum, target))
        if b >= self.nbins:
            return self.r_max
        # Linear interpolation within the crossing bin for sub-bin accuracy: the cumulative
        # mass rises from cum[b-1] at radius b·binw to cum[b] at (b+1)·binw.
        cum_lo = cum[b - 1] if b > 0 else 0.0
        frac = (target - cum_lo) / hist[b] if hist[b] > 0.0 else 0.0
        return float((b + frac) * self.binw)

    def sample(self, parts, grid, dx: float, *, half_mass: bool = True) -> dict:
        """Compute the scalar conservation diagnostics on the device.

        Returns ``M``, ``kinetic``, ``potential`` (grid PE), ``energy`` (KE+PE), ``momentum``
        and ``ang_momentum`` (length-3 fp64 vectors), the center of mass, and (when
        ``half_mass``) the half-mass radius. No full position/velocity/grid array is copied
        to the host — only the small particle buffer, the grid buffer, and the radial
        histogram (f64 reductions on CPU/CUDA; Kahan-fp32 partials finished in host fp64 on
        Metal — RV15/D22).
        """
        if self._fp64:
            _reduce_particle_scalars(
                parts.pos_x, parts.pos_y, parts.pos_z,
                parts.vel_x, parts.vel_y, parts.vel_z, parts.mass, parts.n, self._pbuf,
            )
            s = self._pbuf.to_numpy()
            _reduce_grid_pe(grid.rho, grid.phi, self.grid_size, self._gbuf)
            raw_pe = float(self._gbuf.to_numpy()[0])
        else:
            _reduce_particle_scalars_kahan(
                parts.pos_x, parts.pos_y, parts.pos_z,
                parts.vel_x, parts.vel_y, parts.vel_z, parts.mass, parts.n,
                _KAHAN_LANES, self._ppart,
            )
            s = self._ppart.to_numpy().astype(np.float64).sum(axis=1)  # (11, lanes) → (11,)
            _reduce_grid_pe_kahan(grid.rho, grid.phi, self.grid_size, self._gpart)
            raw_pe = float(self._gpart.to_numpy().astype(np.float64).sum())
        # Apply the ½·dx³ factor here (host-side) — the kernel returns the raw Σ ρΦ.
        pe = 0.5 * float(dx) ** 3 * raw_pe

        total = float(s[0])
        ke = float(s[1])
        momentum = s[2:5].copy()
        ang_momentum = s[5:8].copy()
        com = (s[8:11] / total) if total > 0.0 else np.zeros(3)

        out = {
            "M": total,
            "kinetic": ke,
            "potential": pe,
            "energy": ke + pe,
            "momentum": momentum,
            "ang_momentum": ang_momentum,
            "center_of_mass": com,
        }
        if half_mass:
            out["half_mass_radius"] = self._half_mass_radius(
                parts, float(com[0]), float(com[1]), float(com[2]), total
            )
        return out


__all__ = ["DeviceDiagnostics"]
