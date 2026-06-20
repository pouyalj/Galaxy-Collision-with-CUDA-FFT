"""Cloud-In-Cell (CIC) mass deposition and force gather (Stage 3).

AGENT.md §5.5: trilinear (CIC) weights on *both* deposit and gather. The same
symmetric weights on each side make the scheme momentum-conserving and kill the
blocky "square galaxy" nearest-grid-point artifact the paper noted (legacy bug #5).

**Grid convention (node-centered — the single source of truth for Stage 3).**
Samples live at grid *nodes*, not cell centers. Node ``(i, j, k)`` sits at physical
position ``(i, j, k)·dx`` (``dx`` = node spacing in kpc, 1.0 for the box here), so the
``N`` nodes per axis span ``[0, (N−1)·dx]``. A particle at position ``x`` (kpc) sits at
grid coordinate ``x / dx`` and splits its mass across the 8 surrounding nodes by
trilinear weight; the half-integer point ``(i+0.5)·dx`` is the midpoint *between*
nodes (equidistant from all 8), which is why it deposits in equal eighths. ``rho`` is a
mass *density* (M_sun / kpc^3), so each deposited mass is divided by the cell volume
``dx^3``. **The FFT oracle, multigrid, and the analytic point-mass test must all adopt
this same node placement** — a half-cell disagreement here is the classic off-by-½
that corrupts Φ(r) ≈ −GM/r.

Side effect to keep in mind: with ``N = 256`` / ``dx = 1`` the nodes really span
``[0, 255]`` kpc, and the IC generator centers galaxies at node ``128 = N/2`` — so the
nominal "256 kpc box" is a half-cell loose. That is fine, but solver and tests must
agree on it.

**Open boundary conditions (D5): we never wrap.** Weight that would land outside
the box is dropped — a particle straddling an edge loses that fraction. Physical
mass lives well inside the 256 kpc box, so this is a guard against out-of-bounds
indexing (legacy bug #2), not a physical model. Gather treats out-of-box nodes as
zero acceleration, consistent with Φ → 0 far from the mass.
"""

import taichi as ti

# NOTE: this module deliberately does *not* use ``from __future__ import annotations``.
# PEP 563 would stringize the ``ti.template()`` kernel annotations, which Taichi parses
# at materialize time and cannot resolve from a string. Keep annotations as live objects.


@ti.kernel
def _deposit_cic(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    mass: ti.template(),
    rho: ti.template(),
    n: ti.i32,
    dx: ti.f32,
    grid_size: ti.i32,
):
    inv_vol = 1.0 / (dx * dx * dx)
    for idx in ti.grouped(rho):
        rho[idx] = 0.0
    for p in range(n):
        gx = pos_x[p] / dx
        gy = pos_y[p] / dx
        gz = pos_z[p] / dx
        ix = ti.cast(ti.floor(gx), ti.i32)
        iy = ti.cast(ti.floor(gy), ti.i32)
        iz = ti.cast(ti.floor(gz), ti.i32)
        fx = gx - ti.cast(ix, ti.f32)
        fy = gy - ti.cast(iy, ti.f32)
        fz = gz - ti.cast(iz, ti.f32)
        m = mass[p] * inv_vol
        for cx in ti.static(range(2)):
            for cy in ti.static(range(2)):
                for cz in ti.static(range(2)):
                    jx = ix + cx
                    jy = iy + cy
                    jz = iz + cz
                    if 0 <= jx < grid_size and 0 <= jy < grid_size and 0 <= jz < grid_size:
                        wx = fx if cx == 1 else 1.0 - fx
                        wy = fy if cy == 1 else 1.0 - fy
                        wz = fz if cz == 1 else 1.0 - fz
                        ti.atomic_add(rho[jx, jy, jz], m * wx * wy * wz)


@ti.kernel
def _gather_accel(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    ax_g: ti.template(),
    ay_g: ti.template(),
    az_g: ti.template(),
    acc_x: ti.template(),
    acc_y: ti.template(),
    acc_z: ti.template(),
    n: ti.i32,
    dx: ti.f32,
    grid_size: ti.i32,
):
    for p in range(n):
        gx = pos_x[p] / dx
        gy = pos_y[p] / dx
        gz = pos_z[p] / dx
        ix = ti.cast(ti.floor(gx), ti.i32)
        iy = ti.cast(ti.floor(gy), ti.i32)
        iz = ti.cast(ti.floor(gz), ti.i32)
        fx = gx - ti.cast(ix, ti.f32)
        fy = gy - ti.cast(iy, ti.f32)
        fz = gz - ti.cast(iz, ti.f32)
        ax = 0.0
        ay = 0.0
        az = 0.0
        for cx in ti.static(range(2)):
            for cy in ti.static(range(2)):
                for cz in ti.static(range(2)):
                    jx = ix + cx
                    jy = iy + cy
                    jz = iz + cz
                    if 0 <= jx < grid_size and 0 <= jy < grid_size and 0 <= jz < grid_size:
                        wx = fx if cx == 1 else 1.0 - fx
                        wy = fy if cy == 1 else 1.0 - fy
                        wz = fz if cz == 1 else 1.0 - fz
                        w = wx * wy * wz
                        ax += w * ax_g[jx, jy, jz]
                        ay += w * ay_g[jx, jy, jz]
                        az += w * az_g[jx, jy, jz]
        acc_x[p] = ax
        acc_y[p] = ay
        acc_z[p] = az


def deposit_density(parts, rho, dx: float = 1.0) -> None:
    """Deposit particle masses onto ``rho`` (M_sun/kpc^3) by CIC. Zeroes ``rho`` first."""
    grid_size = rho.shape[0]
    _deposit_cic(
        parts.pos_x, parts.pos_y, parts.pos_z, parts.mass, rho, parts.n, float(dx), grid_size
    )


def gather_acceleration(parts, ax_g, ay_g, az_g, acc_x, acc_y, acc_z, dx: float = 1.0) -> None:
    """CIC-interpolate the node acceleration field (``ax_g``..) to particles (``acc_x``..)."""
    grid_size = ax_g.shape[0]
    _gather_accel(
        parts.pos_x,
        parts.pos_y,
        parts.pos_z,
        ax_g,
        ay_g,
        az_g,
        acc_x,
        acc_y,
        acc_z,
        parts.n,
        float(dx),
        grid_size,
    )


__all__ = ["deposit_density", "gather_acceleration"]
