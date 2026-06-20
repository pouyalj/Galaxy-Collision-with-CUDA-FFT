"""Tests for CIC deposition and force gather (Stage 3, Checkpoint A).

The headline correctness properties:
- **Mass conservation:** total deposited mass (∑ρ·cell_vol) equals ∑particle mass
  for any interior particle distribution.
- **CIC split:** a particle on a node deposits to that one node; a particle at a
  cell center splits evenly across the 8 surrounding nodes.
- **Gather partition-of-unity:** a constant acceleration field gathers back exactly,
  and a linear field is interpolated exactly (CIC is linearly exact).

All need Taichi and self-skip without it.
"""

from __future__ import annotations

import numpy as np
import pytest


def _make_parts(pos, mass):
    """Allocate a ParticleState and load (N,3) positions + (N,) masses into it."""
    from galaxy_collision.data import ParticleState

    pos = np.asarray(pos, dtype=np.float32)
    mass = np.asarray(mass, dtype=np.float32)
    parts = ParticleState(pos.shape[0])
    parts.pos_x.from_numpy(np.ascontiguousarray(pos[:, 0]))
    parts.pos_y.from_numpy(np.ascontiguousarray(pos[:, 1]))
    parts.pos_z.from_numpy(np.ascontiguousarray(pos[:, 2]))
    parts.mass.from_numpy(mass)
    return parts


def test_mass_conservation():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import GridState
    from galaxy_collision.deposit import deposit_density

    ti.init(arch=ti.cpu, random_seed=0)
    n, gsize = 5000, 64
    rng = np.random.default_rng(0)
    # Keep particles strictly interior so no CIC weight is clipped at the boundary.
    pos = rng.uniform(2.0, gsize - 2.0, size=(n, 3))
    mass = rng.uniform(0.5, 2.0, size=n)
    parts = _make_parts(pos, mass)

    grid = GridState(gsize)
    deposit_density(parts, grid.rho, dx=1.0)
    total = grid.rho.to_numpy().sum() * 1.0**3  # cell volume = dx^3
    assert total == pytest.approx(mass.sum(), rel=1e-4)


def test_single_particle_on_node():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import GridState
    from galaxy_collision.deposit import deposit_density

    ti.init(arch=ti.cpu)
    gsize = 16
    parts = _make_parts([[8.0, 8.0, 8.0]], [3.0])
    grid = GridState(gsize)
    deposit_density(parts, grid.rho, dx=1.0)
    rho = grid.rho.to_numpy()
    assert rho[8, 8, 8] == pytest.approx(3.0)
    assert rho.sum() == pytest.approx(3.0)


def test_node_midpoint_splits_eight_ways():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import GridState
    from galaxy_collision.deposit import deposit_density

    ti.init(arch=ti.cpu)
    gsize = 16
    # (8.5, 8.5, 8.5) is the midpoint *between* nodes (node-centered convention),
    # equidistant from all 8 surrounding nodes -> equal eighths.
    parts = _make_parts([[8.5, 8.5, 8.5]], [8.0])
    grid = GridState(gsize)
    deposit_density(parts, grid.rho, dx=1.0)
    rho = grid.rho.to_numpy()
    block = rho[8:10, 8:10, 8:10]
    assert np.allclose(block, 1.0)  # 8 M_sun / 8 nodes
    assert rho.sum() == pytest.approx(8.0)


def test_open_boundary_guard_no_wrap():
    """Particles outside / straddling the box must not wrap or crash (legacy bug #2).

    Open BCs (D5): out-of-box CIC weight is dropped, never wrapped. We assert no
    exception, that deposited mass never exceeds the input mass, and that the face
    opposite a high-edge particle stays empty (a wrap would light it up).
    """
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import GridState
    from galaxy_collision.deposit import deposit_density

    ti.init(arch=ti.cpu)
    gsize = 16
    # Fully below the box, fully above it, and straddling the top x-face (15.5 puts
    # half the cloud on node 15 and half on the nonexistent node 16, which is dropped).
    parts = _make_parts(
        [[-3.0, 8.0, 8.0], [40.0, 8.0, 8.0], [15.5, 8.0, 8.0]],
        [1.0, 1.0, 1.0],
    )
    grid = GridState(gsize)
    deposit_density(parts, grid.rho, dx=1.0)  # must not raise
    rho = grid.rho.to_numpy()

    total = rho.sum()
    assert total <= 3.0 + 1e-5  # no mass invented
    # Only the straddling particle contributes, and only its in-box half (node 15).
    assert rho[15, 8, 8] == pytest.approx(0.5)
    # A wrap of the high-x particles would deposit on the x=0 face — it must be empty.
    assert rho[0, :, :].sum() == pytest.approx(0.0)


def test_gather_constant_field():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.deposit import gather_acceleration

    ti.init(arch=ti.cpu, random_seed=1)
    n, gsize = 1000, 32
    rng = np.random.default_rng(1)
    pos = rng.uniform(1.0, gsize - 2.0, size=(n, 3))
    parts = _make_parts(pos, np.ones(n))

    ax_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    ay_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    az_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    ax_g.fill(2.0)
    ay_g.fill(-3.0)
    az_g.fill(0.5)
    acc_x = ti.field(ti.f32, shape=n)
    acc_y = ti.field(ti.f32, shape=n)
    acc_z = ti.field(ti.f32, shape=n)

    gather_acceleration(parts, ax_g, ay_g, az_g, acc_x, acc_y, acc_z, dx=1.0)
    assert np.allclose(acc_x.to_numpy(), 2.0, atol=1e-5)
    assert np.allclose(acc_y.to_numpy(), -3.0, atol=1e-5)
    assert np.allclose(acc_z.to_numpy(), 0.5, atol=1e-5)


def test_gather_linear_field_is_exact():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.deposit import gather_acceleration

    ti.init(arch=ti.cpu, random_seed=2)
    n, gsize = 500, 32
    rng = np.random.default_rng(2)
    pos = rng.uniform(1.0, gsize - 2.0, size=(n, 3)).astype(np.float32)
    parts = _make_parts(pos, np.ones(n))

    # ax_g[i,j,k] = 0.3*i + 0.1*j - 0.2*k ; CIC interpolates linear fields exactly.
    ii, jj, kk = np.meshgrid(
        np.arange(gsize), np.arange(gsize), np.arange(gsize), indexing="ij"
    )
    field = (0.3 * ii + 0.1 * jj - 0.2 * kk).astype(np.float32)
    ax_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    ay_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    az_g = ti.field(ti.f32, shape=(gsize, gsize, gsize))
    ax_g.from_numpy(field)
    ay_g.fill(0.0)
    az_g.fill(0.0)
    acc_x = ti.field(ti.f32, shape=n)
    acc_y = ti.field(ti.f32, shape=n)
    acc_z = ti.field(ti.f32, shape=n)

    gather_acceleration(parts, ax_g, ay_g, az_g, acc_x, acc_y, acc_z, dx=1.0)
    expected = 0.3 * pos[:, 0] + 0.1 * pos[:, 1] - 0.2 * pos[:, 2]
    assert np.allclose(acc_x.to_numpy(), expected, atol=1e-3)
