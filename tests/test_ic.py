"""Tests for two-galaxy initial-condition generation (Stage 2).

Exit-criteria focus: ICs match the target mass and disk profile, and are exactly
reproducible from a seed. Pure NumPy (no Taichi) except the load-into-field test.
"""

from __future__ import annotations

import numpy as np
import pytest

from galaxy_collision import units
from galaxy_collision.config import GALAXY_MASS_MSUN, SimConfig
from galaxy_collision.ic import (
    V_APPROACH,
    GalaxyModel,
    _circular_velocity,
    _disk_tables,
    _sample_disk,
    build_ic,
)


def _cfg(**kw) -> SimConfig:
    base = dict(ic_preset="two_galaxy_4v", n_particles=100_000, seed=0)
    base.update(kw)
    return SimConfig(**base)


# --- Reproducibility ------------------------------------------------------------


def test_reproducible_same_seed():
    a = build_ic(_cfg(seed=42))
    b = build_ic(_cfg(seed=42))
    assert np.array_equal(a.pos, b.pos)
    assert np.array_equal(a.vel, b.vel)
    assert np.array_equal(a.mass, b.mass)
    assert np.array_equal(a.gid, b.gid)


def test_different_seed_differs():
    a = build_ic(_cfg(seed=1))
    b = build_ic(_cfg(seed=2))
    assert not np.array_equal(a.pos, b.pos)


# --- Counts, mass, IDs ----------------------------------------------------------


def test_particle_count_includes_two_bhs():
    ic = build_ic(_cfg(n_particles=100_000))
    # 50k stars/galaxy + 1 BH each.
    assert ic.n == 100_000 + 2
    assert set(np.unique(ic.gid)) == {0, 1}


def test_total_mass_matches_target():
    ic = build_ic(_cfg(n_particles=100_000))
    for g in (0, 1):
        sel = ic.gid == g
        stars = sel & np.isclose(ic.mass, ic.particle_mass)
        assert np.isclose(ic.mass[stars].sum(), GALAXY_MASS_MSUN, rtol=1e-3)


def test_particle_mass_knob_derives_count():
    pm = 4.64e6
    ic = build_ic(_cfg(n_particles=None, particle_mass=pm))
    n_per_galaxy = round(GALAXY_MASS_MSUN / pm)
    assert ic.n == 2 * n_per_galaxy + 2
    assert np.isclose(ic.particle_mass, pm)


# --- Central black holes (D9) ---------------------------------------------------


def test_central_bh_present_and_massive():
    ic = build_ic(_cfg(n_particles=100_000))
    model = GalaxyModel()
    bh = np.isclose(ic.mass, model.bh_mass)
    assert bh.sum() == 2  # one per galaxy
    # Each BH sits at its galaxy center and carries the configured mass.
    box_center = 256 / 2.0
    centers = {0: box_center - 45.0, 1: box_center + 45.0}  # separation/2 = 45
    for g in (0, 1):
        idx = np.where(bh & (ic.gid == g))[0]
        assert idx.size == 1
        assert np.isclose(ic.pos[idx[0], 0], centers[g], atol=1e-3)
        assert ic.mass[idx[0]] > 10 * ic.particle_mass


# --- Disk spatial profile -------------------------------------------------------


def test_disk_radial_profile_matches_surface_density():
    model = GalaxyModel()
    r_tab, cdf_tab = _disk_tables(model)
    rng = np.random.default_rng(0)
    pos, _ = _sample_disk(rng, 200_000, model, disk_mass=2.0e11, bulge_mass=2.0e10)
    radius = np.hypot(pos[:, 0], pos[:, 1])
    # Empirical CDF must track the analytic Σ(r)·r CDF.
    for rq in (3.0, 6.0, 12.0, 20.0):
        empirical = float((radius <= rq).mean())
        analytic = float(np.interp(rq, r_tab, cdf_tab))
        assert abs(empirical - analytic) < 0.02, (rq, empirical, analytic)
    assert radius.max() <= model.disk_rmax + 1e-6


def test_disk_is_thin():
    model = GalaxyModel()
    rng = np.random.default_rng(0)
    pos, _ = _sample_disk(rng, 100_000, model, disk_mass=2.0e11, bulge_mass=2.0e10)
    # Laplace(0, z0) has std = sqrt(2)*z0.
    assert np.isclose(pos[:, 2].std(), np.sqrt(2.0) * model.disk_z0, rtol=0.1)


# --- Rotation curve -------------------------------------------------------------


def test_circular_velocity_sane_at_solar_radius():
    model = GalaxyModel()
    r_tab, cdf_tab = _disk_tables(model)
    disk_mass = model.disk_fraction * model.total_mass
    bulge_mass = (1 - model.disk_fraction) * model.total_mass
    vc = _circular_velocity(np.array([8.0]), model, disk_mass, bulge_mass, r_tab, cdf_tab)[0]
    vc_kms = units.kpc_per_myr_to_kms(vc)
    assert 150.0 < vc_kms < 280.0  # MW-like rotation speed


# --- Approach kinematics --------------------------------------------------------


def test_galaxies_approach_each_other():
    ic = build_ic(_cfg(ic_preset="two_galaxy_4v", n_particles=100_000))
    vx0 = ic.vel[ic.gid == 0, 0].mean()
    vx1 = ic.vel[ic.gid == 1, 0].mean()
    assert vx0 > 0 and vx1 < 0  # closing
    closing = vx0 - vx1
    assert np.isclose(closing, 4.0 * V_APPROACH, rtol=0.05)


def test_4v_faster_than_2v():
    fast = build_ic(_cfg(ic_preset="two_galaxy_4v"))
    slow = build_ic(_cfg(ic_preset="two_galaxy_2v"))
    cl_fast = fast.vel[fast.gid == 0, 0].mean() - fast.vel[fast.gid == 1, 0].mean()
    cl_slow = slow.vel[slow.gid == 0, 0].mean() - slow.vel[slow.gid == 1, 0].mean()
    assert np.isclose(cl_fast / cl_slow, 2.0, rtol=0.05)


# --- Output shape & guards ------------------------------------------------------


def test_output_dtypes_and_shapes():
    ic = build_ic(_cfg())
    assert ic.pos.dtype == np.float32 and ic.pos.shape == (ic.n, 3)
    assert ic.vel.dtype == np.float32 and ic.vel.shape == (ic.n, 3)
    assert ic.mass.dtype == np.float32 and ic.mass.shape == (ic.n,)
    assert "particles" in ic.summary()


@pytest.mark.parametrize("preset", ["hello", "plummer"])
def test_unimplemented_presets_raise(preset):
    with pytest.raises(NotImplementedError):
        build_ic(_cfg(ic_preset=preset))


# --- Taichi integration ---------------------------------------------------------


def test_load_into_particle_state():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import ParticleState
    from galaxy_collision.ic import load_into_particle_state

    ti.init(arch=ti.cpu)
    ic = build_ic(_cfg(n_particles=10_000))
    parts = ParticleState(ic.n)
    load_into_particle_state(ic, parts)
    ti.sync()
    assert np.isclose(parts.pos_x[0], ic.pos[0, 0], atol=1e-3)
    assert set(np.unique(parts.gid.to_numpy())) == {0, 1}
