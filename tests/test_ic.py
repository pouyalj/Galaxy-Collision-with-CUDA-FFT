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
    COMPONENT_BH,
    COMPONENT_DISK,
    V_APPROACH,
    GalaxyModel,
    _circular_velocity,
    _disk_tables,
    _sample_disk,
    build_ic,
    select_tracer,
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


# --- Sun-like tracer (Stage 4 / 4B) ---------------------------------------------


def test_select_tracer_is_disk_particle_near_target_radius():
    ic = build_ic(_cfg(n_particles=20000))
    idx = select_tracer(ic, target_radius=8.32, gid=0)
    assert ic.gid[idx] == 0
    assert ic.component[idx] == COMPONENT_DISK
    # Its in-disk galactocentric radius about galaxy 0's BH is close to the target. With ~9000
    # disk particles spanning 25 kpc, the nearest one sits well within 1 kpc of 8.32.
    bh = (ic.gid == 0) & (ic.component == COMPONENT_BH)
    c = ic.pos[bh][0]
    r = np.hypot(ic.pos[idx, 0] - c[0], ic.pos[idx, 1] - c[1])
    assert abs(r - 8.32) < 1.0, r


def test_select_tracer_requires_component_tag():
    from galaxy_collision.ic import ICResult

    bare = ICResult(pos=np.zeros((3, 3)), vel=np.zeros((3, 3)), mass=np.ones(3),
                    gid=np.zeros(3, np.int32), preset="plummer", particle_mass=1.0)
    with pytest.raises(ValueError, match="component"):
        select_tracer(bare)


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


def test_plummer_ic_mass_and_count():
    from galaxy_collision.ic import PlummerModel, build_plummer_ic

    cfg = SimConfig(ic_preset="plummer", grid_size=64, n_particles=8000, seed=1)
    model = PlummerModel(total_mass=1.0e10, scale_a=5.0, center=(32.0, 32.0, 32.0))
    ic = build_plummer_ic(cfg, model=model)
    assert ic.n == 8000
    assert ic.preset == "plummer"
    assert ic.mass.sum() == pytest.approx(1.0e10, rel=1e-5)
    assert ic.particle_mass == pytest.approx(1.0e10 / 8000)
    assert set(np.unique(ic.gid)) == {0}  # single galaxy


def test_plummer_ic_half_mass_radius():
    """Sampled half-mass radius matches the Plummer value r_½ ≈ 1.305 a."""
    from galaxy_collision import diagnostics
    from galaxy_collision.ic import PlummerModel, build_plummer_ic

    cfg = SimConfig(ic_preset="plummer", grid_size=128, n_particles=40_000, seed=2)
    a = 5.0
    model = PlummerModel(total_mass=1.0e10, scale_a=a, center=(64.0, 64.0, 64.0))
    ic = build_plummer_ic(cfg, model=model)
    r_half = diagnostics.lagrangian_radius(ic.mass, ic.pos, 0.5, center=np.array(model.center))
    assert r_half == pytest.approx(1.305 * a, rel=0.05)


def test_plummer_ic_is_virialized():
    """Sampled velocities satisfy the virial theorem 2T + W ≈ 0 (T/|W| ≈ ½)."""
    from galaxy_collision import diagnostics
    from galaxy_collision.ic import PlummerModel, build_plummer_ic

    cfg = SimConfig(ic_preset="plummer", grid_size=64, n_particles=3000, seed=5)
    model = PlummerModel(total_mass=1.0e10, scale_a=5.0, center=(32.0, 32.0, 32.0))
    ic = build_plummer_ic(cfg, model=model)
    t = diagnostics.kinetic_energy(ic.mass, ic.vel)
    # Direct softened PE with a tiny softening (≪ a) ≈ the analytic Plummer self-energy.
    w = diagnostics.potential_energy_direct(ic.mass, ic.pos, softening=0.05)
    assert t / abs(w) == pytest.approx(0.5, abs=0.08)


def test_plummer_ic_reproducible():
    from galaxy_collision.ic import build_plummer_ic

    cfg = SimConfig(ic_preset="plummer", grid_size=64, n_particles=2000, seed=11)
    a = build_plummer_ic(cfg)
    b = build_plummer_ic(cfg)
    np.testing.assert_array_equal(a.pos, b.pos)
    np.testing.assert_array_equal(a.vel, b.vel)


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
