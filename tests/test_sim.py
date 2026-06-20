"""Stage-3 exit gate: the full PM pipeline holds a Plummer sphere in equilibrium.

AGENT.md §6 test 3 / §7 exit criteria: a sampled Plummer sphere run through the real
kick→drift→deposit→solve→grad→gather→kick loop must stay in steady state — directly
testing the "density dissipates over time" bug the 2020 code had. Over **several dynamical
times** we check that the sphere neither collapses nor evaporates (half-mass radius stays
in a band), stays virialized (KE/|PE| ≈ ½), conserves energy under threshold, and barely
drifts in linear momentum.

We use a *dense* sphere (small dynamical time t_dyn ≈ 8.3 Myr) so several crossings fit in
a CI-affordable step count — running long enough that a real instability would actually
show. The grid-softened sphere relaxes mildly to the grid's equilibrium (r_half drifts a
few %, the outer Plummer tail beyond the box is dropped), which the bands accommodate; a
collapse/evaporation would blow well past them. The exact FFT oracle is one check, the
multigrid production path the other. (Two-body Kepler lives in tests/test_integrator.py.)
"""

from __future__ import annotations

import numpy as np
import pytest

from galaxy_collision import units
from galaxy_collision.config import SimConfig

# A dense Plummer sphere: t_dyn = sqrt(a^3/GM) ≈ 8.3 Myr, so ~5-7 t_dyn is ~45-120 steps.
_PLUMMER_MASS = 4.0e11
_PLUMMER_A = 5.0
# Stability bands (observed over 5-7 t_dyn: r_half ~1.0-1.11, virial ~0.48-0.52, fft drift
# ~5e-3, mg drift ~1e-2, momentum drift ~1e-4). Bands give margin but stay non-vacuous: a
# collapse drives r_half far below 0.8, evaporation far above 1.2.
R_HALF_LO, R_HALF_HI = 0.80, 1.20
VIRIAL_LO, VIRIAL_HI = 0.45, 0.55
EDRIFT_FFT, EDRIFT_MG = 0.01, 0.03
PDRIFT_TOL = 1.0e-3


def _dense_plummer(grid_size):
    from galaxy_collision.ic import PlummerModel

    c = grid_size / 2.0
    return PlummerModel(total_mass=_PLUMMER_MASS, scale_a=_PLUMMER_A, center=(c, c, c))


def _momentum_scale():
    return _PLUMMER_MASS * np.sqrt(units.G * _PLUMMER_MASS / _PLUMMER_A)


def _plummer_config(solver, **overrides):
    base = dict(
        name="plummer-gate", backend="cpu", ic_preset="plummer", solver=solver,
        grid_size=48, n_particles=2000, dt=0.5, steps=120, output_cadence=0, seed=7,
    )
    base.update(overrides)
    return SimConfig(**base)


def _metrics(history):
    rh0 = history[0]["half_mass_radius"]
    rh = np.array([h["half_mass_radius"] / rh0 for h in history])
    e0 = history[0]["energy"]
    edrift = np.array([abs((h["energy"] - e0) / e0) for h in history])
    virial = np.array([h["kinetic"] / abs(h["potential"]) for h in history])
    p0 = history[0]["momentum"]
    pdrift = np.array([np.linalg.norm(h["momentum"] - p0) for h in history]) / _momentum_scale()
    return rh, edrift, virial, pdrift


def _assert_stable(rh, virial, pdrift, n_tdyn):
    assert n_tdyn >= 5.0, f"run too short to test stability: {n_tdyn:.1f} t_dyn"
    assert rh.min() > R_HALF_LO and rh.max() < R_HALF_HI, (
        f"r_half left band: [{rh.min():.3f}, {rh.max():.3f}]"
    )
    assert VIRIAL_LO < virial.min() and virial.max() < VIRIAL_HI, "virial ratio left equilibrium"
    assert pdrift.max() < PDRIFT_TOL, f"momentum drift {pdrift.max():.2e} exceeds {PDRIFT_TOL}"


def test_plummer_stays_stable_fft():
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision import sim

    cfg = _plummer_config("fft", grid_size=48, n_particles=2000, dt=0.5, steps=120)
    res = sim.run_simulation(
        cfg, write_snapshots=False, history_cadence=8, plummer_model=_dense_plummer(48)
    )
    rh, edrift, virial, pdrift = _metrics(res["history"])
    n_tdyn = cfg.steps * cfg.dt / np.sqrt(_PLUMMER_A**3 / (units.G * _PLUMMER_MASS))
    _assert_stable(rh, virial, pdrift, n_tdyn)
    assert edrift.max() < EDRIFT_FFT, f"FFT energy drift {edrift.max():.2e} exceeds {EDRIFT_FFT}"


def test_plummer_stays_stable_multigrid():
    """The production solver (open-BC multigrid), warm-started, over >5 dynamical times."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision import sim

    cfg = _plummer_config("multigrid", grid_size=48, n_particles=1500, dt=1.0, steps=45)
    res = sim.run_simulation(
        cfg, write_snapshots=False, history_cadence=8,
        solver_kwargs={"n_cycles": 10}, plummer_model=_dense_plummer(48),
    )
    rh, edrift, virial, pdrift = _metrics(res["history"])
    n_tdyn = cfg.steps * cfg.dt / np.sqrt(_PLUMMER_A**3 / (units.G * _PLUMMER_MASS))
    _assert_stable(rh, virial, pdrift, n_tdyn)
    assert edrift.max() < EDRIFT_MG, f"mg energy drift {edrift.max():.2e} exceeds {EDRIFT_MG}"


def test_snapshots_written_and_readable(tmp_path):
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision import sim
    from galaxy_collision.ic import PlummerModel
    from galaxy_collision.io import read_snapshot

    # A compact sphere (a=2) that fits a small 32^3 box; this test checks I/O, not physics.
    # output_cadence=2 with steps=5 -> snapshots at 0, 2, 4, and the forced final step 5.
    cfg = _plummer_config(
        "fft", grid_size=32, n_particles=500, steps=5, output_cadence=2,
        output_dir=str(tmp_path), seed=3,
    )
    compact = PlummerModel(total_mass=1.0e10, scale_a=2.0, center=(16.0, 16.0, 16.0))
    res = sim.run_simulation(cfg, history_cadence=2, plummer_model=compact)
    assert len(res["snapshots"]) == 4  # steps 0, 2, 4, and the forced final (5)
    snap = read_snapshot(res["snapshots"][-1])
    assert snap.step == 5
    assert snap.n == 500
    assert snap.rho is not None and snap.phi is not None
    assert "energy" in snap.diagnostics  # diagnostics persisted with the snapshot


def test_run_rejects_ic_that_does_not_fit_box():
    """The box-fit guard refuses a two-galaxy IC whose centers fall outside a small grid."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision import sim

    cfg = SimConfig(ic_preset="two_galaxy_4v", grid_size=48, n_particles=1000, steps=1, seed=0)
    with pytest.raises(ValueError, match="outside the .* box"):
        sim.run_simulation(cfg, write_snapshots=False)


def test_build_ic_dispatch():
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision.sim import _build_ic

    plummer = _build_ic(_plummer_config("fft", grid_size=64, n_particles=1000))
    assert plummer.preset == "plummer"
    assert plummer.n == 1000

    two = _build_ic(SimConfig(ic_preset="two_galaxy_4v", grid_size=256, n_particles=1000))
    assert set(np.unique(two.gid)) == {0, 1}

    with pytest.raises(NotImplementedError, match="no Stage-3 simulation builder"):
        _build_ic(SimConfig(ic_preset="hello"))
