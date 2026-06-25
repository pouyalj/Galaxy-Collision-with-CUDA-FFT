"""Energy-metric relationships (AGENT.md §6, §11 RV11).

The two-galaxy paper runs are far too large for the O(N²) direct PE, so the grid PE is the
only potential-energy metric available there. These tests *lock* how the grid PE relates to
the physical (direct) PE, justifying the decision to report grid-energy **drift** (not a
self-energy-subtracted absolute) for big runs:

- The grid PE ½ΣρΦ·dV recovers the physical softened pair-sum PE (ε ≈ 1 cell) up to a small,
  *stable* discretization offset — ~2% for the FFT oracle, ~4% for multigrid. That the offset is
  small and seed-stable means absolute grid energy is meaningful (within the offset) and its
  drift is a valid conservation metric.
- The offset is the ½ΣρΦ estimator bias, not the CIC self-energy (empirically ~0.1%), so no
  self-energy subtraction is applied (RV11).
"""

from __future__ import annotations

import numpy as np
import pytest

from galaxy_collision import diagnostics, units


def _plummer_grid_and_direct(solver_name, seed):
    import taichi as ti

    from galaxy_collision.config import SimConfig
    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.deposit import deposit_density
    from galaxy_collision.ic import PlummerModel, build_plummer_ic, load_into_particle_state
    from galaxy_collision.solver import make_solver

    ti.init(arch=ti.cpu)
    gs, dx = 48, 1.0  # N kept modest: the reference direct PE below is O(N²)
    c = gs / 2.0
    cfg = SimConfig(ic_preset="plummer", grid_size=gs, n_particles=4000, seed=seed)
    model = PlummerModel(total_mass=2.0e10, scale_a=4.0, center=(c, c, c))
    icr = build_plummer_ic(cfg, model=model)
    parts = ParticleState(icr.n)
    load_into_particle_state(icr, parts)
    grid = GridState(gs)
    deposit_density(parts, grid.rho, dx)
    phi = ti.field(ti.f32, shape=(gs, gs, gs))
    solver = make_solver(
        solver_name,
        gs,
        dx=dx,
        grav_constant=units.G,
        **({"n_cycles": 40} if solver_name == "multigrid" else {}),
    )
    solver.solve(grid.rho, phi)
    pe_grid = diagnostics.potential_energy_grid(grid.rho.to_numpy(), phi.to_numpy(), dx)
    pe_direct = diagnostics.potential_energy_direct(
        icr.mass.astype(np.float64), icr.pos.astype(np.float64), softening=1.0, grav=units.G
    )
    return pe_grid, pe_direct


@pytest.mark.parametrize("solver_name,hi", [("fft", 1.08), ("multigrid", 1.10)])
def test_grid_pe_recovers_direct_pe_within_offset(solver_name, hi):
    """Grid PE matches the physical (ε≈1 cell) direct PE up to a small, stable offset."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    ratios = []
    for seed in (1, 2, 3):
        pe_grid, pe_direct = _plummer_grid_and_direct(solver_name, seed)
        ratios.append(pe_grid / pe_direct)
    for r in ratios:
        assert 0.97 <= r < hi, f"{solver_name} grid/direct PE ratio {r:.3f} outside [0.97, {hi})"
    # Stable across realizations (offset is a discretization property, not realization-fragile).
    assert max(ratios) - min(ratios) < 0.04, f"{solver_name} offset varies too much: {ratios}"


def test_total_energy_grid_and_virial():
    """The grid total-energy helper sums KE + grid PE; the virial helper is offset-insensitive."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.config import SimConfig
    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.deposit import deposit_density
    from galaxy_collision.ic import PlummerModel, build_plummer_ic, load_into_particle_state
    from galaxy_collision.solver import make_solver

    ti.init(arch=ti.cpu)
    gs, dx = 64, 1.0
    c = gs / 2.0
    cfg = SimConfig(ic_preset="plummer", grid_size=gs, n_particles=8000, seed=7)
    model = PlummerModel(total_mass=1.0e10, scale_a=4.0, center=(c, c, c))
    icr = build_plummer_ic(cfg, model=model)
    parts = ParticleState(icr.n)
    load_into_particle_state(icr, parts)
    grid = GridState(gs)
    deposit_density(parts, grid.rho, dx)
    phi = ti.field(ti.f32, shape=(gs, gs, gs))
    make_solver("fft", gs, dx=dx, grav_constant=units.G).solve(grid.rho, phi)

    m = icr.mass.astype(np.float64)
    vel = icr.vel.astype(np.float64)
    rho_np, phi_np = grid.rho.to_numpy(), phi.to_numpy()
    ke = diagnostics.kinetic_energy(m, vel)
    pe = diagnostics.potential_energy_grid(rho_np, phi_np, dx)
    assert diagnostics.total_energy_grid(m, vel, rho_np, phi_np, dx) == pytest.approx(ke + pe)
    # A sampled Plummer sphere is in virial equilibrium: T/|W| ≈ 0.5.
    assert diagnostics.virial_ratio(ke, pe) == pytest.approx(0.5, abs=0.1)


def test_device_diagnostics_match_host_reference():
    """RV6b: the device-resident reductions must reproduce the host fp64 diagnostics.

    KE / grid-PE / linear & angular momentum are exact reductions (agree to fp64 round-off);
    the half-mass radius comes from a radial histogram + linear interpolation, so it agrees
    with the exact sorted host value to well under a percent.
    """
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.config import SimConfig
    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.deposit import deposit_density
    from galaxy_collision.diagnostics_device import DeviceDiagnostics
    from galaxy_collision.ic import PlummerModel, build_plummer_ic, load_into_particle_state
    from galaxy_collision.solver import make_solver

    ti.init(arch=ti.cpu)
    gs, dx = 64, 1.0
    c = gs / 2.0
    cfg = SimConfig(ic_preset="plummer", grid_size=gs, n_particles=8000, seed=11)
    icr = build_plummer_ic(cfg, model=PlummerModel(total_mass=2.0e10, scale_a=4.0,
                                                   center=(c, c, c)))
    parts = ParticleState(icr.n)
    load_into_particle_state(icr, parts)
    grid = GridState(gs)
    deposit_density(parts, grid.rho, dx)
    make_solver("fft", gs, dx=dx, grav_constant=units.G).solve(grid.rho, grid.phi)

    dev = DeviceDiagnostics(gs).sample(parts, grid, dx)

    m = icr.mass.astype(np.float64)
    pos, vel = icr.pos.astype(np.float64), icr.vel.astype(np.float64)
    p_ref = diagnostics.linear_momentum(m, vel)
    l_ref = diagnostics.angular_momentum(m, pos, vel)

    # On fp64 backends (CPU/CUDA) the device reductions are exact to fp64 round-off. On Metal
    # there is no hardware fp64 (RV15/D22), so they run a Kahan-compensated fp32 path; measured
    # here KE ~3e-10, grid-PE ~4e-9 relative — near-fp64, but the 1e-9 bar is fp64-only.
    from galaxy_collision.backend import supports_fp64

    scalar_rel = 1e-9 if supports_fp64() else 1e-6
    assert dev["kinetic"] == pytest.approx(diagnostics.kinetic_energy(m, vel), rel=scalar_rel)
    assert dev["potential"] == pytest.approx(
        diagnostics.potential_energy_grid(grid.rho.to_numpy(), grid.phi.to_numpy(), dx),
        rel=scalar_rel,
    )
    if supports_fp64():
        # Vector reductions: device vs NumPy differ only by fp64 summation order (parallel
        # tree-reduction vs sequential), which shows up in near-zero, cancellation-heavy
        # components — so compare at 1e-6, not machine epsilon.
        np.testing.assert_allclose(dev["momentum"], p_ref, rtol=1e-6)
        np.testing.assert_allclose(dev["ang_momentum"], l_ref, rtol=1e-6)
    else:
        # Metal/fp32: the Plummer net momentum is ~0 (near-total cancellation), so rtol on the
        # vector is meaningless. Check the deviation against the *summation* scale (Σm|v| for P,
        # Σm|r×v| for L) — i.e. that the reduction reproduces the host to fp32 precision relative
        # to the magnitudes actually summed. Observed ~1e-9 of scale here; 1e-5 is safe headroom.
        p_scale = float((m[:, None] * np.abs(vel)).sum(axis=0).max())
        l_scale = float((m[:, None] * np.abs(np.cross(pos, vel))).sum(axis=0).max())
        assert np.abs(dev["momentum"] - p_ref).max() < 1e-5 * p_scale
        assert np.abs(dev["ang_momentum"] - l_ref).max() < 1e-5 * l_scale
    assert dev["half_mass_radius"] == pytest.approx(
        diagnostics.lagrangian_radius(m, pos, 0.5), rel=5e-3
    )
