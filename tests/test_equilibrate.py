"""Disk-velocity equilibration on the PM grid (§5.5 / D18 / RV-softening).

The analytic IC rotation curve (`_circular_velocity`) carries two approximations that the
grid does not: a 0.3 kpc Plummer softening (vs the ~1 kpc grid/CIC softening) and a disk
treated as *spherically* enclosed mass (vs the real razor-thin geometry). A disk launched on
that curve is out of balance with the force the PM chain actually applies, so it breathes/heats.

`equilibrate_disk_velocities` instead launches the disk on the curve the grid *measures*. These
tests pin the two corrections it makes — fixture-independent in direction:

1. **Inner disk:** the grid softens more than 0.3 kpc, so the measured v_c is *lower* than the
   analytic one inside a few kpc (no spurious over-rotation).
2. **Disk edge:** a thin disk's in-plane force exceeds the spherical-enclosed-mass estimate, so
   the measured v_c is *higher* than analytic near the edge.

and the consequence — over a short time (≪ the cold-disk numerical-heating timescale) the
equilibrated disk holds its half-mass radius far better than the analytic one. (Long-term
cold-disk spreading is a separate resolution/Toomre effect; warm disks are a Stage-8 item.)
"""

from __future__ import annotations

import numpy as np
import pytest

from galaxy_collision import units

_GS = 64
_C = _GS / 2.0


def _compact_galaxy(seed=3, n=12000):
    """One compact galaxy centered in a 64^3 box, so the inner-disk softening mismatch
    (0.3 vs ~1 kpc) actually bites at dx=1. Returns (model, pos, vel_analytic, mass, component)."""
    from galaxy_collision import ic as icmod

    model = icmod.GalaxyModel(
        total_mass=5.0e10,
        disk_fraction=0.9,
        disk_rmax=12.0,
        sigma_rb=1.0,
        sigma_lc=1.5,
        bulge_a=0.6,
        disk_z0=0.3,
        bh_mass=5.0e7,
        softening=0.3,
        center=(_C, _C, _C),
        bulk_velocity=(0.0, 0.0, 0.0),
    )
    rng = np.random.default_rng(seed)
    pos, vel, mass, gid, comp = icmod._build_galaxy(rng, 0, model, 5.0e10 / n, n)
    return model, pos, vel, mass, comp


def _v_phi(pos, vel, comp):
    """In-plane circular speed and cylindrical radius for the disk particles (about the center)."""
    from galaxy_collision import ic as icmod

    disk = comp == icmod.COMPONENT_DISK
    dx_, dy_ = pos[disk, 0] - _C, pos[disk, 1] - _C
    r = np.hypot(dx_, dy_)
    safe = np.maximum(r, 1e-6)
    ux, uy = dx_ / safe, dy_ / safe
    vphi = -vel[disk, 0] * uy + vel[disk, 1] * ux
    return r, vphi


def _equilibrate(pos, vel, mass, comp):
    """Return a velocity array equilibrated on the multigrid grid force."""
    import taichi as ti

    from galaxy_collision import ic as icmod
    from galaxy_collision.data import GridState
    from galaxy_collision.solver import make_solver

    icr = icmod.ICResult(
        pos=pos.astype(np.float32),
        vel=vel.astype(np.float32),
        mass=mass.astype(np.float32),
        gid=np.zeros(pos.shape[0], np.int32),
        preset="two_galaxy_4v",
        particle_mass=float(mass[comp != icmod.COMPONENT_BH][0]),
        component=comp.astype(np.int8),
    )
    grid = GridState(_GS)
    ax = ti.field(ti.f32, shape=(_GS,) * 3)
    ay = ti.field(ti.f32, shape=(_GS,) * 3)
    az = ti.field(ti.f32, shape=(_GS,) * 3)
    solver = make_solver("multigrid", _GS, dx=1.0, grav_constant=units.G, n_cycles=30)
    icmod.equilibrate_disk_velocities(
        icr, dx=1.0, solver=solver, rho=grid.rho, phi=grid.phi, ax=ax, ay=ay, az=az
    )
    return icr.vel.astype(np.float64)


def test_equilibration_corrects_inner_and_edge_rotation():
    """Mechanism (deterministic): measured v_c is lower than analytic inside, higher at the edge."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    ti.init(arch=ti.cpu)
    _, pos, vel_a, mass, comp = _compact_galaxy()
    vel_e = _equilibrate(pos, vel_a, mass, comp)

    r, vpa = _v_phi(pos, vel_a, comp)
    _, vpe = _v_phi(pos, vel_e, comp)

    inner = r < 2.5
    edge = (r > 9.0) & (r < 11.5)
    inner_ratio = np.median(vpe[inner]) / np.median(vpa[inner])
    edge_ratio = np.median(vpe[edge]) / np.median(vpa[edge])
    assert inner_ratio < 0.93, f"inner v_c not softened vs analytic: ratio {inner_ratio:.3f}"
    assert edge_ratio > 1.05, f"edge v_c missing the thin-disk enhancement: ratio {edge_ratio:.3f}"
    # The measured curve is physical: finite and positive across the disk body.
    assert np.all(np.isfinite(vpe)) and np.median(vpe[r < 11.0]) > 0.0


def test_equilibrated_disk_breathes_less_short_time():
    """Consequence: over ~50 Myr the equilibrated disk holds its half-mass radius and heats less
    radially than the analytic-softened launch (which is out of balance and breathes). Each case
    equilibrates and evolves with the *same* solver, mirroring run_simulation."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision import ic as icmod
    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.deposit import deposit_density, gather_acceleration, potential_to_accel
    from galaxy_collision.integrator import kdk_step
    from galaxy_collision.solver import make_solver

    ti.init(arch=ti.cpu)
    _, pos, vel_a, mass, comp = _compact_galaxy(n=8000)
    disk = comp == icmod.COMPONENT_DISK
    n = pos.shape[0]

    def evolve(equilibrate, steps=80, dt=0.6):
        parts = ParticleState(n)
        pm = float(mass[comp != icmod.COMPONENT_BH][0])
        icr = icmod.ICResult(
            pos=pos.astype(np.float32),
            vel=vel_a.astype(np.float32),
            mass=mass.astype(np.float32),
            gid=np.zeros(n, np.int32),
            preset="two_galaxy_4v",
            particle_mass=pm,
            component=comp.astype(np.int8),
        )
        grid = GridState(_GS)
        ax = ti.field(ti.f32, shape=(_GS,) * 3)
        ay = ti.field(ti.f32, shape=(_GS,) * 3)
        az = ti.field(ti.f32, shape=(_GS,) * 3)
        acc_x = ti.field(ti.f32, shape=n)
        acc_y = ti.field(ti.f32, shape=n)
        acc_z = ti.field(ti.f32, shape=n)
        # Warm-started, so a modest cycle count converges each step; the same solver does the
        # equilibration so the launch is in balance with exactly the force the loop applies.
        solver = make_solver("multigrid", _GS, dx=1.0, grav_constant=units.G, n_cycles=12)
        if equilibrate:
            icmod.equilibrate_disk_velocities(
                icr,
                dx=1.0,
                solver=solver,
                rho=grid.rho,
                phi=grid.phi,
                ax=ax,
                ay=ay,
                az=az,
            )
        icmod.load_into_particle_state(icr, parts)
        warm = {"on": False}

        def accel_fn():
            deposit_density(parts, grid.rho, 1.0)
            solver.solve(grid.rho, grid.phi, warm_start=warm["on"])
            potential_to_accel(grid.phi, ax, ay, az, 1.0)
            gather_acceleration(parts, ax, ay, az, acc_x, acc_y, acc_z, 1.0)
            warm["on"] = True

        def disk_stats():
            px, py = parts.pos_x.to_numpy(), parts.pos_y.to_numpy()
            vx, vy = parts.vel_x.to_numpy(), parts.vel_y.to_numpy()
            r = np.hypot(px[disk] - _C, py[disk] - _C)
            safe = np.maximum(r, 1e-6)
            v_r = vx[disk] * (px[disk] - _C) / safe + vy[disk] * (py[disk] - _C) / safe
            return float(np.median(r)), float(np.sqrt(np.mean(v_r**2)))

        accel_fn()
        rh0, _ = disk_stats()
        for _ in range(steps):
            kdk_step(parts, acc_x, acc_y, acc_z, accel_fn, dt)
        rh1, vr1 = disk_stats()
        return abs(rh1 / rh0 - 1.0), vr1

    drift_a, heat_a = evolve(equilibrate=False)
    drift_e, heat_e = evolve(equilibrate=True)

    # Absolute: the equilibrated disk holds its half-radius short-term. The comparisons use a
    # margin (the benefit is large; bare strict '<' could flake on shared solver/shot noise). The
    # deterministic mechanism check is test_equilibration_corrects_inner_and_edge_rotation.
    assert drift_e < 0.04, f"equilibrated disk r_half drifted {drift_e:.1%} short-term"
    assert drift_e < 0.8 * drift_a, f"didn't cut breathing: {drift_e:.3f} vs {drift_a:.3f}"
    assert heat_e < 0.95 * heat_a, f"didn't reduce heating: {heat_e:.4f} vs {heat_a:.4f}"
