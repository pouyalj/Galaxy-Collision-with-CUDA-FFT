"""KDK leapfrog + softened forces + diagnostics (Stage 3, Checkpoint C).

The headline is the **two-body Kepler orbit** (AGENT.md §6 test 4): a circular binary
integrated with KDK + direct softened gravity must conserve energy and angular momentum
and stay on its orbit for a full period — validating the integrator and softening.
Supporting tests cover kick/drift kinematics, the φ → −∇φ "grad" step, the direct force
vs the analytic two-body value, and the diagnostics helpers.
"""

from __future__ import annotations

import numpy as np
import pytest

from galaxy_collision import diagnostics, units


def _alloc_vec_fields(ti, n):
    f = lambda: ti.field(ti.f32, shape=n)  # noqa: E731
    return f(), f(), f(), f(), f(), f()


def test_kick_and_drift_uniform_acceleration():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.integrator import drift, kick

    ti.init(arch=ti.cpu)
    n = 4
    px, py, pz, vx, vy, vz = _alloc_vec_fields(ti, n)
    ax = ti.field(ti.f32, shape=n)
    ay = ti.field(ti.f32, shape=n)
    az = ti.field(ti.f32, shape=n)
    ax.fill(0.5)  # constant acceleration along x
    ay.fill(0.0)
    az.fill(0.0)

    # One KDK step from rest at origin with a = 0.5: kick½ -> drift -> kick½.
    dt = 0.2
    kick(vx, vy, vz, ax, ay, az, n, dt / 2)
    drift(px, py, pz, vx, vy, vz, n, dt)
    kick(vx, vy, vz, ax, ay, az, n, dt / 2)

    # After one KDK step under constant a from rest: x = ½ a dt², v = a dt.
    assert vx.to_numpy()[0] == pytest.approx(0.5 * dt, rel=1e-5)
    assert px.to_numpy()[0] == pytest.approx(0.5 * 0.5 * dt**2, rel=1e-5)


def test_potential_to_accel_central_difference():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.deposit import potential_to_accel

    ti.init(arch=ti.cpu)
    n, dx = 16, 1.0
    # phi = 2x + 3y - z  ->  g = -grad(phi) = (-2, -3, +1) everywhere on the interior.
    ii, jj, kk = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    phi = ti.field(ti.f32, shape=(n, n, n))
    phi.from_numpy((2.0 * ii + 3.0 * jj - 1.0 * kk).astype(np.float32))
    axg = ti.field(ti.f32, shape=(n, n, n))
    ayg = ti.field(ti.f32, shape=(n, n, n))
    azg = ti.field(ti.f32, shape=(n, n, n))
    potential_to_accel(phi, axg, ayg, azg, dx=dx)
    interior = (slice(1, n - 1),) * 3
    assert np.allclose(axg.to_numpy()[interior], -2.0, atol=1e-4)
    assert np.allclose(ayg.to_numpy()[interior], -3.0, atol=1e-4)
    assert np.allclose(azg.to_numpy()[interior], 1.0, atol=1e-4)


def test_direct_accel_matches_two_body_analytic():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.integrator import direct_accel

    ti.init(arch=ti.cpu)
    n = 2
    px, py, pz, vx, vy, vz = _alloc_vec_fields(ti, n)
    mass = ti.field(ti.f32, shape=n)
    acc_x = ti.field(ti.f32, shape=n)
    acc_y = ti.field(ti.f32, shape=n)
    acc_z = ti.field(ti.f32, shape=n)
    r = 10.0
    px.from_numpy(np.array([0.0, r], dtype=np.float32))
    m2 = 5.0e10
    mass.from_numpy(np.array([1.0e10, m2], dtype=np.float32))

    direct_accel(px, py, pz, mass, acc_x, acc_y, acc_z, n, units.G, 0.0)
    # Body 0 feels +x pull of magnitude G·m2/r² toward body 1.
    assert acc_x.to_numpy()[0] == pytest.approx(units.G * m2 / r**2, rel=1e-4)
    assert acc_x.to_numpy()[1] == pytest.approx(-units.G * 1.0e10 / r**2, rel=1e-4)


def test_two_body_kepler_orbit():
    """Circular binary over one period: conserve E and L, stay on orbit (test 4)."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.integrator import direct_accel, drift, kick

    ti.init(arch=ti.cpu)
    n = 2
    m = 1.0e10
    m_tot = 2.0 * m
    r = 10.0  # separation (kpc)
    eps = 0.01  # softening ≪ r, negligible
    v_rel = np.sqrt(units.G * m_tot / r)  # circular relative speed
    v_each = 0.5 * v_rel
    period = 2.0 * np.pi * np.sqrt(r**3 / (units.G * m_tot))

    px, py, pz, vx, vy, vz = _alloc_vec_fields(ti, n)
    mass = ti.field(ti.f32, shape=n)
    acc_x = ti.field(ti.f32, shape=n)
    acc_y = ti.field(ti.f32, shape=n)
    acc_z = ti.field(ti.f32, shape=n)
    mass.from_numpy(np.array([m, m], dtype=np.float32))
    px.from_numpy(np.array([-r / 2, r / 2], dtype=np.float32))
    vy.from_numpy(np.array([-v_each, v_each], dtype=np.float32))  # counter-rotating about CM

    dt = period / 2000.0
    n_steps = 2000  # one full period

    def accel():
        direct_accel(px, py, pz, mass, acc_x, acc_y, acc_z, n, units.G, eps**2)

    def state():
        pos = np.stack([px.to_numpy(), py.to_numpy(), pz.to_numpy()], axis=1).astype(np.float64)
        vel = np.stack([vx.to_numpy(), vy.to_numpy(), vz.to_numpy()], axis=1).astype(np.float64)
        return pos, vel

    m_np = np.array([m, m], dtype=np.float64)
    pos0, vel0 = state()
    e0 = diagnostics.total_energy_direct(m_np, pos0, vel0, softening=eps)
    l0 = diagnostics.angular_momentum(m_np, pos0, vel0)
    p0 = diagnostics.linear_momentum(m_np, vel0)

    sep_min, sep_max = r, r
    for _ in range(n_steps):
        accel()
        kick(vx, vy, vz, acc_x, acc_y, acc_z, n, dt / 2)
        drift(px, py, pz, vx, vy, vz, n, dt)
        accel()
        kick(vx, vy, vz, acc_x, acc_y, acc_z, n, dt / 2)
        pos, _ = state()
        sep = np.linalg.norm(pos[1] - pos[0])
        sep_min, sep_max = min(sep_min, sep), max(sep_max, sep)

    posf, velf = state()
    ef = diagnostics.total_energy_direct(m_np, posf, velf, softening=eps)
    lf = diagnostics.angular_momentum(m_np, posf, velf)
    pf = diagnostics.linear_momentum(m_np, velf)

    # Symplectic KDK: energy drift bounded and tiny; angular momentum ~exactly conserved.
    assert abs((ef - e0) / e0) < 1e-3
    assert np.linalg.norm(lf - l0) / np.linalg.norm(l0) < 1e-4
    assert np.linalg.norm(pf - p0) < 1e-6 * m_tot * v_rel + 1e-9
    # Orbit stays circular (separation barely varies) and returns near the start.
    assert (sep_max - sep_min) / r < 0.02
    assert np.linalg.norm(posf[0] - pos0[0]) / r < 0.02


def _two_body_state(ti, m, r, v_each):
    """A circular two-body binary loaded into a ParticleState (CM at origin)."""
    from galaxy_collision.data import ParticleState

    parts = ParticleState(2)
    parts.pos_x.from_numpy(np.array([-r / 2, r / 2], dtype=np.float32))
    parts.pos_y.from_numpy(np.zeros(2, dtype=np.float32))
    parts.pos_z.from_numpy(np.zeros(2, dtype=np.float32))
    parts.vel_x.from_numpy(np.zeros(2, dtype=np.float32))
    parts.vel_y.from_numpy(np.array([-v_each, v_each], dtype=np.float32))
    parts.vel_z.from_numpy(np.zeros(2, dtype=np.float32))
    parts.mass.from_numpy(np.array([m, m], dtype=np.float32))
    return parts


def _read_state(parts):
    pos = np.stack(
        [parts.pos_x.to_numpy(), parts.pos_y.to_numpy(), parts.pos_z.to_numpy()], axis=1
    ).astype(np.float64)
    vel = np.stack(
        [parts.vel_x.to_numpy(), parts.vel_y.to_numpy(), parts.vel_z.to_numpy()], axis=1
    ).astype(np.float64)
    return pos, vel


def test_kdk_step_one_eval_conserves_and_matches_inline():
    """kdk_step: exactly one force eval/step, conserves E, and == the naive 2-eval loop."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.integrator import direct_accel, drift, kdk_step, kick

    ti.init(arch=ti.cpu)
    m, r, eps = 1.0e10, 10.0, 0.01
    v_rel = np.sqrt(units.G * 2.0 * m / r)
    v_each = 0.5 * v_rel
    period = 2.0 * np.pi * np.sqrt(r**3 / (units.G * 2.0 * m))
    dt = period / 1000.0
    n_steps = 1000
    m_np = np.array([m, m], dtype=np.float64)

    # --- run via kdk_step (cached accel: one eval/step) ---
    parts = _two_body_state(ti, m, r, v_each)
    ax = ti.field(ti.f32, shape=2)
    ay = ti.field(ti.f32, shape=2)
    az = ti.field(ti.f32, shape=2)
    calls = {"n": 0}

    def accel_fn():
        calls["n"] += 1
        direct_accel(
            parts.pos_x, parts.pos_y, parts.pos_z, parts.mass, ax, ay, az, 2, units.G, eps**2
        )

    pos0, vel0 = _read_state(parts)
    e0 = diagnostics.total_energy_direct(m_np, pos0, vel0, softening=eps)
    accel_fn()  # prime acc once before the loop
    for _ in range(n_steps):
        kdk_step(parts, ax, ay, az, accel_fn, dt)
    pos_kdk, vel_kdk = _read_state(parts)
    ef = diagnostics.total_energy_direct(m_np, pos_kdk, vel_kdk, softening=eps)

    # One acceleration evaluation per step, plus the single prime.
    assert calls["n"] == n_steps + 1
    assert abs((ef - e0) / e0) < 1e-3
    assert np.linalg.norm(pos_kdk[0] - pos0[0]) / r < 0.02  # returns near start

    # --- reference: the naive inline 2-eval loop must give the same trajectory ---
    parts2 = _two_body_state(ti, m, r, v_each)
    bx = ti.field(ti.f32, shape=2)
    by = ti.field(ti.f32, shape=2)
    bz = ti.field(ti.f32, shape=2)

    def accel2():
        direct_accel(
            parts2.pos_x, parts2.pos_y, parts2.pos_z, parts2.mass, bx, by, bz, 2, units.G, eps**2
        )

    for _ in range(n_steps):
        accel2()
        kick(parts2.vel_x, parts2.vel_y, parts2.vel_z, bx, by, bz, 2, dt / 2)
        drift(
            parts2.pos_x, parts2.pos_y, parts2.pos_z,
            parts2.vel_x, parts2.vel_y, parts2.vel_z, 2, dt,
        )
        accel2()
        kick(parts2.vel_x, parts2.vel_y, parts2.vel_z, bx, by, bz, 2, dt / 2)
    pos_inline, vel_inline = _read_state(parts2)
    # EXACT equality, not allclose: the cached-accel kdk_step is bit-identical to the naive
    # 2-eval KDK loop (the redundant top-of-loop eval recomputes the same a at unchanged x).
    # A loose tolerance would also pass a DKD reordering (it differs by only ~3e-5) — exact
    # equality is what actually locks the *ordering*, giving this an anti-regression bite.
    assert np.array_equal(pos_kdk, pos_inline)
    assert np.array_equal(vel_kdk, vel_inline)


def test_potential_energy_grid_matches_gaussian_self_energy():
    """½ΣρΦ·dV from the solved grid vs the analytic Gaussian self-energy −GM²/(2σ√π).

    Covers the *production* PE diagnostic end-to-end (deposit-free: ρ filled analytically
    to isolate the grid PE + solver), and exercises a dx≠1 grid.
    """
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver.fft_oracle import FFTPoissonSolver

    ti.init(arch=ti.cpu)
    n, dx, sigma, mass = 64, 2.0, 10.0, 2.0e11  # dx≠1 on purpose
    axis = (np.arange(n) - (n - 1) / 2.0) * dx
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    g = np.exp(-(xx**2 + yy**2 + zz**2) / (2.0 * sigma**2))
    rho_np = (g / (g.sum() * dx**3) * mass).astype(np.float32)  # density, ∑ρ·dV = mass
    rho = ti.field(ti.f32, shape=(n, n, n))
    rho.from_numpy(rho_np)
    phi = ti.field(ti.f32, shape=(n, n, n))
    FFTPoissonSolver(n, dx=dx).solve(rho, phi)

    w_grid = diagnostics.potential_energy_grid(rho.to_numpy(), phi.to_numpy(), dx=dx)
    w_analytic = -units.G * mass**2 / (2.0 * sigma * np.sqrt(np.pi))
    assert w_grid == pytest.approx(w_analytic, rel=0.01)


def test_potential_to_accel_dx_not_one():
    """The grad step is dx-correct: physical gradient, not index gradient."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.deposit import potential_to_accel

    ti.init(arch=ti.cpu)
    n, dx = 16, 2.0
    ii, jj, kk = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    # phi = 2·x_phys + 3·y_phys, with x_phys = i·dx -> g = (−2, −3, 0) everywhere.
    phi = ti.field(ti.f32, shape=(n, n, n))
    phi.from_numpy((2.0 * ii * dx + 3.0 * jj * dx).astype(np.float32))
    axg = ti.field(ti.f32, shape=(n, n, n))
    ayg = ti.field(ti.f32, shape=(n, n, n))
    azg = ti.field(ti.f32, shape=(n, n, n))
    potential_to_accel(phi, axg, ayg, azg, dx=dx)
    interior = (slice(1, n - 1),) * 3
    assert np.allclose(axg.to_numpy()[interior], -2.0, atol=1e-4)
    assert np.allclose(ayg.to_numpy()[interior], -3.0, atol=1e-4)


def test_force_near_boundary_not_dropped():
    """A particle within one cell of a face still gathers the full force.

    Regression guard: the grad step used to zero boundary nodes, which silently halved the
    force of any particle whose CIC stencil touched a face. One-sided differences fix it.
    """
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import ParticleState
    from galaxy_collision.deposit import gather_acceleration, potential_to_accel

    ti.init(arch=ti.cpu)
    n, dx = 16, 1.0
    ii, _, _ = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    phi = ti.field(ti.f32, shape=(n, n, n))
    phi.from_numpy((2.0 * ii).astype(np.float32))  # g = -2 everywhere
    axg = ti.field(ti.f32, shape=(n, n, n))
    ayg = ti.field(ti.f32, shape=(n, n, n))
    azg = ti.field(ti.f32, shape=(n, n, n))
    potential_to_accel(phi, axg, ayg, azg, dx=dx)

    # Particle at grid coord 14.5: its x-stencil spans nodes 14 and 15 (15 is a face node).
    parts = ParticleState(1)
    parts.pos_x.from_numpy(np.array([14.5], dtype=np.float32))
    parts.pos_y.from_numpy(np.array([8.0], dtype=np.float32))
    parts.pos_z.from_numpy(np.array([8.0], dtype=np.float32))
    parts.mass.from_numpy(np.array([1.0], dtype=np.float32))
    acc_x = ti.field(ti.f32, shape=1)
    acc_y = ti.field(ti.f32, shape=1)
    acc_z = ti.field(ti.f32, shape=1)
    gather_acceleration(parts, axg, ayg, azg, acc_x, acc_y, acc_z, dx=dx)
    assert acc_x.to_numpy()[0] == pytest.approx(-2.0, abs=1e-4)  # full force, not -1.0


def test_diagnostics_momenta_known_config():
    m = np.array([2.0, 3.0])
    pos = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    vel = np.array([[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]])
    # Linear momentum: 2·(0,1,0) + 3·(0,-1,0) = (0,-1,0).
    assert np.allclose(diagnostics.linear_momentum(m, vel), [0.0, -1.0, 0.0])
    # Angular momentum about origin: 2·(r×v) + 3·(r×v); each r×v = (1,0,0)×(0,1,0)=(0,0,1)
    # and (-1,0,0)×(0,-1,0)=(0,0,1) -> total (0,0,5).
    assert np.allclose(diagnostics.angular_momentum(m, pos, vel), [0.0, 0.0, 5.0])
    assert diagnostics.kinetic_energy(m, vel) == pytest.approx(0.5 * (2.0 + 3.0))


def test_potential_energy_direct_two_body():
    m = np.array([1.0e10, 2.0e10])
    pos = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]])
    expected = -units.G * m[0] * m[1] / 4.0
    assert diagnostics.potential_energy_direct(m, pos) == pytest.approx(expected, rel=1e-10)
