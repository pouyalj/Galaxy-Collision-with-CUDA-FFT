"""FFT oracle vs the analytic point-mass potential (AGENT.md §6 test 1).

A single mass M deposited on one node must produce Φ(r) = −G·M/r at every other node
(the free-space Green's function the oracle convolves with is exactly −G·dx³/r, so this
also confirms the zero-padding kills circular wraparound and the indexing is right).
Needs Taichi only to hold the grid fields; the solve itself is NumPy.
"""

from __future__ import annotations

import pytest

from galaxy_collision import units


def test_fft_oracle_matches_point_mass_potential():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver.fft_oracle import FFTPoissonSolver

    ti.init(arch=ti.cpu)
    n, dx, mass = 48, 1.0, 1.0e10
    c = n // 2
    rho = ti.field(ti.f32, shape=(n, n, n))
    phi = ti.field(ti.f32, shape=(n, n, n))
    rho.fill(0.0)
    rho[c, c, c] = mass / dx**3  # one node carries all the mass

    FFTPoissonSolver(n, dx=dx).solve(rho, phi)
    phi_np = phi.to_numpy()

    # Sample along +x at a range of distances (well away from the box faces).
    for r in (2, 4, 8, 12):
        analytic = -units.G * mass / (r * dx)
        assert phi_np[c + r, c, c] == pytest.approx(analytic, rel=2e-3)

    # Isotropy: the three axes agree.
    assert phi_np[c + 6, c, c] == pytest.approx(phi_np[c, c + 6, c], rel=1e-4)
    assert phi_np[c + 6, c, c] == pytest.approx(phi_np[c, c, c + 6], rel=1e-4)


def test_fft_oracle_rejects_pad_factor_one():
    from galaxy_collision.solver.fft_oracle import FFTPoissonSolver

    with pytest.raises(ValueError, match="pad_factor"):
        FFTPoissonSolver(16, pad_factor=1)
