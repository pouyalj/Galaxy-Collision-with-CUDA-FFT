"""Open-BC multigrid: convergence + agreement with the FFT oracle (AGENT.md §6 test 2).

The headline check is test 2: the *same* density through multigrid and through the
zero-padded FFT oracle must agree within tolerance. We use a well-resolved Gaussian
blob (width ≫ cell size) centered in the box: there the 7-point discrete Laplacian and
the oracle's exact 1/r kernel both approach the continuum, and the centered mass makes
the monopole/quadrupole boundary expansion accurate — so the two solvers should agree
to ~1%. A V-cycle must also actually reduce the residual.

All need Taichi.
"""

from __future__ import annotations

import numpy as np
import pytest


def _gaussian_rho(ti, n, dx, sigma, total_mass):
    """A normalized Gaussian mass blob centered in the box, as a Taichi rho field."""
    axis = (np.arange(n) - (n - 1) / 2.0) * dx
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    g = np.exp(-(xx**2 + yy**2 + zz**2) / (2.0 * sigma**2))
    mass = g / g.sum() * total_mass  # mass per node
    rho_np = (mass / dx**3).astype(np.float32)  # density
    rho = ti.field(ti.f32, shape=(n, n, n))
    rho.from_numpy(rho_np)
    return rho


def test_vcycle_reduces_residual():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver.multigrid import MultigridPoissonSolver

    ti.init(arch=ti.cpu)
    n, dx = 32, 1.0
    rho = _gaussian_rho(ti, n, dx, sigma=3.0, total_mass=1.0e10)
    phi = ti.field(ti.f32, shape=(n, n, n))

    solver = MultigridPoissonSolver(n, dx=dx, n_cycles=1)
    solver.solve(rho, phi)
    r1 = solver.residual_norm()
    for _ in range(6):
        solver._vcycle(0)
    r2 = solver.residual_norm()
    # Genuine multigrid acceleration: 6 V-cycles cut the residual by >10x. (Pure
    # Gauss-Seidel smoothing alone would barely dent the low-frequency error — a
    # broken coarse-grid correction shows up here as a stall.)
    assert r2 < 0.1 * r1


def test_multigrid_matches_fft_oracle():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver.fft_oracle import FFTPoissonSolver
    from galaxy_collision.solver.multigrid import MultigridPoissonSolver

    ti.init(arch=ti.cpu)
    n, dx = 48, 1.0
    rho = _gaussian_rho(ti, n, dx, sigma=4.0, total_mass=2.0e11)

    phi_mg = ti.field(ti.f32, shape=(n, n, n))
    phi_fft = ti.field(ti.f32, shape=(n, n, n))
    MultigridPoissonSolver(n, dx=dx, n_cycles=40).solve(rho, phi_mg)
    FFTPoissonSolver(n, dx=dx).solve(rho, phi_fft)

    a = phi_mg.to_numpy()
    b = phi_fft.to_numpy()
    # Compare on the interior (exclude a boundary margin where the BC approximation and
    # the oracle's open kernel differ most). Residual difference is dominated by the
    # genuine discrete-Laplacian vs continuum-1/r-kernel gap, not multigrid error.
    m = 6
    core = (slice(m, n - m),) * 3
    da, db = a[core], b[core]
    peak = np.abs(b).max()
    rel_max = np.abs(da - db).max() / peak
    rel_mean = np.abs(da - db).mean() / peak
    assert rel_max < 0.03, f"peak-relative deviation {rel_max:.4f} exceeds 3%"
    assert rel_mean < 0.005, f"mean-relative deviation {rel_mean:.4f} exceeds 0.5%"


def test_solver_factory():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver import make_solver
    from galaxy_collision.solver.fft_oracle import FFTPoissonSolver
    from galaxy_collision.solver.multigrid import MultigridPoissonSolver

    ti.init(arch=ti.cpu)
    assert isinstance(make_solver("multigrid", 16), MultigridPoissonSolver)
    assert isinstance(make_solver("fft", 16), FFTPoissonSolver)
    with pytest.raises(ValueError, match="unknown solver"):
        make_solver("nope", 16)
