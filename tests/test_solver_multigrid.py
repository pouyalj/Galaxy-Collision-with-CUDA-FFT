"""Open-BC multigrid: convergence, an analytic anchor, and agreement with the oracle.

AGENT.md §6 test 2 / §11 RV9. Three layers, ordered by how *principled* (fixture-independent)
the check is:

1. **Convergence (fixture-independent).** A V-cycle reduces the residual, and the production
   cycle count drives the residual norm well below the RHS norm — i.e. multigrid actually
   solves *its* discrete equation. This is solver correctness, independent of the source.
2. **Analytic far field (derived tolerance).** A point mass must produce Φ(r) → −GM/r. The
   7-point discrete Laplacian's lattice Green's function approaches the continuum 1/r from
   below with a measured **O((dx/r)²) ≈ 0.26/r²** error (1.98% at r=4 → 0.06% at r=20),
   *decreasing monotonically* with r. We assert exactly that law (cap 0.4/r²) — a tolerance
   *derived* from the discretization, not tuned. The oracle, convolving the exact 1/r kernel,
   matches to machine precision (cross-checks test 1).
3. **Multigrid ≈ FFT oracle (loose sanity).** On the same source the two solvers differ only
   by the discretization gap of (2): for a σ=4 Gaussian it *plateaus at ~1.7% peak* once
   converged, independent of cycle count. This is fixture-sensitive (sharper source → larger
   core gap), so it is kept as a loose sanity cap; the tight checks are (1) and (2).

All need Taichi.
"""

from __future__ import annotations

import numpy as np
import pytest

from galaxy_collision import units


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


def _resid_ratio(solver, rho_np, grav, dx):
    """Finest-grid residual norm relative to the RHS norm ‖4πGρ‖ — a scale-free convergence
    measure. < 1e-3 means the discrete equation is solved to 0.1% of its own source."""
    rhs_norm = 4.0 * np.pi * grav * float(np.sqrt((rho_np.astype(np.float64) ** 2).sum()))
    return solver.residual_norm() / rhs_norm


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


def test_multigrid_point_mass_far_field_monopole():
    """Analytic anchor (RV9): multigrid's potential → −GM/r with the *derived* O((dx/r)²)
    discrete-Green's-function error, decreasing monotonically. This pins multigrid to ground
    truth directly (not just against the oracle), with a tolerance read off the discretization
    rather than tuned."""
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver.multigrid import MultigridPoissonSolver

    ti.init(arch=ti.cpu)
    n, dx, mass, grav = 96, 1.0, 1.0e10, units.G
    c = n // 2
    rho = ti.field(ti.f32, shape=(n, n, n))
    phi = ti.field(ti.f32, shape=(n, n, n))
    rho.fill(0.0)
    rho[c, c, c] = mass / dx**3  # all mass on one node — the hardest case for the lattice

    solver = MultigridPoissonSolver(n, dx=dx, grav_constant=grav, n_cycles=80)
    solver.solve(rho, phi)
    # Converged to its own equation (so the residual deviation below is purely discretization).
    assert _resid_ratio(solver, rho.to_numpy(), grav, dx) < 1e-4

    phi_np = phi.to_numpy()
    prev = np.inf
    for r in (4, 6, 8, 12, 16, 20):
        analytic = -grav * mass / (r * dx)
        # The 6 samples all lie on coordinate axes, so by cubic symmetry they are equal up to
        # fp32 read noise — averaging cancels that noise. (Axis is the worst-case direction for
        # the lattice Green's function, so the measured ~0.26/r² is the largest coefficient.)
        vals = [phi_np[c + r, c, c], phi_np[c - r, c, c], phi_np[c, c + r, c],
                phi_np[c, c - r, c], phi_np[c, c, c + r], phi_np[c, c, c - r]]
        rel = abs(float(np.mean(vals)) - analytic) / abs(analytic)
        assert rel < 0.4 / r**2, f"r={r}: deviation {rel:.4f} exceeds the 0.4/r² discretization law"
        # Monotone decreasing, with an fp32-noise-scaled relative slack.
        assert rel < prev * (1.0 + 1e-4) + 1e-7, f"r={r}: deviation {rel:.4f} not decreasing"
        prev = rel


def test_multigrid_matches_fft_oracle():
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver.fft_oracle import FFTPoissonSolver
    from galaxy_collision.solver.multigrid import MultigridPoissonSolver

    ti.init(arch=ti.cpu)
    n, dx = 48, 1.0
    rho = _gaussian_rho(ti, n, dx, sigma=4.0, total_mass=2.0e11)
    rho_np = rho.to_numpy()

    phi_mg = ti.field(ti.f32, shape=(n, n, n))
    phi_fft = ti.field(ti.f32, shape=(n, n, n))
    mg = MultigridPoissonSolver(n, dx=dx, n_cycles=40)
    mg.solve(rho, phi_mg)
    FFTPoissonSolver(n, dx=dx).solve(rho, phi_fft)

    # Principled gate: multigrid has actually converged (residual ≪ source). The residual
    # deviation from the oracle below is therefore the discretization gap, not iteration error.
    assert _resid_ratio(mg, rho_np, units.G, dx) < 1e-3

    a = phi_mg.to_numpy()
    b = phi_fft.to_numpy()
    # Loose sanity cap on the interior (exclude a boundary margin). The peak deviation is the
    # discrete-Laplacian vs continuum-1/r gap (≈1.7% here); it is fixture-dependent (sharper
    # source → larger core gap), so this cap is generous on purpose — the tight, principled
    # checks live in test_multigrid_point_mass_far_field_monopole and the convergence gate above.
    m = 6
    core = (slice(m, n - m),) * 3
    da, db = a[core], b[core]
    peak = np.abs(b).max()
    rel_max = np.abs(da - db).max() / peak
    rel_mean = np.abs(da - db).mean() / peak
    assert rel_max < 0.03, f"peak-relative deviation {rel_max:.4f} exceeds 3%"
    assert rel_mean < 0.005, f"mean-relative deviation {rel_mean:.4f} exceeds 0.5%"


def test_multigrid_matches_fft_oracle_off_center():
    """Off-center source stresses the multipole open-BC faces (dipole/quadrupole non-zero).

    The centered-blob test leaves the boundary near-spherically-symmetric (monopole nearly
    exact); an off-center mass is the harder case for the face expansion, so this guards the
    open-BC quality the production solver relies on.
    """
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.solver.fft_oracle import FFTPoissonSolver
    from galaxy_collision.solver.multigrid import MultigridPoissonSolver

    ti.init(arch=ti.cpu)
    n, dx = 64, 1.0
    axis = np.arange(n) * dx
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    cx, cy, cz = 40.0, 28.0, 34.0  # deliberately off the box center
    g = np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2 + (zz - cz) ** 2) / (2.0 * 4.0**2))
    rho_np = (g / (g.sum() * dx**3) * 2.0e11).astype(np.float32)
    rho = ti.field(ti.f32, shape=(n, n, n))
    rho.from_numpy(rho_np)

    phi_mg = ti.field(ti.f32, shape=(n, n, n))
    phi_fft = ti.field(ti.f32, shape=(n, n, n))
    MultigridPoissonSolver(n, dx=dx, n_cycles=50).solve(rho, phi_mg)
    FFTPoissonSolver(n, dx=dx).solve(rho, phi_fft)

    a, b = phi_mg.to_numpy(), phi_fft.to_numpy()
    m = 6
    core = (slice(m, n - m),) * 3
    rel_max = np.abs(a[core] - b[core]).max() / np.abs(b).max()
    assert rel_max < 0.03, f"off-center peak-relative deviation {rel_max:.4f} exceeds 3%"


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
