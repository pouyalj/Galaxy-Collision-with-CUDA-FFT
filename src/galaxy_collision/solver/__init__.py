"""Poisson solvers (pluggable, per AGENT.md D5).

``base.PoissonSolver`` is the interface; ``multigrid.MultigridPoissonSolver`` is the
open-BC portable default and ``fft_oracle.FFTPoissonSolver`` is the zero-padded
isolated-FFT validation oracle. Implementations are imported lazily by name to avoid
importing Taichi/NumPy at package import time.
"""

from galaxy_collision.solver.base import PoissonSolver


def make_solver(name: str, grid_size: int, **kwargs) -> PoissonSolver:
    """Construct a solver by config name ('multigrid' or 'fft')."""
    if name == "multigrid":
        from galaxy_collision.solver.multigrid import MultigridPoissonSolver

        return MultigridPoissonSolver(grid_size, **kwargs)
    if name == "fft":
        from galaxy_collision.solver.fft_oracle import FFTPoissonSolver

        return FFTPoissonSolver(grid_size, **kwargs)
    raise ValueError(f"unknown solver {name!r} (expected 'multigrid' or 'fft')")


__all__ = ["PoissonSolver", "make_solver"]
