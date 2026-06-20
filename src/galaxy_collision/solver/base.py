"""The pluggable Poisson-solver interface (AGENT.md D5).

A solver maps a mass density ``rho`` to the gravitational potential ``phi`` by solving

    ∇²Φ = 4πG ρ

with **open (isolated) boundary conditions** — Φ → 0 far from the mass, *not* periodic
(rejecting the legacy plain-FFT periodic solve, bug #6). Two implementations sit behind
this interface:

- :class:`~galaxy_collision.solver.multigrid.MultigridPoissonSolver` — the portable
  default (real-space V-cycles, open BCs via a multipole boundary expansion).
- :class:`~galaxy_collision.solver.fft_oracle.FFTPoissonSolver` — the zero-padded
  isolated-Green's-function FFT oracle used to validate multigrid and reproduce the
  paper's spectral results.

Both operate on **node-centered** grids (see ``deposit.py``): node ``(i,j,k)`` is at
physical ``(i,j,k)·dx``. ``rho`` and ``phi`` are square Taichi fields of side
``grid_size``; ``solve`` reads ``rho`` and overwrites ``phi`` in place.
"""

from __future__ import annotations

from abc import ABC, abstractmethod


class PoissonSolver(ABC):
    """Maps mass density → gravitational potential (∇²Φ = 4πGρ, open BCs)."""

    grid_size: int
    dx: float

    @abstractmethod
    def solve(self, rho, phi, warm_start: bool = False) -> None:
        """Read density field ``rho``, write potential field ``phi`` (both node-centered).

        When ``warm_start`` is True, iterative solvers may use the current contents of
        ``phi`` as the initial guess (cheap for a time loop where Φ changes little per
        step). Direct solvers ignore it.
        """
        raise NotImplementedError


__all__ = ["PoissonSolver"]
