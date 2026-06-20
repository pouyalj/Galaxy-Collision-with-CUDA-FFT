"""Conservation diagnostics: energy, momentum, angular momentum, tracer (Stage 3).

AGENT.md §5.5/§6: energy drift is the headline correctness metric, with linear and
angular momentum as supporting checks. **Conservation expectations differ by force path:**

- **Linear momentum** is *approximately* conserved on the PM path. Zero net self-force
  needs *two* symmetries together: the symmetric CIC deposit/gather pair **and** the
  antisymmetric central-difference gradient (interior). Drop either — e.g. a one-sided
  gradient — and a spurious self-force appears. It is *exactly* conserved on the direct
  pairwise path (Newton's third law: equal-and-opposite forces).
- **Angular momentum** is conserved *exactly only on the direct path*, where the force is
  central. On the PM **grid** it is conserved only *approximately*: grid anisotropy makes
  the interpolated inter-particle force slightly non-central, so a small torque leaks in.

So the strict ΔL≈0 check the two-body Kepler test passes uses the direct central-force
path. Validating that the production **PM force chain** (deposit→solve→grad→gather)
conserves energy/momentum to PM tolerance is a separate, looser check — the Stage-E exit
gate (a stable Plummer sphere with bounded energy drift), not something Checkpoint C's
direct-path tests cover.

Per the precision policy (D6) diagnostics are computed in **fp64** — here on the host in
NumPy, fed by ``field.to_numpy()`` at the output cadence. They are not in the per-step hot
path, so host fp64 is both correct and clear (a Metal-safe fp32/Kahan path is a later
concern; Metal has no hardware fp64).

Two flavors of potential energy:
- :func:`potential_energy_direct` — exact softened pair sum −½G Σ_{i≠j} mᵢmⱼ/√(r²+ε²),
  O(N²). Used to validate the integrator on small systems (e.g. two-body Kepler).
- :func:`potential_energy_grid` — the PM estimate ½ Σ ρΦ·dV from the solved grid, O(grid).
  This is what a production run tracks.

All position/velocity arrays are ``(N, 3)``; mass is ``(N,)``. Inputs are cast to fp64.
"""

from __future__ import annotations

import numpy as np

from galaxy_collision import units


def kinetic_energy(mass: np.ndarray, vel: np.ndarray) -> float:
    """Total kinetic energy ½ Σ mᵢ |vᵢ|² (M_sun·kpc²/Myr²)."""
    mass = np.asarray(mass, dtype=np.float64)
    vel = np.asarray(vel, dtype=np.float64)
    return float(0.5 * np.sum(mass * np.sum(vel**2, axis=1)))


def linear_momentum(mass: np.ndarray, vel: np.ndarray) -> np.ndarray:
    """Total linear momentum Σ mᵢ vᵢ, as a length-3 fp64 vector."""
    mass = np.asarray(mass, dtype=np.float64)
    vel = np.asarray(vel, dtype=np.float64)
    return (mass[:, None] * vel).sum(axis=0)


def angular_momentum(
    mass: np.ndarray, pos: np.ndarray, vel: np.ndarray, center: np.ndarray | None = None
) -> np.ndarray:
    """Total angular momentum Σ mᵢ (rᵢ−c) × vᵢ about ``center`` (default: origin)."""
    mass = np.asarray(mass, dtype=np.float64)
    pos = np.asarray(pos, dtype=np.float64)
    vel = np.asarray(vel, dtype=np.float64)
    r = pos if center is None else pos - np.asarray(center, dtype=np.float64)
    return (mass[:, None] * np.cross(r, vel)).sum(axis=0)


def potential_energy_direct(
    mass: np.ndarray, pos: np.ndarray, softening: float = 0.0, grav: float | None = None
) -> float:
    """Exact Plummer-softened potential energy −½G Σ_{i≠j} mᵢmⱼ/√(rᵢⱼ²+ε²) (O(N²))."""
    g = units.G if grav is None else grav
    mass = np.asarray(mass, dtype=np.float64)
    pos = np.asarray(pos, dtype=np.float64)
    diff = pos[:, None, :] - pos[None, :, :]
    r2 = (diff**2).sum(axis=-1) + softening**2
    np.fill_diagonal(r2, np.inf)  # drop self-interaction (1/√inf = 0, no warning)
    inv_r = 1.0 / np.sqrt(r2)
    mm = mass[:, None] * mass[None, :]
    return float(-0.5 * g * np.sum(mm * inv_r))


def potential_energy_grid(rho: np.ndarray, phi: np.ndarray, dx: float = 1.0) -> float:
    """PM potential energy estimate ½ Σ ρΦ·dV from the solved grid."""
    rho = np.asarray(rho, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    # dx**3 is the cell volume of the ½∫ρΦ integral. It is unrelated to (not a
    # double-count of) the dx**3 the FFT solver folds into its Green's kernel: that one
    # turns ρ into per-cell mass inside the convolution; Φ here is already the potential.
    return float(0.5 * np.sum(rho * phi) * dx**3)


def total_energy_direct(
    mass: np.ndarray,
    pos: np.ndarray,
    vel: np.ndarray,
    softening: float = 0.0,
    grav: float | None = None,
) -> float:
    """KE + direct softened PE (the conserved quantity for the two-body Kepler test)."""
    return kinetic_energy(mass, vel) + potential_energy_direct(mass, pos, softening, grav)


__all__ = [
    "kinetic_energy",
    "linear_momentum",
    "angular_momentum",
    "potential_energy_direct",
    "potential_energy_grid",
    "total_energy_direct",
]
