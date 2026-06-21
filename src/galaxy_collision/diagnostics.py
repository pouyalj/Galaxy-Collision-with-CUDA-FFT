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
conserves energy/momentum to PM tolerance is a separate, looser check — the **Stage-3 exit
gate** (Checkpoint E: a stable Plummer sphere with bounded energy drift), not something the
direct-path Kepler/Plummer-IC tests cover.

**On the PM energy metric (important):** ``potential_energy_grid`` = ½ΣρΦ·dV is the
standard PM energy, but it is *not* the exact discrete Hamiltonian whose gradient (CIC
gather of −∇Φ) generates the integrated force. So PM ``energy_drift`` has an irreducible
grid-discretization floor that does *not* shrink with ``dt`` — it is a PM-tolerance check,
not a proof of symplectic conservation. (The trajectory *is* symplectic; the Kepler test
proves that on the direct path, where energy and force share one Hamiltonian.)

**Big-run energy diagnostics (RV11).** For the two-galaxy runs (10–100M particles) the O(N²)
direct PE is infeasible, so the grid PE is the only available potential-energy metric. How to use
it: the grid PE ½ΣρΦ·dV recovers the physical (direct, ε≈1 cell) PE up to a *fixed* discretization
offset of the estimator — ~2% for the FFT oracle, ~4% for multigrid, stable across realizations
(this offset is locked by ``tests/test_diagnostics.py``). That offset, not a self-energy term, is
the bulk of the grid↔physical gap: the FFT oracle zeros only the *on-node* point self-term
(``g(r=0)=0``), while the CIC *cloud* self-energy — a particle's mass spread over up to 8 nodes —
survives for *both* solvers and is sub-cell-position dependent (so it fluctuates as particles move
across a cell); probes put it well under a percent of the grid PE. Subtracting it would shift the
absolute number by ~0.1% and would *not* remove the run-to-run scatter, so we **don't** subtract it.

The big-run conservation metric is therefore the **drift** of the grid energy used as a *monitoring*
signal — never a pass/fail gate. RV7 retired grid-PE drift as the Stage-3 exit gate precisely
because its realization scatter crosses ~1% on some seeds (the Plummer exit gate uses the direct PE
instead); that scatter, not the small self-energy, is the caveat to keep in mind here. The
offset-insensitive **virial** ratio T/|W| (recorded each sample) is the robust supporting check:
≈½ for an isolated virialized object, and for a *collision* a smooth quantity to track for
consistency — not a fixed-½ gate, since the global T includes the bulk approach KE (subtract the
bulk frame for a single-remnant virial).

Per the precision policy (D6) diagnostics are computed in **fp64** — here on the host in
NumPy, fed by ``field.to_numpy()`` at the output cadence. They are not in the per-step hot
path, so host fp64 is both correct and clear (a Metal-safe fp32/Kahan path is a later
concern; Metal has no hardware fp64).

Two flavors of potential energy:
- :func:`potential_energy_direct` — exact softened pair sum −½G Σ_{i≠j} mᵢmⱼ/√(r²+ε²),
  O(N²). Used to validate the integrator on small systems (e.g. two-body Kepler).
- :func:`potential_energy_grid` — the PM estimate ½ Σ ρΦ·dV from the solved grid, O(grid).
  This is what a production run tracks (a PM-tolerance metric; see the note above).

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


def center_of_mass(mass: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Mass-weighted mean position, as a length-3 fp64 vector."""
    mass = np.asarray(mass, dtype=np.float64)
    pos = np.asarray(pos, dtype=np.float64)
    return (mass[:, None] * pos).sum(axis=0) / mass.sum()


def lagrangian_radius(
    mass: np.ndarray, pos: np.ndarray, fraction: float = 0.5, center: np.ndarray | None = None
) -> float:
    """Radius enclosing ``fraction`` of the total mass about ``center`` (default: the CM).

    The half-mass radius (fraction=0.5) is the headline structural diagnostic for the
    Plummer-stability test: a sphere that collapses or evaporates moves it.
    """
    mass = np.asarray(mass, dtype=np.float64)
    pos = np.asarray(pos, dtype=np.float64)
    c = center_of_mass(mass, pos) if center is None else np.asarray(center, dtype=np.float64)
    r = np.linalg.norm(pos - c, axis=1)
    order = np.argsort(r)
    cum = np.cumsum(mass[order])
    target = fraction * mass.sum()
    idx = int(np.searchsorted(cum, target))
    return float(r[order][min(idx, len(r) - 1)])


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
    """PM potential energy estimate ½ Σ ρΦ·dV from the solved grid.

    A PM-tolerance metric, not the exact integrated Hamiltonian — see the module docstring:
    energy drift built on it has a dt-independent grid-discretization floor.
    """
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


def total_energy_grid(
    mass: np.ndarray, vel: np.ndarray, rho: np.ndarray, phi: np.ndarray, dx: float = 1.0
) -> float:
    """KE + PM grid PE — the big-run total energy (RV11). Absolute value carries the documented
    ~2–4% grid-PE discretization offset; track its *drift*, not the absolute number."""
    return kinetic_energy(mass, vel) + potential_energy_grid(rho, phi, dx)


def virial_ratio(kinetic: float, potential: float) -> float:
    """T/|W| — ≈ 0.5 for an isolated system in virial equilibrium. A robust big-run check (RV11):
    dimensionless and insensitive to the grid-PE offset. For a *collision* the global T includes
    the bulk approach KE, so it is not ≈0.5 there — track its evolution as a consistency monitor
    (subtract the bulk frame for a meaningful single-remnant virial)."""
    return kinetic / abs(potential) if potential != 0.0 else float("nan")


__all__ = [
    "kinetic_energy",
    "linear_momentum",
    "angular_momentum",
    "center_of_mass",
    "lagrangian_radius",
    "potential_energy_direct",
    "potential_energy_grid",
    "total_energy_direct",
    "total_energy_grid",
    "virial_ratio",
]
