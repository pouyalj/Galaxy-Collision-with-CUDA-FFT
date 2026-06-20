"""Structure-of-Arrays (SoA) particle & grid state, plus a memory estimator.

SoA is mandatory (AGENT.md D3): every field is a flat, contiguous Taichi field,
never an array-of-pointers (the legacy ``float**`` layout — bug #8 — is exactly what
made the headline 1e8 config unrunnable). Positions/velocities live in internal
units (kpc, kpc/Myr); mass in M_sun; ``gid`` tags the galaxy (0 = MW, 1 = Andromeda).

The :func:`estimate_memory` helper is pure Python (no Taichi needed) so callers can
budget VRAM for a chosen N or particle mass *before* allocating anything.
"""

from __future__ import annotations

from dataclasses import dataclass

# Bytes per dtype in the SoA layout.
_F32 = 4
_I32 = 4

# Particle fields: pos(x,y,z) + vel(x,y,z) + mass, all f32; gid as i32.
BYTES_PER_PARTICLE = 6 * _F32 + _F32 + _I32  # 32 B (28 B without mass; see AGENT.md §5.2)

# Grid fields for the production (multigrid) path, counted in units of G^3 f32 arrays
# (AGENT.md §5.2, §11 RV4b — resolved at Stage 3):
#   - rho + phi                          → 2
#   - acceleration field g = −∇Φ (x,y,z) → 3
#   - multigrid hierarchy (phi/rhs/res across levels) ≈ 3 × Σ_l 8^{−l} = 3 × 8/7 ≈ 3.43
# The zero-padded FFT *oracle* adds a (2G)^3 complex pad (~8× a real G^3 grid) but is
# NVIDIA-validation-only, so it is intentionally excluded from the production estimate.
_GRID_RHO_PHI = 2
_GRID_ACCEL = 3
_GRID_MG_HIERARCHY = 3.0 * 8.0 / 7.0
_GRID_FIELDS = _GRID_RHO_PHI + _GRID_ACCEL + _GRID_MG_HIERARCHY  # ≈ 8.43 G^3 f32 arrays


@dataclass(frozen=True)
class MemoryEstimate:
    """Predicted memory footprint for a given configuration, in bytes."""

    n_particles: int
    grid_size: int
    particle_bytes: int
    grid_bytes: int

    @property
    def total_bytes(self) -> int:
        return self.particle_bytes + self.grid_bytes

    @property
    def total_gb(self) -> float:
        return self.total_bytes / 1e9

    def summary(self) -> str:
        return (
            f"N={self.n_particles:,} particles, grid={self.grid_size}^3 | "
            f"particles {self.particle_bytes / 1e9:.2f} GB + "
            f"grid {self.grid_bytes / 1e9:.2f} GB = {self.total_gb:.2f} GB"
        )


def estimate_memory(n_particles: int, grid_size: int = 256) -> MemoryEstimate:
    """Estimate the SoA memory footprint for the production (multigrid) path.

    Counts particles plus the grid buffers actually allocated by a Stage-3 run on the
    multigrid solver: rho, phi, the acceleration field, and the multigrid hierarchy
    (see ``_GRID_FIELDS``). It does **not** include the zero-padded FFT *oracle* buffer,
    which is NVIDIA-validation-only and allocated only when that solver is selected.
    """
    # Positive-only, consistent with SimConfig and ParticleState (AGENT.md §11 RV4).
    if n_particles <= 0:
        raise ValueError(f"n_particles must be positive, got {n_particles}")
    if grid_size <= 0:
        raise ValueError(f"grid_size must be positive, got {grid_size}")
    particle_bytes = n_particles * BYTES_PER_PARTICLE
    grid_bytes = int(_GRID_FIELDS * grid_size**3 * _F32)
    return MemoryEstimate(
        n_particles=n_particles,
        grid_size=grid_size,
        particle_bytes=particle_bytes,
        grid_bytes=grid_bytes,
    )


class ParticleState:
    """SoA particle fields as Taichi scalar fields. Requires ``ti.init`` first."""

    def __init__(self, n: int):
        import taichi as ti

        if n <= 0:
            raise ValueError(f"particle count must be positive, got {n}")
        self.n = n
        self.pos_x = ti.field(dtype=ti.f32, shape=n)
        self.pos_y = ti.field(dtype=ti.f32, shape=n)
        self.pos_z = ti.field(dtype=ti.f32, shape=n)
        self.vel_x = ti.field(dtype=ti.f32, shape=n)
        self.vel_y = ti.field(dtype=ti.f32, shape=n)
        self.vel_z = ti.field(dtype=ti.f32, shape=n)
        self.mass = ti.field(dtype=ti.f32, shape=n)
        self.gid = ti.field(dtype=ti.i32, shape=n)

    def fields(self):
        """Return all fields in a stable order (useful for bulk ops/tests)."""
        return (
            self.pos_x,
            self.pos_y,
            self.pos_z,
            self.vel_x,
            self.vel_y,
            self.vel_z,
            self.mass,
            self.gid,
        )


class GridState:
    """SoA grid fields (density rho, potential phi). Requires ``ti.init`` first."""

    def __init__(self, grid_size: int = 256):
        import taichi as ti

        if grid_size <= 0:
            raise ValueError(f"grid_size must be positive, got {grid_size}")
        self.grid_size = grid_size
        self.rho = ti.field(dtype=ti.f32, shape=(grid_size, grid_size, grid_size))
        self.phi = ti.field(dtype=ti.f32, shape=(grid_size, grid_size, grid_size))


__all__ = [
    "BYTES_PER_PARTICLE",
    "MemoryEstimate",
    "estimate_memory",
    "ParticleState",
    "GridState",
]
