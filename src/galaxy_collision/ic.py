"""Initial-condition generation for the two-galaxy collision (Stage 2).

Builds physically-motivated ICs in internal units (kpc, Myr, M_sun), reproducible
from a seed. Per AGENT.md §5.5 and the decisions scoped 2026-06-16:

- **Particle count derived from mass (D8):** N_per_galaxy = M_galaxy / m_particle.
  Star particles share a uniform mass; m_particle is the resolution knob.
- **Disk:** sampled from the paper's Gaia-derived surface-density *shape* (Eq. 1,
  Σ(r) ∝ ℓc / √((r−rb)² + ℓc²)), normalized to the disk's mass — so total mass is
  exact by construction. Razor-thin profile with an exponential vertical scale height.
- **Bulge (~10% of mass):** a Hernquist sphere with isotropic velocity dispersion.
- **Central black hole (D9):** one massive, softened particle per galaxy carrying
  real mass (fixing the 2020 force-only kludge).
- **Disk velocities:** cold — circular orbits at the local v_c from the enclosed-mass
  rotation curve (Plummer-softened), plus the bulk approach velocity.
- **Two galaxies:** placed symmetrically about the box center, given closing
  velocities for the paper's 4v / 2v runs (v ≈ 111.7 km/s).

Modeling approximations (to be validated/tuned at Stage 3 against the equilibrium
tests, not exit-critical here): the disk's contribution to the rotation curve uses
its cylindrically-enclosed mass treated spherically (not the exact thin-disk Bessel
solution); the bulge dispersion is a virial estimate (σ_1d ≈ v_c/√3). The headline
Stage-2 guarantees — exact target mass, correct spatial profiles, reproducibility —
are independent of these.

A dark-matter halo is deferred (D10); ``halo=None`` is the documented hook.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from galaxy_collision import units
from galaxy_collision.config import GALAXY_MASS_MSUN, SimConfig

# Paper's measured approach speed v (≈125,000 mph ≈ 402,000 km/h ≈ 111.7 km/s),
# in internal units (kpc/Myr). The 4v / 2v runs scale this.
V_APPROACH_KMS = 402_000.0 / 3600.0
V_APPROACH = units.kms_to_kpc_per_myr(V_APPROACH_KMS)

_PRESET_APPROACH_FACTOR = {"two_galaxy_4v": 4.0, "two_galaxy_2v": 2.0}


@dataclass(frozen=True)
class GalaxyModel:
    """Physical parameters for one galaxy (internal units: kpc, Myr, M_sun)."""

    total_mass: float = GALAXY_MASS_MSUN  # stellar mass (disk + bulge), M_sun
    disk_fraction: float = 0.9  # remainder goes to the bulge

    # Disk surface-density shape (paper Eq. 1) + geometry.
    sigma_rb: float = 2.0  # kpc
    sigma_lc: float = 2.5  # kpc
    disk_rmax: float = 25.0  # kpc (paper: mass quoted within 25 kpc)
    disk_z0: float = 0.3  # kpc, exponential vertical scale height

    # Hernquist bulge scale radius.
    bulge_a: float = 1.0  # kpc

    # Central black hole (D9) and force softening.
    bh_mass: float = 1.0e8  # M_sun
    softening: float = 0.3  # kpc (~fraction of a 1 kpc cell)

    # Placement & kinematics (filled in when assembling the pair).
    center: tuple[float, float, float] = (128.0, 128.0, 128.0)
    bulk_velocity: tuple[float, float, float] = (0.0, 0.0, 0.0)

    # Deferred hook (D10): a static NFW halo would attach here.
    halo: None = None


@dataclass
class ICResult:
    """Generated initial conditions as SoA numpy arrays (internal units)."""

    pos: np.ndarray  # (N, 3) kpc
    vel: np.ndarray  # (N, 3) kpc/Myr
    mass: np.ndarray  # (N,) M_sun
    gid: np.ndarray  # (N,) int32 — 0 = MW, 1 = Andromeda
    preset: str
    particle_mass: float

    @property
    def n(self) -> int:
        return self.pos.shape[0]

    def summary(self) -> str:
        lines = [
            f"IC '{self.preset}': {self.n:,} particles, "
            f"m_particle={self.particle_mass:.3g} M_sun"
        ]
        for g in (0, 1):
            sel = self.gid == g
            stars = sel & (self.mass < self.particle_mass * 10)
            lines.append(
                f"  galaxy {g}: {sel.sum():,} particles, "
                f"M_star={self.mass[stars].sum():.3g} M_sun, "
                f"BH×{(sel & ~stars).sum()}"
            )
        return "\n".join(lines)


# --- Disk profile tables --------------------------------------------------------


def _disk_tables(model: GalaxyModel, n_grid: int = 4096):
    """Tabulate the disk radial CDF (fraction of disk mass within r)."""
    r = np.linspace(0.0, model.disk_rmax, n_grid)
    sigma = model.sigma_lc / np.sqrt((r - model.sigma_rb) ** 2 + model.sigma_lc**2)
    pdf = sigma * r  # mass per unit radius ∝ Σ(r)·2πr
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(r))])
    cdf = cum / cum[-1]
    return r, cdf


def _circular_velocity(
    radius: np.ndarray,
    model: GalaxyModel,
    disk_mass: float,
    bulge_mass: float,
    r_tab: np.ndarray,
    cdf_tab: np.ndarray,
) -> np.ndarray:
    """Plummer-softened circular speed from the enclosed mass at cylindrical radius."""
    m_disk = disk_mass * np.interp(radius, r_tab, cdf_tab)
    m_bulge = bulge_mass * radius**2 / (radius + model.bulge_a) ** 2
    m_enc = m_disk + m_bulge + model.bh_mass
    e2 = model.softening**2
    vc2 = units.G * m_enc * radius**2 / (radius**2 + e2) ** 1.5
    return np.sqrt(vc2)


# --- Component sampling ---------------------------------------------------------


def _sample_disk(rng, n: int, model: GalaxyModel, disk_mass: float, bulge_mass: float):
    r_tab, cdf_tab = _disk_tables(model)
    u = rng.random(n)
    radius = np.interp(u, cdf_tab, r_tab)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    z = rng.laplace(0.0, model.disk_z0, n)

    x = radius * np.cos(phi)
    y = radius * np.sin(phi)
    pos = np.stack([x, y, z], axis=1)

    # Cold circular orbits in the disk (xy) plane: v = v_c * tangential_hat.
    vc = _circular_velocity(radius, model, disk_mass, bulge_mass, r_tab, cdf_tab)
    safe = np.maximum(radius, 1e-6)
    tang = np.stack([-y / safe, x / safe, np.zeros_like(radius)], axis=1)
    vel = vc[:, None] * tang
    return pos, vel


def _sample_bulge(rng, n: int, model: GalaxyModel, disk_mass: float, bulge_mass: float):
    # Hernquist radial sampling: M(<r)/M = r²/(r+a)² = u  →  r = a√u/(1−√u).
    u = rng.random(n)
    sqrt_u = np.sqrt(u)
    radius = model.bulge_a * sqrt_u / (1.0 - sqrt_u)
    radius = np.minimum(radius, 50.0 * model.bulge_a)  # clip the rare extreme tail

    direction = rng.normal(size=(n, 3))
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    pos = radius[:, None] * direction

    # Isotropic virial-estimate dispersion (approximate; tuned at Stage 3).
    r_tab, cdf_tab = _disk_tables(model)
    vc = _circular_velocity(radius, model, disk_mass, bulge_mass, r_tab, cdf_tab)
    sigma_1d = vc / np.sqrt(3.0)
    vel = rng.normal(size=(n, 3)) * sigma_1d[:, None]
    return pos, vel


def _build_galaxy(rng, gid: int, model: GalaxyModel, particle_mass: float, n_total: int):
    """Build one galaxy's star particles + central BH (centered at origin first)."""
    n_disk = int(round(model.disk_fraction * n_total))
    n_bulge = n_total - n_disk
    if n_disk < 1 or n_bulge < 1:
        raise ValueError(
            f"n_total={n_total} too small to populate both disk and bulge "
            f"(got n_disk={n_disk}, n_bulge={n_bulge})"
        )
    disk_mass = n_disk * particle_mass
    bulge_mass = n_bulge * particle_mass

    d_pos, d_vel = _sample_disk(rng, n_disk, model, disk_mass, bulge_mass)
    b_pos, b_vel = _sample_bulge(rng, n_bulge, model, disk_mass, bulge_mass)

    pos = np.concatenate([d_pos, b_pos], axis=0)
    vel = np.concatenate([d_vel, b_vel], axis=0)
    mass = np.full(n_total, particle_mass, dtype=np.float64)

    # Central black hole: one massive, softened particle at the center, at rest
    # relative to the galaxy (bulk velocity added below).
    pos = np.concatenate([pos, np.zeros((1, 3))], axis=0)
    vel = np.concatenate([vel, np.zeros((1, 3))], axis=0)
    mass = np.concatenate([mass, [model.bh_mass]])

    # Translate to the galaxy center and add the bulk approach velocity.
    pos += np.asarray(model.center)
    vel += np.asarray(model.bulk_velocity)
    gid_arr = np.full(pos.shape[0], gid, dtype=np.int32)
    return pos, vel, mass, gid_arr


# --- Top-level assembly ---------------------------------------------------------


def build_ic(
    config: SimConfig,
    galaxy: GalaxyModel | None = None,
    separation: float = 90.0,
    impact_parameter: float = 0.0,
) -> ICResult:
    """Build two-galaxy initial conditions from a :class:`SimConfig`.

    ``separation`` (kpc) is the initial center-to-center distance along x;
    ``impact_parameter`` (kpc) offsets the galaxies perpendicular (along y) for a
    grazing encounter. Galaxies are placed symmetrically about the box center and
    given closing velocities set by the preset (4v / 2v).
    """
    if config.ic_preset not in _PRESET_APPROACH_FACTOR:
        raise NotImplementedError(
            f"IC preset {config.ic_preset!r} is not implemented in Stage 2 "
            f"(available: {sorted(_PRESET_APPROACH_FACTOR)})"
        )
    base = galaxy or GalaxyModel()

    # Resolution: derive per-galaxy N and the uniform particle mass (D8).
    if config.particle_mass is not None:
        particle_mass = config.particle_mass
        n_per_galaxy = int(round(base.total_mass / particle_mass))
    elif config.n_particles is not None:
        n_per_galaxy = config.n_particles // 2
        particle_mass = base.total_mass / n_per_galaxy
    else:
        n_per_galaxy = 50_000
        particle_mass = base.total_mass / n_per_galaxy

    # Geometry: symmetric placement about the box center.
    box_center = config.grid_size / 2.0
    half = separation / 2.0
    off = impact_parameter / 2.0
    center0 = (box_center - half, box_center - off, box_center)
    center1 = (box_center + half, box_center + off, box_center)

    # Kinematics: closing speed = factor × v, split between the two galaxies.
    factor = _PRESET_APPROACH_FACTOR[config.ic_preset]
    each = factor * V_APPROACH / 2.0
    vel0 = (each, 0.0, 0.0)  # left galaxy moves +x
    vel1 = (-each, 0.0, 0.0)  # right galaxy moves -x

    model0 = replace(base, center=center0, bulk_velocity=vel0)
    model1 = replace(base, center=center1, bulk_velocity=vel1)

    rng = np.random.default_rng(config.seed)
    p0, v0, m0, g0 = _build_galaxy(rng, 0, model0, particle_mass, n_per_galaxy)
    p1, v1, m1, g1 = _build_galaxy(rng, 1, model1, particle_mass, n_per_galaxy)

    return ICResult(
        pos=np.concatenate([p0, p1], axis=0).astype(np.float32),
        vel=np.concatenate([v0, v1], axis=0).astype(np.float32),
        mass=np.concatenate([m0, m1]).astype(np.float32),
        gid=np.concatenate([g0, g1]),
        preset=config.ic_preset,
        particle_mass=particle_mass,
    )


def load_into_particle_state(ic: ICResult, parts) -> None:
    """Copy an :class:`ICResult` into an allocated Taichi ``ParticleState``."""
    if parts.n != ic.n:
        raise ValueError(f"ParticleState has {parts.n} slots but IC has {ic.n} particles")
    parts.pos_x.from_numpy(np.ascontiguousarray(ic.pos[:, 0]))
    parts.pos_y.from_numpy(np.ascontiguousarray(ic.pos[:, 1]))
    parts.pos_z.from_numpy(np.ascontiguousarray(ic.pos[:, 2]))
    parts.vel_x.from_numpy(np.ascontiguousarray(ic.vel[:, 0]))
    parts.vel_y.from_numpy(np.ascontiguousarray(ic.vel[:, 1]))
    parts.vel_z.from_numpy(np.ascontiguousarray(ic.vel[:, 2]))
    parts.mass.from_numpy(np.ascontiguousarray(ic.mass))
    parts.gid.from_numpy(np.ascontiguousarray(ic.gid.astype(np.int32)))


__all__ = [
    "GalaxyModel",
    "ICResult",
    "build_ic",
    "load_into_particle_state",
    "V_APPROACH",
    "V_APPROACH_KMS",
]
