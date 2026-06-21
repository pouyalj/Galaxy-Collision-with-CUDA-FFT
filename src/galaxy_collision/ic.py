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

# Particle component tags (ICResult.component), used by the disk-velocity equilibration
# (§5.5 / D18) and the Stage-4 tracer selection.
COMPONENT_DISK, COMPONENT_BULGE, COMPONENT_BH = 0, 1, 2


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
    # (N,) int8 component tag (disk/bulge/BH); None for presets without that structure
    # (e.g. the single-sphere Plummer IC). See COMPONENT_* and equilibrate_disk_velocities.
    component: np.ndarray | None = None

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
    component = np.concatenate(
        [np.full(n_disk, COMPONENT_DISK), np.full(n_bulge, COMPONENT_BULGE)]
    )

    # Central black hole: one massive, softened particle at the center, at rest
    # relative to the galaxy (bulk velocity added below).
    pos = np.concatenate([pos, np.zeros((1, 3))], axis=0)
    vel = np.concatenate([vel, np.zeros((1, 3))], axis=0)
    mass = np.concatenate([mass, [model.bh_mass]])
    component = np.concatenate([component, [COMPONENT_BH]]).astype(np.int8)

    # Translate to the galaxy center and add the bulk approach velocity.
    pos += np.asarray(model.center)
    vel += np.asarray(model.bulk_velocity)
    gid_arr = np.full(pos.shape[0], gid, dtype=np.int32)
    return pos, vel, mass, gid_arr, component


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
    p0, v0, m0, g0, c0 = _build_galaxy(rng, 0, model0, particle_mass, n_per_galaxy)
    p1, v1, m1, g1, c1 = _build_galaxy(rng, 1, model1, particle_mass, n_per_galaxy)

    return ICResult(
        pos=np.concatenate([p0, p1], axis=0).astype(np.float32),
        vel=np.concatenate([v0, v1], axis=0).astype(np.float32),
        mass=np.concatenate([m0, m1]).astype(np.float32),
        gid=np.concatenate([g0, g1]),
        preset=config.ic_preset,
        particle_mass=particle_mass,
        component=np.concatenate([c0, c1]).astype(np.int8),
    )


# --- Plummer sphere (equilibrium test IC, AGENT.md §6 test 3) -------------------


@dataclass(frozen=True)
class PlummerModel:
    """An isolated Plummer sphere: ρ(r) = (3M/4πa³)(1 + r²/a²)^(−5/2)."""

    total_mass: float = 1.0e10  # M_sun
    scale_a: float = 5.0  # kpc (Plummer scale radius)
    center: tuple[float, float, float] = (128.0, 128.0, 128.0)


def _sample_plummer_q(rng, n: int) -> np.ndarray:
    """Sample the dimensionless speed ratio q=v/v_esc from g(q)=q²(1−q²)^{7/2} by rejection.

    g peaks at ≈0.092 (q²=2/9), so the constant 0.1 envelope (Aarseth–Hénon–Wielen 1974)
    bounds it. Vectorized: oversample, keep accepted, repeat until n are filled.
    """
    out = np.empty(n, dtype=np.float64)
    filled = 0
    while filled < n:
        m = n - filled
        x = rng.random(m)
        y = rng.random(m)
        accept = 0.1 * y < x**2 * (1.0 - x**2) ** 3.5
        k = int(np.count_nonzero(accept))
        out[filled : filled + k] = x[accept]
        filled += k
    return out


def _sample_plummer(rng, n: int, total_mass: float, scale_a: float, grav: float):
    """Sample a Plummer sphere in equilibrium: positions + isotropic velocities (at origin)."""
    # Positions: invert M(<r)/M = r³/(r²+a²)^{3/2} = u  →  r = a / √(u^{−2/3} − 1).
    # Clamp u away from 0 so u^{−2/3} stays finite (u=0 would give r=0 anyway).
    u = np.clip(rng.random(n), 1e-12, 1.0)
    radius = scale_a / np.sqrt(u ** (-2.0 / 3.0) - 1.0)
    radius = np.minimum(radius, 20.0 * scale_a)  # clip the rare far tail
    pdir = rng.normal(size=(n, 3))
    pdir /= np.linalg.norm(pdir, axis=1, keepdims=True)
    pos = radius[:, None] * pdir

    # Velocities: speed = q · v_esc(r), v_esc = √(2GM) (r²+a²)^{−1/4}; isotropic direction.
    q = _sample_plummer_q(rng, n)
    v_esc = np.sqrt(2.0 * grav * total_mass) * (radius**2 + scale_a**2) ** (-0.25)
    vdir = rng.normal(size=(n, 3))
    vdir /= np.linalg.norm(vdir, axis=1, keepdims=True)
    vel = (q * v_esc)[:, None] * vdir

    # Recenter: with finite N the sampled sphere carries a small net CM offset and bulk
    # momentum. Subtract both so the sphere sits at the origin and is genuinely at rest —
    # otherwise it drifts across the box and total momentum is a nonzero sampling artifact.
    pos -= pos.mean(axis=0)
    vel -= vel.mean(axis=0)
    return pos, vel


def build_plummer_ic(
    config: SimConfig,
    model: PlummerModel | None = None,
    grav: float | None = None,
) -> ICResult:
    """Build an isolated Plummer-sphere IC (for the equilibrium-stability exit-gate test).

    A single sphere of equal-mass particles at rest at the box center, sampled in virial
    equilibrium with the analytic Plummer distribution function. The particle count is taken
    from ``config.n_particles`` if set, else derived from ``config.particle_mass`` against
    *this sphere's* ``total_mass`` (not the two-galaxy mass that ``resolve_n_particles``
    assumes), else defaults to 10,000.
    """
    if config.ic_preset != "plummer":
        raise NotImplementedError(
            f"build_plummer_ic requires ic_preset='plummer', got {config.ic_preset!r}"
        )
    m = model or PlummerModel(center=(config.grid_size / 2.0,) * 3)
    g = units.G if grav is None else grav
    # Resolve N against the Plummer total mass — do NOT use config.resolve_n_particles(),
    # which derives N from the two-galaxy GALAXY_MASS_MSUN and is wrong for one sphere.
    if config.n_particles is not None:
        n = config.n_particles
    elif config.particle_mass is not None:
        n = round(m.total_mass / config.particle_mass)
    else:
        n = 10_000

    rng = np.random.default_rng(config.seed)
    pos, vel = _sample_plummer(rng, n, m.total_mass, m.scale_a, g)
    pos = pos + np.asarray(m.center)
    particle_mass = m.total_mass / n
    return ICResult(
        pos=pos.astype(np.float32),
        vel=vel.astype(np.float32),
        mass=np.full(n, particle_mass, dtype=np.float32),
        gid=np.zeros(n, dtype=np.int32),
        preset="plummer",
        particle_mass=particle_mass,
    )


# --- Disk-velocity equilibration on the PM grid (§5.5 / D18) --------------------


class _RawParticles:
    """Minimal duck-typed particle holder for the deposit/gather kernels (pos + opt. mass).

    Lazily imports Taichi so ``ic`` stays importable (and the pure-NumPy IC tests keep
    running) without a Taichi runtime.
    """

    def __init__(self, pos: np.ndarray, mass: np.ndarray | None = None):
        import taichi as ti

        self.n = int(pos.shape[0])
        self.pos_x = ti.field(ti.f32, shape=self.n)
        self.pos_y = ti.field(ti.f32, shape=self.n)
        self.pos_z = ti.field(ti.f32, shape=self.n)
        self.pos_x.from_numpy(np.ascontiguousarray(pos[:, 0], dtype=np.float32))
        self.pos_y.from_numpy(np.ascontiguousarray(pos[:, 1], dtype=np.float32))
        self.pos_z.from_numpy(np.ascontiguousarray(pos[:, 2], dtype=np.float32))
        if mass is not None:
            self.mass = ti.field(ti.f32, shape=self.n)
            self.mass.from_numpy(np.ascontiguousarray(mass, dtype=np.float32))


def equilibrate_disk_velocities(
    icr: ICResult,
    *,
    dx: float,
    solver,
    rho,
    phi,
    ax,
    ay,
    az,
    n_radii: int = 64,
    n_azimuth: int = 64,
    parts=None,
) -> dict:
    """Launch each disk in equilibrium with the PM grid force it will actually feel (§5.5/D18).

    The analytic IC rotation curve (:func:`_circular_velocity`) uses a 0.3 kpc Plummer
    softening, but the deposit→solve→grad→gather chain softens the force at the ~1 kpc cell
    scale — so a disk built from the analytic curve is launched too fast in its inner few kpc
    and "breathes"/heats. Here we instead *measure* the curve the grid produces and launch on
    it:

    For each galaxy, its own particles (only) are deposited on the production grid, solved with
    the production ``solver``, and differentiated to the acceleration field. A ring of probes in
    the disk midplane samples the azimuthally-averaged inward radial acceleration ⟨−a_R⟩(R), so
    ``v_c(R) = √(R·⟨−a_R⟩)`` is exactly the circular speed the grid supports. Each disk particle's
    in-plane velocity is set to ``v_c(R)`` (its galaxy's bulk approach velocity preserved). Bulge
    (pressure-supported) and BH velocities are left untouched.

    This balances each disk against its galaxy's **self-gravity only** — the companion's tidal
    field is neglected. That is the right approximation at the wide initial separation (~90 kpc),
    but it is an *initial-condition* equilibrium: it does **not** hold once the galaxies interact
    (through pericenter the tidal field is large and the disk is *supposed* to be driven out of
    equilibrium — that is the physics we are simulating, not an IC error).

    Because the *same* solver/grid/dx are used here and in the run, the disk is in equilibrium
    with the run's (self-gravity) force to within the solver's own convergence — no
    analytic-softening guess.
    ``rho``/``phi``/``ax``/``ay``/``az`` are caller-owned scratch grid fields (the production grid
    fields work — the main loop re-primes them). Mutates ``icr.vel`` in place and, if ``parts`` is
    given, writes the new disk velocities there too. Returns ``{gid: (r_table, vc_table)}``.
    """
    import taichi as ti

    from galaxy_collision.deposit import deposit_density, gather_acceleration, potential_to_accel

    if icr.component is None:
        raise ValueError("equilibrate_disk_velocities needs ICResult.component (two-galaxy IC)")
    # G is carried by `solver` (its grav_constant) — used implicitly via the solve below.
    pos = icr.pos.astype(np.float64)
    vel = icr.vel.astype(np.float64)
    curves: dict = {}

    for gid in np.unique(icr.gid):
        sel = icr.gid == gid
        idx_g = np.where(sel)[0]
        comp_g = icr.component[sel]
        pos_g = pos[idx_g]

        bh = comp_g == COMPONENT_BH
        if bh.any():
            center = pos_g[bh][0]
            bulk = vel[idx_g[bh][0]]
        else:  # defensive: no BH — fall back to the stellar centroid, at rest
            center = pos_g.mean(axis=0)
            bulk = np.zeros(3)

        disk = comp_g == COMPONENT_DISK
        if not disk.any():
            continue
        dxp = pos_g[disk, 0] - center[0]
        dyp = pos_g[disk, 1] - center[1]
        r_disk = np.hypot(dxp, dyp)
        r_max = float(r_disk.max())

        # Galaxy self-gravity potential → acceleration field on the production grid.
        deposit_density(_RawParticles(pos_g, icr.mass[idx_g]), rho, dx)
        solver.solve(rho, phi, warm_start=False)
        potential_to_accel(phi, ax, ay, az, dx)

        # Probe rings in the midplane; azimuthally average the inward radial acceleration.
        r_tab = np.linspace(dx, r_max, n_radii)
        theta = np.linspace(0.0, 2.0 * np.pi, n_azimuth, endpoint=False)
        rr, tt = np.meshgrid(r_tab, theta, indexing="ij")  # (n_radii, n_azimuth)
        cos_t, sin_t = np.cos(tt), np.sin(tt)
        px = (center[0] + rr * cos_t).ravel()
        py = (center[1] + rr * sin_t).ravel()
        pz = np.full(px.size, center[2])
        probes = _RawParticles(np.stack([px, py, pz], axis=1))
        pax = ti.field(ti.f32, shape=probes.n)
        pay = ti.field(ti.f32, shape=probes.n)
        paz = ti.field(ti.f32, shape=probes.n)
        gather_acceleration(probes, ax, ay, az, pax, pay, paz, dx)
        a_rad = pax.to_numpy().reshape(n_radii, n_azimuth) * cos_t + \
            pay.to_numpy().reshape(n_radii, n_azimuth) * sin_t
        # Azimuthally average the radial acceleration, *then* clamp the mean. (Clamping each
        # sample first, np.maximum(−a_rad,0).mean, would bias ⟨−a_R⟩ upward on a ring that is
        # net-inward but has a few noisy/outward azimuths.)
        inward = np.maximum(-a_rad.mean(axis=1), 0.0)  # ⟨−a_R⟩(R)
        vc_tab = np.sqrt(r_tab * inward)
        curves[int(gid)] = (r_tab, vc_tab)

        # Assign disk velocities: v_c(R) tangential + the bulk approach. Interpolate the angular
        # speed ω = v_c/R (smooth and ≈constant through the solid-body core) rather than v_c, so a
        # particle inside the innermost ring (R < dx) gets v_c = R·ω → 0 as R → 0, instead of
        # np.interp clamping it to v_c(dx) and launching the sub-kpc disk over-rotating.
        omega_disk = np.interp(r_disk, r_tab, vc_tab / r_tab)
        vc_disk = r_disk * omega_disk
        safe = np.maximum(r_disk, 1e-6)
        tang = np.stack([-dyp / safe, dxp / safe, np.zeros_like(r_disk)], axis=1)
        vel[idx_g[disk]] = vc_disk[:, None] * tang + bulk[None, :]

    icr.vel = vel.astype(np.float32)
    if parts is not None:
        parts.vel_x.from_numpy(np.ascontiguousarray(icr.vel[:, 0]))
        parts.vel_y.from_numpy(np.ascontiguousarray(icr.vel[:, 1]))
        parts.vel_z.from_numpy(np.ascontiguousarray(icr.vel[:, 2]))
    return curves


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
    "PlummerModel",
    "build_plummer_ic",
    "load_into_particle_state",
    "equilibrate_disk_velocities",
    "COMPONENT_DISK",
    "COMPONENT_BULGE",
    "COMPONENT_BH",
    "V_APPROACH",
    "V_APPROACH_KMS",
]
