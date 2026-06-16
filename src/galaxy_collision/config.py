"""Run-configuration schema and YAML loader.

A :class:`SimConfig` is the single source of truth for a run: backend, particle
budget (set N directly *or* derive it from a per-particle mass, per AGENT.md D8),
timestep, step count, initial-condition preset, Poisson solver, and output cadence.

Stage 0 only needs the schema to load, validate, and round-trip cleanly through
YAML. Later stages consume these fields without reparsing.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any

import yaml

# --- Allowed enumerations -------------------------------------------------------

BACKENDS = ("cpu", "cuda", "metal")
SOLVERS = ("multigrid", "fft")
IC_PRESETS = ("hello", "two_galaxy_4v", "two_galaxy_2v", "plummer")

# Paper value: mass within 25 kpc of each galaxy (M_sun). Used to derive N from a
# per-particle mass (D8). See AGENT.md §5.2.
GALAXY_MASS_MSUN = 2.32e11


class ConfigError(ValueError):
    """Raised when a configuration is structurally or semantically invalid."""


@dataclass
class SimConfig:
    """A single simulation run configuration.

    Exactly one of ``n_particles`` or ``particle_mass`` may be given as the
    resolution knob; if both are omitted the IC preset's own default is used.
    """

    name: str = "hello-sim"
    backend: str = "cpu"

    # Resolution knob (mutually exclusive). particle_mass is M_sun per particle.
    n_particles: int | None = None
    particle_mass: float | None = None

    # Integration.
    dt: float = 0.01  # Myr (the 2020 code's dt=1 == 10^4 yr == 0.01 Myr)
    steps: int = 1

    # Physics setup.
    ic_preset: str = "hello"
    solver: str = "multigrid"
    grid_size: int = 256  # cells per side; box is grid_size kpc, 1 kpc^3 cells

    # Output.
    output_cadence: int = 0  # steps between snapshots; 0 disables output
    output_dir: str = "outputs"

    # Reproducibility.
    seed: int = 0

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.backend not in BACKENDS:
            raise ConfigError(f"backend must be one of {BACKENDS}, got {self.backend!r}")
        if self.solver not in SOLVERS:
            raise ConfigError(f"solver must be one of {SOLVERS}, got {self.solver!r}")
        if self.ic_preset not in IC_PRESETS:
            raise ConfigError(f"ic_preset must be one of {IC_PRESETS}, got {self.ic_preset!r}")
        if self.n_particles is not None and self.particle_mass is not None:
            raise ConfigError(
                "set only one resolution knob: n_particles OR particle_mass, not both"
            )
        if self.n_particles is not None and self.n_particles <= 0:
            raise ConfigError(f"n_particles must be positive, got {self.n_particles}")
        if self.particle_mass is not None and self.particle_mass <= 0:
            raise ConfigError(f"particle_mass must be positive, got {self.particle_mass}")
        if self.dt <= 0:
            raise ConfigError(f"dt must be positive, got {self.dt}")
        if self.steps < 0:
            raise ConfigError(f"steps must be non-negative, got {self.steps}")
        if self.grid_size <= 0:
            raise ConfigError(f"grid_size must be positive, got {self.grid_size}")
        if self.output_cadence < 0:
            raise ConfigError(f"output_cadence must be non-negative, got {self.output_cadence}")

    def resolve_n_particles(self) -> int | None:
        """Return the particle count, deriving it from particle_mass if needed (D8).

        N = (mass of both galaxies) / m_particle. Returns ``None`` when neither knob
        is set (the IC preset supplies its own default in a later stage).
        """
        if self.n_particles is not None:
            return self.n_particles
        if self.particle_mass is not None:
            return round(2 * GALAXY_MASS_MSUN / self.particle_mass)
        return None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def from_dict(data: dict[str, Any]) -> SimConfig:
    """Build a :class:`SimConfig` from a plain dict, rejecting unknown keys."""
    if not isinstance(data, dict):
        raise ConfigError(f"config must be a mapping, got {type(data).__name__}")
    known = {f.name for f in fields(SimConfig)}
    unknown = set(data) - known
    if unknown:
        raise ConfigError(f"unknown config keys: {sorted(unknown)}")
    return SimConfig(**data)


def load_config(path: str | Path) -> SimConfig:
    """Load and validate a YAML run configuration from ``path``."""
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"config file not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    return from_dict(data)


def dump_config(config: SimConfig, path: str | Path) -> None:
    """Write a config to YAML (round-trips through :func:`load_config`)."""
    path = Path(path)
    with path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(config.to_dict(), fh, sort_keys=False)


__all__ = [
    "SimConfig",
    "ConfigError",
    "load_config",
    "dump_config",
    "from_dict",
    "BACKENDS",
    "SOLVERS",
    "IC_PRESETS",
]
