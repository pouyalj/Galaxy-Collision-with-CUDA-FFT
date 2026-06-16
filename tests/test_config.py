"""Tests for the run-configuration schema and loader (Stage 0)."""

from __future__ import annotations

from pathlib import Path

import pytest

from galaxy_collision.config import (
    GALAXY_MASS_MSUN,
    ConfigError,
    SimConfig,
    dump_config,
    from_dict,
    load_config,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SMOKE_CONFIG = REPO_ROOT / "configs" / "smoke.yaml"


def test_smoke_config_loads():
    config = load_config(SMOKE_CONFIG)
    assert config.name == "hello-sim"
    assert config.backend == "cpu"
    assert config.steps == 1


def test_defaults_are_valid():
    # __post_init__ validates; constructing with defaults must not raise.
    SimConfig()


def test_round_trip(tmp_path):
    original = SimConfig(name="rt", backend="cuda", n_particles=1000, steps=5, dt=0.02)
    path = tmp_path / "rt.yaml"
    dump_config(original, path)
    reloaded = load_config(path)
    assert reloaded == original


def test_resolve_n_from_mass():
    config = SimConfig(particle_mass=4640.0)
    n = config.resolve_n_particles()
    assert n == round(2 * GALAXY_MASS_MSUN / 4640.0)
    # ~100M working-range target from AGENT.md §5.2.
    assert 90_000_000 < n < 110_000_000


def test_resolve_n_direct():
    assert SimConfig(n_particles=1234).resolve_n_particles() == 1234


def test_resolve_n_none():
    assert SimConfig().resolve_n_particles() is None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"backend": "gpu"},
        {"solver": "spectral"},
        {"ic_preset": "blobs"},
        {"n_particles": 10, "particle_mass": 5.0},
        {"n_particles": -1},
        {"particle_mass": 0.0},
        {"dt": 0.0},
        {"steps": -1},
        {"grid_size": 0},
        {"output_cadence": -1},
    ],
)
def test_invalid_configs_rejected(kwargs):
    with pytest.raises(ConfigError):
        SimConfig(**kwargs)


def test_unknown_keys_rejected():
    with pytest.raises(ConfigError):
        from_dict({"name": "x", "bogus_key": 1})


def test_missing_file():
    with pytest.raises(ConfigError):
        load_config("/no/such/config.yaml")
