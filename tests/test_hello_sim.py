"""Stage-0 smoke test: the no-op simulation runs on the CPU backend.

Skipped cleanly if Taichi is not installed, so the config tests still run in a
minimal environment. CI installs Taichi, so this gates the green baseline there.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from galaxy_collision.config import load_config

SMOKE_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "smoke.yaml"

pytest.importorskip("taichi", reason="Taichi not installed")


def test_hello_sim_runs_from_smoke_config():
    from galaxy_collision.sim import run_hello_sim

    config = load_config(SMOKE_CONFIG)
    config.backend = "cpu"  # force CPU regardless of available accelerators
    result = run_hello_sim(config)

    assert result["status"] == "ok"
    assert result["backend"] == "cpu"
    assert result["steps"] == config.steps
    assert result["n_particles"] >= 1
