"""Shared pytest configuration.

**Re-running the kernel suite on a GPU backend (Stage 5/6).** AGENT.md §7 requires that
Stages 5 and 6 re-pass the Stage 3–4 test suite on their backend. The physics/solver tests
pin ``ti.init(arch=ti.cpu)`` for portable CI; rather than thread a backend parameter through
all 29 call sites, set the env var ``GALAXY_TEST_ARCH=cuda`` (or ``metal``) and the autouse
fixture below redirects those inits onto that backend for the whole run:

    GALAXY_TEST_ARCH=cuda pytest -q

Tests whose *purpose* is a specific backend's reporting/fallback semantics — or that compare
two backends against each other — opt out with ``@pytest.mark.fixed_arch`` so their explicit
arch choices are honored.
"""

from __future__ import annotations

import os

import pytest

_OVERRIDE = os.environ.get("GALAXY_TEST_ARCH", "").lower()


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "fixed_arch: test pins its own Taichi arch; exempt from the GALAXY_TEST_ARCH override.",
    )


@pytest.fixture(autouse=True)
def _force_test_arch(request, monkeypatch):
    """Redirect ``ti.init(arch=…)`` onto ``GALAXY_TEST_ARCH`` when that env var selects a GPU.

    No-op when the var is unset/``cpu`` (the default CI path) or when the test is marked
    ``fixed_arch``. The patch rewrites only the ``arch`` kwarg, so ``random_seed`` /
    ``offline_cache`` and friends pass through unchanged; it patches the ``taichi`` module
    attribute, so both direct test inits and the ones inside ``sim.init_backend`` follow it.
    """
    if not _OVERRIDE or _OVERRIDE == "cpu":
        return
    if request.node.get_closest_marker("fixed_arch"):
        return
    import taichi as ti

    target = getattr(ti, _OVERRIDE)
    real_init = ti.init

    def patched(*args, **kwargs):
        kwargs.pop("arch", None)
        return real_init(*args, arch=target, **kwargs)

    monkeypatch.setattr(ti, "init", patched)
