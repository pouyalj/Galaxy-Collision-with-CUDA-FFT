"""Backend detection & hardware reporting (Stage 3 follow-up).

The headline safety property: Taichi *silently* falls back to CPU when a GPU backend is
unavailable (e.g. ``cuda`` on a Mac, reported as arch ``arm64``). ``_backend_honored`` must
flag that so a run never looks like it used a GPU it didn't. That logic is pure (no Taichi),
so it runs everywhere; the live ``init_backend`` device report self-skips without Taichi.
"""

from __future__ import annotations

import pytest

from galaxy_collision.sim import _backend_honored


@pytest.mark.parametrize(
    "requested,resolved,expected",
    [
        ("cpu", "arm64", True),  # CPU request honored by any native CPU arch
        ("cpu", "x64", True),
        ("cpu", "cpu", True),
        ("cuda", "cuda", True),
        ("metal", "metal", True),
        ("cuda", "arm64", False),  # the silent-fallback trap: asked CUDA, got CPU
        ("cuda", "x64", False),
        ("metal", "arm64", False),
        ("metal", "cuda", False),
    ],
)
def test_backend_honored(requested, resolved, expected):
    assert _backend_honored(requested, resolved) is expected


def test_init_backend_reports_cpu_device():
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision.config import SimConfig
    from galaxy_collision.sim import init_backend

    info = init_backend(SimConfig(backend="cpu"))
    assert info["requested"] == "cpu"
    assert info["resolved"] in {"arm64", "x64", "cpu"}
    assert info["fell_back"] is False  # CPU is always available
    assert info["device"] and "CPU" in info["device"]
