"""Tests for the SoA data model and memory estimator (Stage 1).

The memory-estimator tests are pure Python. The allocation tests need Taichi and
self-skip without it. The 100M-particle allocation (a Stage-1 exit criterion) scales
down to fit available RAM so it runs everywhere, logging the actual count used.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from galaxy_collision.data import BYTES_PER_PARTICLE, estimate_memory


def _total_ram_bytes() -> int | None:
    """Best-effort total physical RAM in bytes, or None if undeterminable."""
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, AttributeError, OSError):
        pass
    try:  # macOS
        out = subprocess.run(
            ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=True
        )
        return int(out.stdout.strip())
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


# --- Memory estimator (pure Python) ---------------------------------------------


def test_bytes_per_particle():
    # 6 f32 (pos+vel) + 1 f32 (mass) + 1 i32 (gid) = 32 bytes.
    assert BYTES_PER_PARTICLE == 32


def test_estimate_100m_matches_agent_md():
    est = estimate_memory(100_000_000, grid_size=256)
    # AGENT.md §5.2: ~2.8 GB particles (at 28 B) → 3.2 GB at our 32 B incl. mass.
    assert 3.0 < est.particle_bytes / 1e9 < 3.4
    # Grid: 2 × 256^3 × 4 B ≈ 0.134 GB.
    assert 0.12 < est.grid_bytes / 1e9 < 0.15
    assert "GB" in est.summary()


def test_estimate_scales_linearly():
    a = estimate_memory(10_000_000)
    b = estimate_memory(20_000_000)
    assert b.particle_bytes == 2 * a.particle_bytes


@pytest.mark.parametrize("bad", [-1])
def test_estimate_rejects_negative(bad):
    with pytest.raises(ValueError):
        estimate_memory(bad)


def test_estimate_rejects_bad_grid():
    with pytest.raises(ValueError):
        estimate_memory(1000, grid_size=0)


# --- Taichi allocation (needs the runtime) --------------------------------------

pytest.importorskip("taichi", reason="Taichi not installed")


def test_particle_and_grid_allocate_and_free():
    import taichi as ti

    from galaxy_collision.data import GridState, ParticleState

    ti.init(arch=ti.cpu)
    parts = ParticleState(1024)
    grid = GridState(64)

    @ti.kernel
    def fill():
        for i in parts.pos_x:
            parts.pos_x[i] = ti.f32(i)
            parts.mass[i] = 840.0
            parts.gid[i] = i % 2

    fill()
    ti.sync()
    assert parts.pos_x[1023] == 1023.0
    assert parts.mass[0] == 840.0
    assert grid.rho.shape == (64, 64, 64)

    ti.reset()  # frees all fields; must not raise


def test_allocate_100m_particles_scaled_to_ram():
    import taichi as ti

    from galaxy_collision.data import ParticleState

    target = 100_000_000
    total = _total_ram_bytes()
    n = target
    if total is not None:
        # Keep the SoA under ~35% of physical RAM to leave headroom for Taichi/OS.
        max_fit = int(0.35 * total / BYTES_PER_PARTICLE)
        n = min(target, max_fit)
    n = max(n, 1_000_000)

    ti.init(arch=ti.cpu)
    parts = ParticleState(n)

    @ti.kernel
    def touch():
        parts.pos_x[0] = 1.0
        parts.pos_x[parts.pos_x.shape[0] - 1] = 2.0

    touch()
    ti.sync()
    assert parts.pos_x[0] == 1.0
    assert parts.pos_x[n - 1] == 2.0
    print(f"[stage1] allocated {n:,} particles "
          f"(~{estimate_memory(n).particle_bytes / 1e9:.2f} GB); target was {target:,}")

    ti.reset()  # frees cleanly
