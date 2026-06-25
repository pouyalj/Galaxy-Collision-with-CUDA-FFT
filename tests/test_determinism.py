"""Cross-backend determinism — CPU vs CUDA vs Metal (AGENT.md §6 test 6).

§5.3 requires that an identical IC + seed produce matching trajectories across backends
within an fp32 tolerance. Stage 5 added the CPU↔CUDA leg; Stage 6 completes the 3-way parity
with CPU↔Metal (RV15 made the device reductions Metal-safe). Each leg guards that a GPU force
chain (device-resident deposit → multigrid → grad → gather → KDK) computes the *same physics*
as the trusted CPU reference, not a silently different kernel.

Each GPU leg self-skips when its device is unavailable, so on CPU-only CI both are no-ops; on
the NVIDIA workstation the CUDA leg runs, and on the Apple dev box the Metal leg runs.
"""

from __future__ import annotations

import numpy as np
import pytest


def _backend_available(backend: str) -> bool:
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision.config import SimConfig
    from galaxy_collision.sim import init_backend

    return not init_backend(SimConfig(backend=backend))["fell_back"]


def _run(backend: str, steps: int):
    from galaxy_collision.config import SimConfig
    from galaxy_collision.sim import run_simulation

    cfg = SimConfig(
        name="determinism", backend=backend, ic_preset="plummer", n_particles=4000,
        dt=0.5, steps=steps, solver="multigrid", grid_size=64, output_cadence=0, seed=3,
    )
    return run_simulation(
        cfg, history_cadence=1, write_snapshots=False, tracer_indices=[0, 137, 999, 2500]
    )


@pytest.mark.fixed_arch  # explicitly compares cpu vs cuda — must not be redirected to one arch
def test_cpu_cuda_trajectories_agree():
    """Same IC+seed on CPU and CUDA: tracer trajectories agree to fp32 tolerance.

    The guarantee asserted here is *physics-level* agreement — ‖Δpos‖ ≤ 1e-3 kpc over a short,
    pre-chaotic horizon — not bitwise equality. We happened to observe exact (Δ=0) agreement on
    this RTX 3070 / x86 pair at this scale, but that is **not** relied on: parallel atomic deposit
    + tree reductions have a non-deterministic summation order, and CPU vs GPU sqrt/div/FMA can
    round differently, so exact agreement will generally break at larger N, on other hardware, or
    across Taichi versions. The 1e-3 kpc cap is deliberately loose headroom for that fp32 round-off
    — a genuinely divergent CUDA kernel would miss by kpc, not microns; the short horizon keeps
    chaotic N-body amplification from inflating round-off.
    """
    if not _backend_available("cuda"):
        pytest.skip("no CUDA device available")

    steps = 15
    rc, rg = _run("cpu", steps), _run("cuda", steps)

    # Sanity: each run actually used the backend it claimed (no silent CPU fallback).
    assert rg["backend_resolved"] == "cuda"
    assert rc["backend_resolved"] in {"x64", "arm64", "cpu"}

    dpos = np.abs(rc["tracer_path"] - rg["tracer_path"])  # (samples, n_tracers, 3)
    assert dpos.max() < 1e-3  # kpc


@pytest.mark.fixed_arch  # explicitly compares cpu vs metal — must not be redirected to one arch
def test_cpu_metal_trajectories_agree():
    """Same IC+seed on CPU and Metal: tracer trajectories agree to fp32 tolerance (Stage 6).

    The Metal leg of test 6. Like the CUDA leg, this asserts *physics-level* agreement — not
    bitwise equality — over a short, pre-chaotic horizon. Metal additionally runs the
    Kahan-compensated fp32 boundary-moment reduction (RV15/D22) in place of CPU's f64 one, so a
    small extra round-off difference is expected; measured CPU↔Metal divergence here is ~2e-6 kpc,
    far under the 1e-3 kpc cap (the same loose headroom the CUDA leg uses — a genuinely divergent
    kernel would miss by kpc, not microns).
    """
    if not _backend_available("metal"):
        pytest.skip("no Metal device available")

    steps = 15
    rc, rm = _run("cpu", steps), _run("metal", steps)

    # Sanity: each run actually used the backend it claimed (no silent CPU fallback).
    assert rm["backend_resolved"] == "metal"
    assert rc["backend_resolved"] in {"x64", "arm64", "cpu"}

    dpos = np.abs(rc["tracer_path"] - rm["tracer_path"])  # (samples, n_tracers, 3)
    assert dpos.max() < 1e-3  # kpc
