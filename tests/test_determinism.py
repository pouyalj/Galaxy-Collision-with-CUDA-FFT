"""Cross-backend determinism — CPU vs CUDA (Stage 5, partial AGENT.md §6 test 6).

§5.3 requires that an identical IC + seed produce matching trajectories across backends
within an fp32 tolerance. Stage 5 can check the CPU↔CUDA half of that now (the full 3-way
parity including Metal is Stage 6). The check is a Stage-5 guard that the CUDA force chain
(device-resident deposit → multigrid → grad → gather → KDK) computes the *same physics* as
the trusted CPU reference, not a silently different kernel.

This test self-skips when no CUDA device is available, so it is a no-op on CPU-only CI and a
real check on the GPU workstation.
"""

from __future__ import annotations

import numpy as np
import pytest


def _cuda_available() -> bool:
    pytest.importorskip("taichi", reason="Taichi not installed")
    from galaxy_collision.config import SimConfig
    from galaxy_collision.sim import init_backend

    return not init_backend(SimConfig(backend="cuda"))["fell_back"]


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

    Observed on the RTX 3070 dev box: the fp32 tracer paths are *bit-identical* to the CPU
    run over this horizon (Taichi's CIC scatter + multigrid reductions are deterministic
    here), so the 1e-3 kpc cap is loose headroom for fp32 reduction-order differences — a
    genuinely divergent CUDA kernel would miss by kpc, not microns.
    """
    if not _cuda_available():
        pytest.skip("no CUDA device available")

    steps = 15
    rc, rg = _run("cpu", steps), _run("cuda", steps)

    # Sanity: each run actually used the backend it claimed (no silent CPU fallback).
    assert rg["backend_resolved"] == "cuda"
    assert rc["backend_resolved"] in {"x64", "arm64", "cpu"}

    dpos = np.abs(rc["tracer_path"] - rg["tracer_path"])  # (samples, n_tracers, 3)
    assert dpos.max() < 1e-3  # kpc
