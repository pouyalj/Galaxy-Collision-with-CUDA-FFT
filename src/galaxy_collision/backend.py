"""Backend capability helpers (Stage 6 / RV15, D6/D22).

Some Taichi GPU archs — notably Apple **Metal**, and the SPIR-V family (Vulkan/OpenGL) —
have **no hardware fp64**: declaring an ``ti.f64`` field or doing f64 arithmetic in a
kernel fails to compile (Metal reports ``Type f64 not supported``). The device reductions
in :mod:`galaxy_collision.solver.multigrid` and :mod:`galaxy_collision.diagnostics_device`
therefore choose their accumulator path from :func:`supports_fp64`:

* **fp64 available (CPU/CUDA):** the original thread-local f64 reductions, unchanged.
* **no fp64 (Metal):** a **column-hybrid fp32** reduction — fp32 partials per ``(i,j)``
  column (or per grid-stride lane), each Kahan-compensated, finished in host fp64 (a tiny
  fixed copy, never a full-grid round-trip). Two things buy the accuracy here, and it's worth
  not conflating them: (1) the *column structure* shortens each fp32 sum from ~1.7e7 terms to
  ~n, and the fp64 cross-column combine is exact — this alone gets the moments to ~1e-8
  (their position factors are constant per column, so most components reduce to
  const·column-mass); (2) **Kahan** refines that and matters most for the grid-PE and
  L2-norm reductions, where ρ·Φ / the residual vary along the column and *don't* factor. The
  ~1e-3 baseline both techniques improve on is a single global fp32 accumulator — no hybrid,
  no Kahan (see ``docs/stage6_plan.md`` 6A).

This module deliberately reads the *currently initialized* arch, so it must be called
after ``ti.init`` (i.e. after :func:`galaxy_collision.sim.init_backend`).
"""

from __future__ import annotations

import taichi as ti

# Taichi arch *names* (``current_cfg().arch.name``) that lack hardware fp64. CPU archs
# (``x64``/``arm64``) and ``cuda`` support it; Metal and the SPIR-V backends do not.
_NO_FP64_ARCHS = frozenset({"metal", "vulkan", "opengl", "gles", "dx11"})


def supports_fp64() -> bool:
    """Whether the currently initialized Taichi arch has hardware fp64.

    Returns ``True`` (the safe default) if Taichi is not yet initialized or the arch is
    unrecognized — the f64 path is always correct where it compiles.
    """
    try:
        arch = ti.lang.impl.current_cfg().arch.name
    except Exception:
        return True
    return arch not in _NO_FP64_ARCHS


__all__ = ["supports_fp64"]
