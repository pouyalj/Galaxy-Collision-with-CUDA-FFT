"""Zero-padded isolated-Green's-function FFT solver — the validation oracle (D5).

This is the **spectrally exact open-boundary** Poisson solve (Hockney & Eastwood): the
potential is the discrete convolution of the density with the free-space Green's
function of the Laplacian. Doing the convolution by FFT would normally be *circular*
(periodic, the very artifact we reject); zero-padding the box to ``pad_factor·N`` per
axis turns it back into the *linear* convolution we want, so the result is the isolated
(Φ → 0 at infinity) potential with no periodic ghost images.

Role (AGENT.md §5.4, §6): it is the ground truth multigrid is checked against (test 2)
and the engine for reproducing the paper's spectral figures. The FFT runs on **CuPy/cuFFT
when available** (Stage 5/5C, D21) and falls back to NumPy's pocketfft otherwise — the math
is identical, selected by the array module ``self._xp``. cuFFT makes the 512³ transform
~19× faster (≈0.1 s vs ≈2 s), so multigrid can be validated against the oracle at the full
production 256³ grid on the GPU, not just the small grids CI affords on CPU. It is *not* the
production solver: at 256³ the persistent Green's-function transform is ~1 GB, but the
**per-solve working set is ~3 GB** (zero-padded ρ + kernel-FT + spectrum, all 512³), and on GPU
CuPy's pool reserves ~5–6 GB once the cuFFT plan workspace is included — so the oracle runs
*alone* for validation, never co-resident with a 100M production run (AGENT.md §5.2).

**Derivation of the kernel.** ∇²Φ = 4πGρ with the free-space Laplacian Green's function
∇²(−1/(4π r)) = δ gives Φ(r) = −G ∫ ρ(r')/|r−r'| d³r'. Discretized on node-centered
cells of volume dx³ (mass m_j = ρ_j·dx³):

    Φ_i = Σ_j ρ_j · (−G·dx³ / r_ij),   so the convolution kernel is  g_k = −G·dx³ / r_k.

The self term (r=0) is set to 0: it shifts only the potential *at* a source node, never
the inter-node values the analytic −GM/r test checks.
"""

from __future__ import annotations

import numpy as np

from galaxy_collision import units
from galaxy_collision.solver.base import PoissonSolver

try:  # optional GPU FFT backend (D21); falls back to NumPy when absent
    import cupy as _cupy
except Exception:  # pragma: no cover - import guard
    _cupy = None


def _cupy_usable() -> bool:
    """True iff CuPy is importable and a CUDA device is present."""
    if _cupy is None:
        return False
    try:
        return _cupy.cuda.runtime.getDeviceCount() > 0
    except Exception:  # pragma: no cover - no driver / no device
        return False


class FFTPoissonSolver(PoissonSolver):
    """Open-BC Poisson solve by zero-padded FFT convolution (the oracle).

    The FFT runs on CuPy/cuFFT when ``use_gpu`` (default: auto-detect a CUDA device + CuPy),
    else NumPy. ``rho``/``phi`` still bridge through the host (Taichi field ↔ NumPy); only the
    transform itself is offloaded — fine for a validation-only solver.
    """

    def __init__(
        self,
        grid_size: int,
        dx: float = 1.0,
        pad_factor: int = 2,
        grav_constant: float | None = None,
        use_gpu: bool | None = None,
    ):
        if pad_factor < 2:
            raise ValueError(f"pad_factor must be >= 2 to avoid circular wrap, got {pad_factor}")
        self.grid_size = grid_size
        self.dx = float(dx)
        self.pad = pad_factor * grid_size
        self.G = units.G if grav_constant is None else grav_constant
        self.on_gpu = _cupy_usable() if use_gpu is None else bool(use_gpu)
        if self.on_gpu and _cupy is None:
            raise RuntimeError("use_gpu=True but CuPy is not installed")
        self._xp = _cupy if self.on_gpu else np
        self._greens_ft = self._build_greens_ft()

    def _build_greens_ft(self):
        """FFT of the periodically-tiled free-space kernel g_k = −G·dx³ / r_k (on ``self._xp``)."""
        xp, m, dx = self._xp, self.pad, self.dx
        # Symmetric index distance so the circular convolution over `m` reproduces the
        # linear convolution: index i and m−i are the same distance from the origin.
        d = xp.minimum(xp.arange(m), m - xp.arange(m)).astype(xp.float64)
        d2 = d * d
        # Broadcast the per-axis squared distances rather than meshgrid (which would
        # materialize three full (m,m,m) arrays — ~3.2 GB transient at m=512). The sum
        # broadcasts into a single (m,m,m) array.
        r2 = d2[:, None, None] + d2[None, :, None] + d2[None, None, :]
        r2[0, 0, 0] = 1.0  # avoid 0-division at the self term; overwritten below
        green = (-self.G * dx**3 / dx) / xp.sqrt(r2)  # g = −G·dx³ / (dx·|index|)
        green[0, 0, 0] = 0.0  # self term (see module docstring)
        return xp.fft.rfftn(green)

    def solve(self, rho, phi, warm_start: bool = False) -> None:
        # warm_start is ignored: this is a direct (non-iterative) spectral solve.
        xp, n, m = self._xp, self.grid_size, self.pad
        rho_np = rho.to_numpy().astype(np.float64)
        padded = xp.zeros((m, m, m), dtype=xp.float64)
        padded[:n, :n, :n] = xp.asarray(rho_np)  # host→device when on GPU; no-op on NumPy
        axes = (0, 1, 2)
        spectrum = xp.fft.rfftn(padded, axes=axes) * self._greens_ft
        conv = xp.fft.irfftn(spectrum, s=(m, m, m), axes=axes)
        out = conv[:n, :n, :n]
        out = _cupy.asnumpy(out) if self.on_gpu else out  # device→host for the Taichi field
        phi.from_numpy(np.ascontiguousarray(out).astype(np.float32))


__all__ = ["FFTPoissonSolver"]
