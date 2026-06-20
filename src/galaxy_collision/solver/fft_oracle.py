"""Zero-padded isolated-Green's-function FFT solver — the validation oracle (D5).

This is the **spectrally exact open-boundary** Poisson solve (Hockney & Eastwood): the
potential is the discrete convolution of the density with the free-space Green's
function of the Laplacian. Doing the convolution by FFT would normally be *circular*
(periodic, the very artifact we reject); zero-padding the box to ``pad_factor·N`` per
axis turns it back into the *linear* convolution we want, so the result is the isolated
(Φ → 0 at infinity) potential with no periodic ghost images.

Role (AGENT.md §5.4, §6): it is the ground truth multigrid is checked against (test 2)
and the engine for reproducing the paper's spectral figures. On CPU it runs on NumPy's
pocketfft; on NVIDIA the same math backs onto cuFFT (Stage 5). It is *not* the
production solver — at 256³ the padded 512³ complex buffer is ~1 GB (AGENT.md §5.2).

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


class FFTPoissonSolver(PoissonSolver):
    """Open-BC Poisson solve by zero-padded FFT convolution (the oracle)."""

    def __init__(
        self,
        grid_size: int,
        dx: float = 1.0,
        pad_factor: int = 2,
        grav_constant: float | None = None,
    ):
        if pad_factor < 2:
            raise ValueError(f"pad_factor must be >= 2 to avoid circular wrap, got {pad_factor}")
        self.grid_size = grid_size
        self.dx = float(dx)
        self.pad = pad_factor * grid_size
        self.G = units.G if grav_constant is None else grav_constant
        self._greens_ft = self._build_greens_ft()

    def _build_greens_ft(self) -> np.ndarray:
        """FFT of the periodically-tiled free-space kernel g_k = −G·dx³ / r_k."""
        m, dx = self.pad, self.dx
        # Symmetric index distance so the circular convolution over `m` reproduces the
        # linear convolution: index i and m−i are the same distance from the origin.
        d = np.minimum(np.arange(m), m - np.arange(m)).astype(np.float64)
        dxx, dyy, dzz = np.meshgrid(d, d, d, indexing="ij")
        r = dx * np.sqrt(dxx**2 + dyy**2 + dzz**2)
        green = np.zeros((m, m, m), dtype=np.float64)
        nonzero = r > 0.0
        green[nonzero] = -self.G * dx**3 / r[nonzero]
        green[0, 0, 0] = 0.0  # self term (see module docstring)
        return np.fft.rfftn(green)

    def solve(self, rho, phi) -> None:
        n, m = self.grid_size, self.pad
        rho_np = rho.to_numpy().astype(np.float64)
        padded = np.zeros((m, m, m), dtype=np.float64)
        padded[:n, :n, :n] = rho_np
        axes = (0, 1, 2)
        spectrum = np.fft.rfftn(padded, axes=axes) * self._greens_ft
        conv = np.fft.irfftn(spectrum, s=(m, m, m), axes=axes)
        phi.from_numpy(np.ascontiguousarray(conv[:n, :n, :n]).astype(np.float32))


__all__ = ["FFTPoissonSolver"]
