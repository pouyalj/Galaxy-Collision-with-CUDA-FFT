"""Open-boundary geometric multigrid Poisson solver — the portable default (D5).

Solves ∇²Φ = 4πGρ on the node-centered grid with **open boundary conditions**,
implemented as real-space V-cycles in Taichi kernels (one source → CPU/CUDA/Metal,
no FFT dependency). This is what runs in production on every backend; the FFT oracle
exists only to validate it.

**Where correctness comes from.** The converged Φ is fixed entirely by three things on
the *finest* grid: the 7-point discrete Laplacian, the residual it defines, and the
Dirichlet values on the six box faces. Those faces are set once from a **multipole
expansion** (monopole + quadrupole about the mass's center of mass) of the enclosed
mass — that is what imposes "Φ → 0 far away" and fixes the legacy periodic-BC artifact
(bug #6). The multigrid hierarchy below is *only* an accelerator: it transfers
*corrections* with homogeneous BCs, so its geometric details (coarsening rule, transfer
stencils) affect how fast we converge, never what we converge to.

**Coarsening.** Vertex-centered, factor 2: coarse node ``J`` ↔ fine node ``2J``, so a
level of ``n`` nodes coarsens to ``n//2``. The coarse grid's far face therefore lands a
half-cell inside the fine one; harmless, per the paragraph above (the error there is
~0 under homogeneous Dirichlet BCs). Smoother: red-black Gauss-Seidel. Restriction:
full weighting. Prolongation: trilinear.

The boundary moments (M, center of mass, quadrupole) are reduced **on the device** by
``_accumulate_moments`` (Stage 5 / RV6a): a single pass over ``rho`` accumulates the ten
raw sums (M, Σm·r, Σm·rₐr_b) into a tiny f64 buffer, and only those 10 scalars cross to
the host for the cheap center-of-mass / second-moment arithmetic. This removes the
per-step full-grid ``rho.to_numpy()`` copy the Stage-3 reference did (legacy bug #13).
The f64 accumulation is CUDA/CPU-only; a Metal-safe fp32/Kahan reduction is a Stage-6
concern (Metal has no hardware fp64). The pure-NumPy ``_moments`` is retained as the
reference the device path is unit-tested against.
"""

import numpy as np
import taichi as ti

from galaxy_collision import units
from galaxy_collision.solver.base import PoissonSolver

# See deposit.py: PEP 563 would stringize the ti.template() annotations below.

_FOUR_PI = 4.0 * np.pi


@ti.kernel
def _gauss_seidel_rb(phi: ti.template(), rhs: ti.template(), n: ti.i32, h2: ti.f32, color: ti.i32):
    """One red-black Gauss-Seidel sweep over interior nodes of one parity ``color``."""
    for i, j, k in ti.ndrange((1, n - 1), (1, n - 1), (1, n - 1)):
        if (i + j + k) & 1 == color:
            phi[i, j, k] = (
                phi[i + 1, j, k]
                + phi[i - 1, j, k]
                + phi[i, j + 1, k]
                + phi[i, j - 1, k]
                + phi[i, j, k + 1]
                + phi[i, j, k - 1]
                - h2 * rhs[i, j, k]
            ) / 6.0


@ti.kernel
def _residual(
    phi: ti.template(), rhs: ti.template(), res: ti.template(), n: ti.i32, inv_h2: ti.f32
):
    """res = rhs − ∇²Φ on the interior; zero on the boundary."""
    for idx in ti.grouped(res):
        res[idx] = 0.0
    for i, j, k in ti.ndrange((1, n - 1), (1, n - 1), (1, n - 1)):
        lap = (
            phi[i + 1, j, k]
            + phi[i - 1, j, k]
            + phi[i, j + 1, k]
            + phi[i, j - 1, k]
            + phi[i, j, k + 1]
            + phi[i, j, k - 1]
            - 6.0 * phi[i, j, k]
        ) * inv_h2
        res[i, j, k] = rhs[i, j, k] - lap


@ti.kernel
def _restrict_full_weight(res_f: ti.template(), rhs_c: ti.template(), nc: ti.i32):
    """Full-weighting restriction of the fine residual onto the coarse RHS."""
    for idx in ti.grouped(rhs_c):
        rhs_c[idx] = 0.0
    for jc, kc, lc in ti.ndrange((1, nc - 1), (1, nc - 1), (1, nc - 1)):
        fi, fj, fk = 2 * jc, 2 * kc, 2 * lc
        acc = 0.0
        for a in ti.static(range(-1, 2)):
            for b in ti.static(range(-1, 2)):
                for c in ti.static(range(-1, 2)):
                    wa = 2.0 if a == 0 else 1.0
                    wb = 2.0 if b == 0 else 1.0
                    wc = 2.0 if c == 0 else 1.0
                    acc += wa * wb * wc * res_f[fi + a, fj + b, fk + c]
        rhs_c[jc, kc, lc] = acc / 64.0


@ti.kernel
def _prolong_add(phi_c: ti.template(), phi_f: ti.template(), nf: ti.i32):
    """Trilinear prolongation of the coarse correction, added into the fine interior."""
    for i, j, k in ti.ndrange((1, nf - 1), (1, nf - 1), (1, nf - 1)):
        ci, cj, ck = i // 2, j // 2, k // 2
        oi, oj, ok = i & 1, j & 1, k & 1
        val = 0.0
        # Trilinear weight per axis: an even fine node coincides with one coarse node
        # (weight 1); an odd fine node sits halfway between two (weight 0.5 each).
        for di in ti.static(range(2)):
            wi = 0.0
            if oi == 1:
                wi = 0.5
            elif ti.static(di == 0):
                wi = 1.0
            for dj in ti.static(range(2)):
                wj = 0.0
                if oj == 1:
                    wj = 0.5
                elif ti.static(dj == 0):
                    wj = 1.0
                for dk in ti.static(range(2)):
                    wk = 0.0
                    if ok == 1:
                        wk = 0.5
                    elif ti.static(dk == 0):
                        wk = 1.0
                    w = wi * wj * wk
                    if w > 0.0:
                        val += w * phi_c[ci + di, cj + dj, ck + dk]
        phi_f[i, j, k] += val


@ti.kernel
def _accumulate_moments(rho: ti.template(), n: ti.i32, dx: ti.f32, out: ti.template()):
    """One-pass device reduction of the 10 raw mass moments into ``out`` (f64, length 10).

    Layout: ``out = [M, Σmx, Σmy, Σmz, Σmx², Σmy², Σmz², Σmxy, Σmxz, Σmyz]`` where
    ``m = rho·dx³`` is the cell mass and ``(x,y,z) = (i,j,k)·dx`` the node position. The
    host turns these into the center of mass and the second moments about it (parallel-axis).
    Reductions target ``out[c]`` at compile-time-constant indices so Taichi can thread-local
    them; f64 keeps the 256³-cell sum accurate (Σ over ~1.7e7 cells), matching the NumPy
    reference. f64 is CUDA/CPU only — see the module docstring (Metal = Stage 6).
    """
    for c in range(10):
        out[c] = 0.0
    vol = ti.cast(dx, ti.f64) ** 3
    for i, j, k in ti.ndrange(n, n, n):
        m = ti.cast(rho[i, j, k], ti.f64) * vol
        x = ti.cast(i, ti.f64) * dx
        y = ti.cast(j, ti.f64) * dx
        z = ti.cast(k, ti.f64) * dx
        out[0] += m
        out[1] += m * x
        out[2] += m * y
        out[3] += m * z
        out[4] += m * x * x
        out[5] += m * y * y
        out[6] += m * z * z
        out[7] += m * x * y
        out[8] += m * x * z
        out[9] += m * y * z


@ti.kernel
def _zero(field: ti.template()):
    for idx in ti.grouped(field):
        field[idx] = 0.0


@ti.kernel
def _scale_into(src: ti.template(), dst: ti.template(), factor: ti.f32):
    for idx in ti.grouped(src):
        dst[idx] = factor * src[idx]


@ti.kernel
def _copy(src: ti.template(), dst: ti.template()):
    for idx in ti.grouped(src):
        dst[idx] = src[idx]


@ti.kernel
def _set_dirichlet_faces(
    phi: ti.template(),
    n: ti.i32,
    dx: ti.f32,
    g_const: ti.f32,
    mass: ti.f32,
    cmx: ti.f32,
    cmy: ti.f32,
    cmz: ti.f32,
    sxx: ti.f32,
    syy: ti.f32,
    szz: ti.f32,
    sxy: ti.f32,
    sxz: ti.f32,
    syz: ti.f32,
):
    """Stamp the six box faces with the multipole potential (monopole + quadrupole)."""
    tr_s = sxx + syy + szz
    for i, j, k in ti.ndrange(n, n, n):
        if i == 0 or i == n - 1 or j == 0 or j == n - 1 or k == 0 or k == n - 1:
            rx = i * dx - cmx
            ry = j * dx - cmy
            rz = k * dx - cmz
            d2 = rx * rx + ry * ry + rz * rz + 1e-12
            d = ti.sqrt(d2)
            nx, ny, nz = rx / d, ry / d, rz / d
            n_s_n = (
                nx * nx * sxx
                + ny * ny * syy
                + nz * nz * szz
                + 2.0 * (nx * ny * sxy + nx * nz * sxz + ny * nz * syz)
            )
            quad = 3.0 * n_s_n - tr_s
            phi[i, j, k] = -g_const * (mass / d + quad / (2.0 * d2 * d))


class MultigridPoissonSolver(PoissonSolver):
    """Open-BC Poisson solve by red-black-GS multigrid V-cycles (the production default)."""

    def __init__(
        self,
        grid_size: int,
        dx: float = 1.0,
        grav_constant: float | None = None,
        n_cycles: int = 20,
        pre_sweeps: int = 2,
        post_sweeps: int = 2,
        coarse_sweeps: int = 40,
    ):
        self.grid_size = grid_size
        self.dx = float(dx)
        self.G = units.G if grav_constant is None else grav_constant
        self.n_cycles = n_cycles
        self.nu1 = pre_sweeps
        self.nu2 = post_sweeps
        self.coarse_sweeps = coarse_sweeps

        # Build the level sizes (finest first), halving while a coarser grid still has
        # an interior (≥ 3 nodes/axis) and the parent is even.
        sizes = [grid_size]
        while sizes[-1] % 2 == 0 and sizes[-1] // 2 >= 3:
            sizes.append(sizes[-1] // 2)
        self.sizes = sizes
        self.dxs = [self.dx * (2**level) for level in range(len(sizes))]

        # One phi/rhs/res field per level.
        self.phi = [ti.field(ti.f32, shape=(s, s, s)) for s in sizes]
        self.rhs = [ti.field(ti.f32, shape=(s, s, s)) for s in sizes]
        self.res = [ti.field(ti.f32, shape=(s, s, s)) for s in sizes]

        # Device-resident accumulator for the 10 raw boundary moments (RV6a). f64 so the
        # 256³ sum stays accurate; CUDA/CPU only (Metal = Stage 6).
        self._moment_buf = ti.field(ti.f64, shape=10)

    # --- boundary moments -------------------------------------------------------

    def _moments_device(self, rho):
        """Device reduction of the boundary moments from the live ``rho`` field (RV6a).

        Returns the same ``(M, (cx,cy,cz), (sxx,syy,szz,sxy,sxz,syz))`` tuple as
        :meth:`_moments`, but reduces on the device and copies only 10 scalars to the host —
        no per-step full-grid ``to_numpy``.
        """
        _accumulate_moments(rho, self.grid_size, self.dx, self._moment_buf)
        s = self._moment_buf.to_numpy()  # 10 f64 scalars
        total = float(s[0])
        if total <= 0.0:
            return 0.0, (0.0, 0.0, 0.0), (0.0,) * 6
        cx, cy, cz = s[1] / total, s[2] / total, s[3] / total
        # Second moments about the CM via parallel-axis: S_ab = Σm·a·b − M·c_a·c_b.
        sxx = float(s[4] - total * cx * cx)
        syy = float(s[5] - total * cy * cy)
        szz = float(s[6] - total * cz * cz)
        sxy = float(s[7] - total * cx * cy)
        sxz = float(s[8] - total * cx * cz)
        syz = float(s[9] - total * cy * cz)
        return total, (float(cx), float(cy), float(cz)), (sxx, syy, szz, sxy, sxz, syz)

    def _moments(self, rho_np: np.ndarray):
        """Reference (host NumPy) moments — kept to unit-test :meth:`_moments_device`.

        Returns (M, CM, S_ab) of the mass about its center of mass (internal units)."""
        n, dx = self.grid_size, self.dx
        mass = rho_np * dx**3
        total = float(mass.sum())
        axis = np.arange(n, dtype=np.float64) * dx
        if total <= 0.0:
            zero6 = (0.0,) * 6
            return 0.0, (0.0, 0.0, 0.0), zero6
        m_i = mass.sum(axis=(1, 2))
        m_j = mass.sum(axis=(0, 2))
        m_k = mass.sum(axis=(0, 1))
        cx = float((m_i * axis).sum() / total)
        cy = float((m_j * axis).sum() / total)
        cz = float((m_k * axis).sum() / total)
        # Second moments about the CM: S_ab = Σ m (a−c_a)(b−c_b).
        sxx = float((m_i * axis**2).sum() - total * cx**2)
        syy = float((m_j * axis**2).sum() - total * cy**2)
        szz = float((m_k * axis**2).sum() - total * cz**2)
        m_ij = mass.sum(axis=2)
        m_ik = mass.sum(axis=1)
        m_jk = mass.sum(axis=0)
        sxy = float(axis @ m_ij @ axis - total * cx * cy)
        sxz = float(axis @ m_ik @ axis - total * cx * cz)
        syz = float(axis @ m_jk @ axis - total * cy * cz)
        return total, (cx, cy, cz), (sxx, syy, szz, sxy, sxz, syz)

    # --- V-cycle ----------------------------------------------------------------

    def _smooth(self, level: int, sweeps: int) -> None:
        n = self.sizes[level]
        h2 = self.dxs[level] ** 2
        for _ in range(sweeps):
            _gauss_seidel_rb(self.phi[level], self.rhs[level], n, h2, 0)
            _gauss_seidel_rb(self.phi[level], self.rhs[level], n, h2, 1)

    def _vcycle(self, level: int) -> None:
        if level == len(self.sizes) - 1:
            self._smooth(level, self.coarse_sweeps)
            return
        self._smooth(level, self.nu1)
        n = self.sizes[level]
        _residual(self.phi[level], self.rhs[level], self.res[level], n, 1.0 / self.dxs[level] ** 2)
        _restrict_full_weight(self.res[level], self.rhs[level + 1], self.sizes[level + 1])
        _zero(self.phi[level + 1])
        self._vcycle(level + 1)
        _prolong_add(self.phi[level + 1], self.phi[level], n)
        self._smooth(level, self.nu2)

    # --- public API -------------------------------------------------------------

    def solve(self, rho, phi, warm_start: bool = False) -> None:
        # Boundary moments reduced on the device (RV6a) — only 10 scalars cross to the host,
        # not the full 256³ grid as the Stage-3 reference did.
        total, (cx, cy, cz), (sxx, syy, szz, sxy, sxz, syz) = self._moments_device(rho)

        # Finest RHS = 4πG ρ. Initial interior guess: reuse `phi` from the previous step
        # (warm start, fast in a time loop) or cold-start at 0. Faces are always (re)set
        # from the current multipole moments, since the mass distribution has moved.
        _scale_into(rho, self.rhs[0], _FOUR_PI * self.G)
        if warm_start:
            _copy(phi, self.phi[0])
        else:
            _zero(self.phi[0])
        _set_dirichlet_faces(
            self.phi[0], self.grid_size, self.dx, self.G,
            total, cx, cy, cz, sxx, syy, szz, sxy, sxz, syz,
        )

        for _ in range(self.n_cycles):
            self._vcycle(0)

        _copy(self.phi[0], phi)

    def residual_norm(self) -> float:
        """L2 norm of the finest-grid residual (diagnostic; valid after a solve)."""
        n = self.grid_size
        _residual(self.phi[0], self.rhs[0], self.res[0], n, 1.0 / self.dx**2)
        r = self.res[0].to_numpy()
        return float(np.sqrt((r**2).sum()))


__all__ = ["MultigridPoissonSolver"]
