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
The f64 accumulation runs on CPU/CUDA; Metal has no hardware fp64, so there it uses a
**Kahan-compensated fp32** reduction instead (``_accumulate_moments_kahan``, selected via
``backend.supports_fp64`` — Stage 6 / RV15, D22). The pure-NumPy ``_moments`` is retained
as the reference both device paths are unit-tested against.
"""

import numpy as np
import taichi as ti

from galaxy_collision import units
from galaxy_collision.backend import supports_fp64
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
def _accumulate_moments(
    rho: ti.template(), n: ti.i32, dx: ti.f32,
    m0: ti.template(), m1: ti.template(), m2: ti.template(), m3: ti.template(),
    m4: ti.template(), m5: ti.template(), m6: ti.template(), m7: ti.template(),
    m8: ti.template(), m9: ti.template(),
):
    """One-pass device reduction of the 10 raw mass moments into ten **0-D** f64 fields.

    Moments: ``[M, Σmx, Σmy, Σmz, Σmx², Σmy², Σmz², Σmxy, Σmxz, Σmyz]`` where ``m = rho·dx³``
    is the cell mass and ``(x,y,z) = (i,j,k)·dx`` the node position. The host turns these into
    the center of mass and the second moments about it (parallel-axis).

    **Why ten separate 0-D fields, not one length-10 field (Stage 5/5B perf):** Taichi's
    thread-local-storage reduction only fires for reductions into a *0-D scalar* field
    (``m[None] += …``). Reducing into an *indexed* element (``out[c] += …``) falls back to
    global atomics, and 256³ cells contending on a handful of addresses cost ~175 ms here
    (vs ~4 ms with TLS) — measured. f64 keeps the ~1.7e7-cell sum accurate (CPU/CUDA; the
    Metal fp32 path is :func:`_accumulate_moments_kahan`, see the module docstring).
    """
    m0[None] = 0.0
    m1[None] = 0.0
    m2[None] = 0.0
    m3[None] = 0.0
    m4[None] = 0.0
    m5[None] = 0.0
    m6[None] = 0.0
    m7[None] = 0.0
    m8[None] = 0.0
    m9[None] = 0.0
    vol = ti.cast(dx, ti.f64) ** 3
    for i, j, k in ti.ndrange(n, n, n):
        m = ti.cast(rho[i, j, k], ti.f64) * vol
        x = ti.cast(i, ti.f64) * dx
        y = ti.cast(j, ti.f64) * dx
        z = ti.cast(k, ti.f64) * dx
        m0[None] += m
        m1[None] += m * x
        m2[None] += m * y
        m3[None] += m * z
        m4[None] += m * x * x
        m5[None] += m * y * y
        m6[None] += m * z * z
        m7[None] += m * x * y
        m8[None] += m * x * z
        m9[None] += m * y * z


@ti.kernel
def _sum_squares(field: ti.template(), out: ti.template()):
    """Σ field² as a single f64 reduction into the **0-D** field ``out`` (device-side L2 norm).

    ``out[None]`` (not an indexed element) so Taichi's thread-local reduction fires — see
    ``_accumulate_moments`` for why this matters (~0.4 ms vs ~28 ms at 256³)."""
    out[None] = 0.0
    for idx in ti.grouped(field):
        out[None] += ti.cast(field[idx], ti.f64) ** 2


@ti.kernel
def _accumulate_moments_kahan(rho: ti.template(), n: ti.i32, dx: ti.f32, part: ti.template()):
    """Kahan-fp32 reduction of the 10 raw moments, for fp64-less backends (Metal, RV15/D22).

    Same moments as :func:`_accumulate_moments`, but f64 is illegal on Metal. Parallelizes
    over the ``(i,j)`` column (``n²`` lanes, division-free); each lane Kahan-sums its ``n``
    cells along ``k`` into ten fp32 partials at ``part[c,i,j]``, finished by the host in fp64
    (a tiny fixed copy, never a full-grid round-trip). For the *moments* specifically the
    column structure does most of the work — x,y are constant down a column, so eight of the
    ten components reduce to const·(column mass), and the fp64 cross-column combine is exact;
    Kahan mainly protects ``Σm`` and the z-moments (which vary along k). Measured ~5e-8 rel.
    at 256³ (the column-hybrid alone is already ~1e-8; the ~1e-3 baseline is a single global
    fp32 accumulator). See ``backend`` / ``docs/stage6_plan.md`` 6A."""
    vol = dx * dx * dx
    for i, j in ti.ndrange(n, n):
        s = ti.Vector.zero(ti.f32, 10)
        c = ti.Vector.zero(ti.f32, 10)
        x = ti.cast(i, ti.f32) * dx
        y = ti.cast(j, ti.f32) * dx
        for k in range(n):
            z = ti.cast(k, ti.f32) * dx
            m = rho[i, j, k] * vol
            v = ti.Vector([m, m * x, m * y, m * z, m * x * x, m * y * y,
                           m * z * z, m * x * y, m * x * z, m * y * z])
            yk = v - c       # Kahan: subtract carried compensation
            t = s + yk
            c = (t - s) - yk  # new compensation = (sum error this step)
            s = t
        for cc in ti.static(range(10)):
            part[cc, i, j] = s[cc]


@ti.kernel
def _sum_squares_kahan(field: ti.template(), n: ti.i32, part: ti.template()):
    """Kahan-fp32 Σ field² over a cubic ``(n,n,n)`` field, for fp64-less backends (Metal).

    The fp32 counterpart of :func:`_sum_squares`. Parallelizes over the ``(i,j)`` column;
    each lane Kahan-sums the ``n`` squares along ``k`` into ``part[i,j]``; the host finishes
    the ``n²→1`` sum in fp64. Both call sites operate on the finest grid, so cubic is safe."""
    for i, j in ti.ndrange(n, n):
        s = 0.0
        c = 0.0
        for k in range(n):
            v = field[i, j, k] * field[i, j, k]
            yk = v - c
            t = s + yk
            c = (t - s) - yk
            s = t
        part[i, j] = s


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
        tol: float = 1e-3,
        min_cycles: int = 1,
    ):
        self.grid_size = grid_size
        self.dx = float(dx)
        self.G = units.G if grav_constant is None else grav_constant
        self.n_cycles = n_cycles  # max V-cycles per solve (cap)
        self.nu1 = pre_sweeps
        self.nu2 = post_sweeps
        self.coarse_sweeps = coarse_sweeps
        # Adaptive cycling (Stage 5/5B, RV5): stop V-cycles once the finest-grid residual norm
        # falls below `tol` × ‖rhs‖, but never fewer than `min_cycles`. A *warm-started* step
        # begins near the solution and reaches the discretization floor (~6e-5) in one cycle, so
        # it stops at `min_cycles`; a cold step (residual ~2e-2 even at 20 cycles) can't beat
        # `tol` and runs the full `n_cycles` cap — so cold-start behavior (and RV10) is unchanged.
        # `last_cycles` records how many cycles the most recent solve actually ran.
        self.tol = float(tol)
        self.min_cycles = max(1, int(min_cycles))
        self.last_cycles = 0

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

        # Boundary-moment + L2-norm reductions pick their path from the backend's fp64
        # support (Stage 6 / RV15, D22): f64 thread-local reductions on CPU/CUDA (unchanged),
        # a Kahan-compensated fp32 reduction on Metal (no hardware fp64). CONTRACT: the f64
        # kernels (_accumulate_moments, _sum_squares) must never be *called* on Metal — Taichi
        # JITs a kernel lazily on first call, so an uncalled f64 kernel costs nothing, but
        # calling one on Metal fails to compile ("Type f64 not supported"). _fp64 gates that.
        self._fp64 = supports_fp64()
        if self._fp64:
            # Ten 0-D f64 accumulators for the raw boundary moments (RV6a). 0-D (not a
            # length-10 field) so each reduction gets Taichi's thread-local optimization —
            # see _accumulate_moments. Plus a 0-D f64 scratch for the L2-norm reductions.
            self._mom = tuple(ti.field(ti.f64, shape=()) for _ in range(10))
            self._sq_buf = ti.field(ti.f64, shape=())
        else:
            # Metal: (i,j)-column fp32 partials, finished in host fp64 (_accumulate_moments_kahan
            # / _sum_squares_kahan). part shape is grid-bound, not N-bound — a tiny fixed copy.
            self._mom_part = ti.field(ti.f32, shape=(10, grid_size, grid_size))
            self._sq_part = ti.field(ti.f32, shape=(grid_size, grid_size))

    # --- boundary moments -------------------------------------------------------

    def _moments_device(self, rho):
        """Device reduction of the boundary moments from the live ``rho`` field (RV6a).

        Returns the same ``(M, (cx,cy,cz), (sxx,syy,szz,sxy,sxz,syz))`` tuple as
        :meth:`_moments`, but reduces on the device and copies only 10 scalars to the host —
        no per-step full-grid ``to_numpy``.
        """
        if self._fp64:
            _accumulate_moments(rho, self.grid_size, self.dx, *self._mom)
            s = [float(m[None]) for m in self._mom]  # 10 f64 scalars to host
        else:
            # Metal: fp32 column partials, summed across columns in host fp64 (RV15/D22).
            _accumulate_moments_kahan(rho, self.grid_size, self.dx, self._mom_part)
            s = self._mom_part.to_numpy().reshape(10, -1).astype(np.float64).sum(axis=1).tolist()
        total = s[0]
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

    def _l2_device(self, field) -> float:
        """L2 norm of ``field`` via a device sum-of-squares reduction.

        f64 thread-local reduction (one scalar to host) on CPU/CUDA; a Kahan-fp32 column
        reduction finished in host fp64 on Metal (RV15/D22)."""
        if self._fp64:
            _sum_squares(field, self._sq_buf)
            return float(self._sq_buf[None]) ** 0.5
        _sum_squares_kahan(field, self.grid_size, self._sq_part)
        return float(self._sq_part.to_numpy().astype(np.float64).sum()) ** 0.5

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

        # Adaptive V-cycling: stop once ‖residual‖ < tol·‖rhs‖ (but ≥ min_cycles), capped at
        # n_cycles. ‖rhs‖ is computed once; the per-cycle residual reduction is device-side.
        n = self.grid_size
        inv_h2 = 1.0 / self.dx**2
        rhs_norm = self._l2_device(self.rhs[0])
        self.last_cycles = self.n_cycles
        for c in range(self.n_cycles):
            self._vcycle(0)
            if c + 1 >= self.min_cycles and rhs_norm > 0.0:
                _residual(self.phi[0], self.rhs[0], self.res[0], n, inv_h2)
                if self._l2_device(self.res[0]) < self.tol * rhs_norm:
                    self.last_cycles = c + 1
                    break

        _copy(self.phi[0], phi)

    def residual_norm(self) -> float:
        """L2 norm of the finest-grid residual (diagnostic; valid after a solve)."""
        n = self.grid_size
        _residual(self.phi[0], self.rhs[0], self.res[0], n, 1.0 / self.dx**2)
        r = self.res[0].to_numpy()
        return float(np.sqrt((r**2).sum()))


__all__ = ["MultigridPoissonSolver"]
