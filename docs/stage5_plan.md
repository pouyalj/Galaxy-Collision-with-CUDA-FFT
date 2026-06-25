# Stage 5 — CUDA & scale-up (DRAFT plan)

> **Status:** draft for owner review (2026-06-25), not yet folded into `AGENT.md` §7.
> Scoped on the provisioned box: **NVIDIA RTX 3070, 8 GB, driver 595.71.05 / CUDA 13.2**,
> Python 3.12, Taichi 1.7.4 (CUDA backend verified live, `fell_back=False`; 118 tests pass).

## 0. Goal & exit gate (from AGENT.md §7)

> **Outcome:** device-resident state, deposition tuning, 100M+ runs.
> **Exit gate:** a 100M-particle run; benchmark + per-stage profile; **re-pass the Stage 3–4 test
> suite** on CUDA (and keep it green on CPU).

Three owner decisions scoped this draft (2026-06-25):

- **D19 (proposed) — Deposition tuning:** benchmark **all three** strategies (plain `atomic_add`,
  per-block privatization, sort-by-cell coalescing); keep the fastest.
- **D20 (proposed) — Exit-gate N:** **100M** on the 3070 (est. ~5.5 GB peak / 8 GB — feasible
  headless, ~2 GB margin; see §4).
- **D21 (proposed) — FFT oracle:** add a **CuPy cuFFT** GPU path so the oracle validates multigrid
  at larger N on-device (optional `cupy-cuda12x` dep; ~1 GB for the 512³ pad, validation-only).

Stage 5 also closes deferred review items **RV5, RV6, RV10** and partially advances test 6
(cross-backend determinism: CPU↔CUDA now, Metal at Stage 6).

## 1. What's already true (baseline, measured on this box)

- Taichi resolves `arch=cuda` on the RTX 3070, no silent CPU fallback.
- Full PM pipeline runs end-to-end on CUDA: Plummer 50k/64³ (drift 3.6e-4) and `two_galaxy_4v`
  2M/256³ multigrid.
- **Baseline perf:** ~0.44 s/step at 2M / 256³ multigrid (after ~24 s one-time kernel JIT);
  peak VRAM ~955 MiB at 2M.
- 118/118 tests pass (CPU).

## 2. The gaps Stage 5 must close (grounded in the code)

| Gap | Where | Why it matters at scale |
|---|---|---|
| Per-step `rho.to_numpy()` + host-side multipole moments | `solver/multigrid.py:solve` / `_moments` | A full 256³ host↔device copy **every step** (legacy bug #13 territory; RV6a). |
| Diagnostics pull full fields to host via `to_numpy()` | `sim.py:measure` | Host hop each history sample; fp64 reductions (RV6b). Less hot (cadence-gated) but still a transfer. |
| Split grid ownership | `GridState` (ρ,Φ) vs ad-hoc `ax_g/ay_g/az_g` in `run_simulation` vs MG hierarchy in solver | No single device-resident owner; hard to reason about residency/memory (RV6c). |
| Deposit uses plain `ti.atomic_add` only | `deposit.py:_deposit_cic` | Atomic contention is the expected hot-kernel bottleneck at 10–100M (§5.6). |
| Multigrid converges ~0.5–0.7/cycle (n_cycles=20) | `solver/multigrid.py` | Suboptimal vertex-centered coarsening; correct but slow (RV5). |
| No 256³ residual-convergence assertion | `tests/test_solver_*.py` | Only small grids tested for convergence (RV10). |

## 3. Checkpoint structure (mirrors Stage 4's 4A/4B/4C)

### 5A — Device-resident force chain & correctness re-pass ✅ **Done (2026-06-25)**
*Foundation. No performance claims yet — make the GPU path correct and resident, then optimize.*

1. **Moments on device (RV6a).** Replace `_moments` NumPy reduction with a Taichi reduction kernel
   (M, CM, second moments S_ab) reading `rho` in place. Remove `rho.to_numpy()` from the per-step
   solve. Validate the kernel reproduces the NumPy moments within fp tolerance on a fixed ρ.
2. **Grid-ownership consolidation (RV6c).** Give `GridState` the acceleration fields
   (`ax_g/ay_g/az_g`) alongside ρ/Φ; keep the MG hierarchy inside the solver. `run_simulation`
   stops allocating grids ad-hoc. One device-resident owner per field.
3. **Device-side diagnostics (RV6b).** Compute KE / momentum / grid-PE via Taichi reductions
   (Kahan/compensated fp32, fp64 where the backend allows) so the cadence-gated `measure` no longer
   round-trips full fields. Pull only scalars to host. *Keep the existing fp64 CPU path as the
   numeric reference the device path is tested against.*
4. **RV10 — 256³ residual regression test.** Assert the production multigrid (n_cycles=20 @ 256³)
   drives the finest-grid residual norm below a derived bound. Now cheap on GPU.
5. **CPU↔CUDA determinism (partial test 6).** Same IC + seed on CPU and CUDA must agree within an
   fp32 trajectory tolerance over a short run. (Full 3-way parity incl. Metal is Stage 6.)

**Exit 5A:** no `.to_numpy()` in the per-step hot path (solve + integrate); Stage 3–4 suite green on
**both** CPU and CUDA; RV10 + CPU↔CUDA tests added and passing.

> **Delivered.** (1) `multigrid._moments_device` — 10-scalar device reduction, per-step
> `rho.to_numpy()` removed. (2) `GridState` owns ρ/Φ + accel grids (RV6c). (3)
> `diagnostics_device.DeviceDiagnostics` — KE/PE/momentum/ang-mom + half-mass (radial histogram)
> reduced on-device; `measure()` copies no full field on a history-only sample. (4)
> `test_multigrid_256_residual_regression` (RV10). (5) `tests/test_determinism.py` — CPU↔CUDA
> tracer paths are bit-identical over the test horizon. (6) `GALAXY_TEST_ARCH` conftest override →
> **122 tests pass on CPU and on CUDA**; ruff clean. RV6/RV10 closed; RV14/RV15 logged in AGENT.md §11.

### 5B — Performance: profiling + deposition + solver tuning
*Earn the "performance non-negotiable" requirement.*

1. **Per-stage profiler.** Wrap deposit / solve / grad / gather / integrate with Taichi's kernel
   profiler (`ti.init(kernel_profiler=True)` + `ti.profiler.print_kernel_profiler_info()`) plus
   wall-clock per phase. Output a per-stage breakdown table — this *is* part of the exit gate.
2. **Deposition tuning (D19) — benchmark all three:**
   - (a) baseline `ti.atomic_add` (exists);
   - (b) per-block/per-cell privatization to cut global-atomic contention;
   - (c) **sort particles by cell index** (radix) for coalesced, low-contention scatter.
   Profile each at 1M/10M/30M; keep the winner. *Memory note (see §4): variant (c)'s sort scratch
   (~8 B/particle) competes for VRAM at 100M — budget it or fall back to the resident variant for the
   headline run.*
3. **Multigrid convergence tuning (RV5).** Try proper 2^L+1 nesting and/or W-cycles / Galerkin
   coarsening to lift the ~0.5–0.7/cycle rate, reducing `n_cycles` at fixed accuracy. Gate any change
   against the FFT oracle + the RV10 residual test so correctness never regresses.
4. **Benchmark matrix.** steps/sec at N ∈ {1M, 10M, 30M, 100M} on CUDA, with the per-stage breakdown.

**Exit 5B:** benchmark matrix produced; deposit is no longer the dominant stage (or documented why);
solver tuning measured against the oracle.

### 5C — GPU oracle + 100M headline run + write-up

1. **CuPy cuFFT oracle (D21).** Add an optional GPU FFT path to `solver/fft_oracle.py` (CuPy when
   `cupy-cuda12x` is present, NumPy/pocketfft otherwise). Use it to validate multigrid at a larger N
   than the CPU oracle allowed; the 512³ pad (~1 GB) only loads at validation N.
2. **100M headline run (D20).** Run a `two_galaxy` collision at 100M / 256³ on the 3070; capture peak
   VRAM (memory-ceiling validation) and the per-stage profile.
3. **Write-up + doc sync.** A `docs/performance.md` benchmark/profile report; update `AGENT.md` §7
   (Stage 5 → ✅, mark RV5/RV6/RV10 done), `README.md` status, and §11.

**Exit 5C (= Stage-5 exit gate):** 100M run completes with documented benchmark + per-stage profile;
Stage 3–4 suite green on CUDA; docs synced per the maintenance contract.

## 4. Memory budget on the 8 GB RTX 3070

Per-particle device footprint is **44 B**: 32 B SoA (pos×3, vel×3, mass, gid) **+ 12 B** for the
`acc_x/y/z` fields. Grid is N-independent at **~0.57 GB** (ρ, Φ, 3× accel, MG hierarchy). Add
~0.5 GB CUDA context.

| N | Particles (44 B) | + grid + ctx | Peak VRAM (approx.) | Headroom (8 GB) |
|---|---|---|---|---|
| 1M | 0.04 GB | ~1.1 GB | **~1.1 GB** | ample |
| 10M | 0.44 GB | ~1.5 GB | **~1.5 GB** | ample |
| 30M | 1.32 GB | ~2.4 GB | **~2.4 GB** | comfortable |
| 100M | 4.40 GB | ~5.5 GB | **~5.5 GB** | ~2 GB — **the real ceiling** |

**Risks at 100M:** (i) the sort-by-cell deposit variant's scratch (~0.8 GB at 100M) pushes peak to
~6.3 GB — fits, but if it doesn't, the headline run uses the resident non-sort variant; (ii) the
512³ FFT-oracle pad (~1 GB) must **not** be co-resident with a 100M run — oracle validation runs at
small N only.

## 5. Risks & mitigations (Stage-5-specific)

- **Atomic-scatter throughput on Taichi/CUDA underwhelms** → the three-variant benchmark (5B) is the
  mitigation; sort-by-cell is the known fallback winner.
- **100M doesn't fit with chosen deposit variant** → fall back to resident atomics for the headline
  run; document max-N per variant.
- **Multigrid tuning regresses accuracy** → every solver change gated against the FFT oracle + RV10.
- **fp32 over long 100M integration** (R3) → energy/virial drift monitoring; the 100M run is a
  scale/perf demonstration, not a conservation gate (consistent with §6/Stage-4 scope).

## 6. Out of scope (stays later)

Metal / cross-backend 3-way parity (Stage 6); realtime GGUI + movie (Stage 7); warm-disk ICs /
research campaigns (Stage 8); RV13 tracer-cadence decoupling (optional, not Stage 5).
