# Stage 6 — Apple / Metal (DRAFT plan)

> **Status:** draft for owner review (2026-06-25), not yet folded into `AGENT.md` §7.
> Scoped on the primary dev box: **Apple M5 Pro, 64 GB unified memory** (D7/D17 target tier),
> Python 3.12, Taichi `arch=metal` confirmed live. Stage 5 (CUDA) is ✅; the Stage 3–4 suite is
> green on CPU **and** CUDA (123 tests).

## 0. Goal & exit gate (from AGENT.md §7)

> **Outcome:** first-class Apple GPU as a production backend; the fp32 compute policy (D6) realized
> end-to-end on Metal.
> **Exit gate:** **cross-backend parity (test 6)** — CPU vs CUDA vs Metal agree within an fp32
> tolerance on a fixed IC; a **perf benchmark vs CUDA**; **re-pass the Stage 3–4 test suite** on
> Metal (and keep it green on CPU + CUDA).

Three owner decisions scoped this draft (2026-06-25):

- **D22 (proposed) — Metal precision path:** the two device reductions that accumulate in `ti.f64`
  (`multigrid._moments_device`, `diagnostics_device`) get a **Kahan-fp32 (compensated-sum)** path
  on-device. Both stay device-resident — no CPU offload — so Metal does not reintroduce the
  per-step host hop that 5A removed. Matches D6's first option; **closes RV15**.
- **D23 (proposed) — Deposition on Metal:** open Stage 6 with a **Phase-0 CIC-deposit spike** on
  Metal (the R2 risk: large atomic-scatter throughput). If atomics underperform badly, the
  documented fallback is **Vulkan-compute + VkFFT** (one GPU source across NVIDIA + Apple);
  Vulkan is a contingency, not in scope unless the spike forces it.
- **D24 (proposed) — Headline benchmark:** **match the CUDA matrix up to 100M** on the M5 Pro for a
  direct CUDA-vs-Metal comparison. The 100M sim state is ~5.5 GB against **64 GB unified memory** —
  ample headroom (see §4), so memory is not the Metal ceiling that it is on the 8 GB 3070.

Stage 6 closes deferred review item **RV15** and **completes test 6** (3-way cross-backend
determinism — CPU↔CUDA landed in 5A; Metal is the third leg here). RV17 (gather tuning) may resurface
in profiling but stays optional.

## 1. What's already true (baseline on this box)

- **Backend plumbing exists.** `sim._ARCH_NAMES` includes `metal`; `_hardware_description` reports
  `Apple M5 Pro — Apple GPU (Metal)`; `init_backend` warns loudly on a silent CPU fallback. A
  `metal` run will *announce* itself correctly today.
- **One portable kernel source (D2).** Deposit/gather, the multigrid stencil, integrator, and the
  diagnostics reductions are all `@ti.kernel` — they *compile* for `arch=metal` without a rewrite.
- **`GALAXY_TEST_ARCH` override** (added 5A) already lets the full suite re-run on an arbitrary
  backend — the mechanism the Metal re-pass rides on.
- **The one hard blocker is fp64.** `multigrid._moments_device`/`_sum_squares` and every accumulator
  in `diagnostics_device` use `ti.f64`. Metal has **no hardware fp64** (D6) — these kernels will
  fail to compile or silently degrade on `arch=metal` until they get an fp32 path (RV15).

## 2. The gaps Stage 6 must close (grounded in the code)

| Gap | Where | Why it blocks Metal |
|---|---|---|
| f64 boundary-moment reduction | `solver/multigrid.py:_accumulate_moments` / `_sum_squares` | `ti.cast(..., ti.f64)` per cell; runs **every solve**. No fp64 on Metal (RV15). |
| f64 diagnostics accumulators | `diagnostics_device.py` (`_PSCALARS` buffer, `_reduce_*`, half-mass hist) | KE/PE/momentum/L all accumulate in `ti.f64`; cadence-gated but still Metal-illegal. |
| CIC atomic scatter unprofiled on Metal | `deposit.py:_deposit_cic` | R2: Metal global-atomic throughput is the known unknown; 5B's winner was tuned on CUDA. |
| Parity test stops at CUDA | `tests/test_determinism.py` | Asserts CPU↔CUDA only; test 6 wants the 3-way agreement incl. Metal. |
| No Metal perf numbers | `docs/performance.md` | The exit gate wants a Metal benchmark vs the CUDA matrix. |

## 3. Checkpoint structure (mirrors Stage 5's 5A/5B/5C)

### 6A — Metal bring-up & correctness re-pass ✅ **Done (2026-06-25)**
*Foundation. Make the GPU path **correct** on Metal first — Kahan-fp32 reductions + a clean suite
re-pass — before any perf claim.*

> **Delivered.** (1) Confirmed f64 is a **hard compile error** on Taichi/Metal
> (`Type f64 not supported`) — RV15 was real, not theoretical. (2) `backend.supports_fp64()`
> selects the reduction path; the f64 thread-local path is **unchanged on CPU/CUDA**, and Metal
> uses a **Kahan-compensated fp32** path — `multigrid._accumulate_moments_kahan` /
> `_sum_squares_kahan` (boundary moments + adaptive-cycling L2 norm) and
> `diagnostics_device._reduce_{particle_scalars,grid_pe}_kahan` + an fp32 radial histogram —
> partials finished in host fp64 (a tiny fixed copy, no full-grid round-trip). (3) **Measured
> accuracy** vs the fp64 numpy reference: moments ~5e-8, KE ~3e-10, grid-PE ~4e-9 rel. Two
> distinct things buy this, and the credit splits: the **column-hybrid structure** (short fp32
> per-column sums + an exact fp64 cross-column combine) does most of the work for the *moments*
> — their position factors are constant per column, so most components reduce to const·(column
> mass), and column-hybrid *alone* already lands ~1e-8 there — while **Kahan** is the dominant
> contributor for the **grid-PE and L2-norm** reductions, where ρ·Φ / the residual vary along
> the column and don't factor. The ~1e-3 figure both improve on is a single global fp32
> accumulator (no hybrid, no Kahan). (4) Whole PM pipeline runs end-to-end on `arch=metal`; no
> other Metal-specific kernel issue surfaced. (5) `test_device_diagnostics_match_host_reference`
> made arch-aware (tight on CPU/CUDA; fp32-scaled on Metal — near-zero momentum compared against
> the Σm|v| summation scale, not rtol). (6) New `test_cpu_metal_trajectories_agree` adds the
> Metal leg of **test 6**: CPU↔Metal tracer paths agree to ~2e-6 kpc over 15 steps (cap 1e-3).
> Note test 6 is **two pairwise legs sharing the CPU anchor** (CPU↔CUDA on the 3070, CPU↔Metal
> on the Mac), not a simultaneous 3-way run — no box has both GPUs, so CUDA↔Metal holds only
> transitively through CPU. **Suite green on CPU and Metal (123); ruff clean.** RV15 closed.

1. **Kahan-fp32 reduction path (D22, RV15).** Add a compensated-summation fp32 reduction for both
   `multigrid._moments_device` and the `diagnostics_device` accumulators, selected by backend
   (fp64 on CPU/CUDA, Kahan-fp32 on Metal). Keep the existing host fp64 `diagnostics` + NumPy
   `_moments` as the numeric reference the fp32 path is unit-tested against; assert the fp32 path
   reproduces the fp64 reduction within a derived tolerance on a fixed ρ / fixed particle set.
2. **Metal pipeline bring-up.** Run the full PM pipeline on `arch=metal`: Plummer (small N) and a
   `two_galaxy_4v` smoke at 256³. Fix whatever Metal-specific kernel issues surface (atomic types,
   template specializations, unsupported ops).
3. **Parity test (test 6, RV-complete).** Extend `tests/test_determinism.py` with a CPU↔Metal leg
   so the same IC + seed agrees within an fp32 trajectory tolerance over a short run. This makes
   test 6 **two pairwise legs sharing the CPU anchor** (CPU↔CUDA + CPU↔Metal), not a simultaneous
   3-way run — no single box has both GPUs, so CUDA↔Metal holds only transitively via CPU. (Each
   leg skips gracefully when its arch is absent.)
4. **Suite re-pass on Metal.** `GALAXY_TEST_ARCH=metal pytest` green; ruff clean. Document any test
   that must widen its tolerance for fp32-on-Metal (and why).

**Exit 6A:** full suite green on **CPU + CUDA + Metal**; no fp64 in any kernel that runs on Metal;
RV15 closed; the 3-way parity test added and passing.

### 6B — Performance: deposition spike + profile + benchmark ✅ **Done (2026-06-25)**
*Earn "first-class GPU," not just "it runs."*

> **Delivered (`docs/performance.md` Stage-6 section).** Metal benchmark matrix at 1M/10M/30M/100M
> on the M5 Pro: **19.1 / 5.6 / 2.1 / 0.65 steps/s**. **100M fits** (21.5 GiB peak RSS / 64 GB —
> unified memory removes the 8 GB-card ceiling; `bench` gained a backend-aware memory probe:
> nvidia-smi on CUDA, peak process RSS on Metal/CPU). **R2 confirmed and characterized (D23):**
> deposit is the throughput ceiling at scale (69% @10M → 86% @100M, ~50× CUDA's deposit). The spike
> *ruled out* the two cheap fixes — it is **not** atomic contention (racy non-atomic `+=` is
> identical, 1.00×, even with 23.7k particles/hot-cell) and **not** particle order (cell-sort
> doesn't help). Root cause: **global-memory write-scatter, geometry-bound** — gather-reads are ~19×
> cheaper than deposit-writes for the same stencil, and a compact gaussian deposits 16× faster than
> the thin-disk galaxy. **Decision: accept + document** (owner, 2026-06-25); logged **RV20** for a
> possible future block-local-privatization / grid-relayout attempt; Vulkan fallback not triggered.
> The solve (grid-bound) stays competitive (~40–60 ms vs CUDA ~28 ms), so at 1M Metal is 0.57× CUDA.

1. **Phase-0 CIC-deposit spike (D23, R2).** Micro-benchmark the deposit on Metal across the 5B
   variants (resident `atomic_add`; per-block privatization; sort-by-cell) at 1M/10M/30M. If Metal
   atomic-scatter underperforms badly, characterize it and decide: accept, switch variant, or escalate
   to the Vulkan fallback (out of scope unless forced — flag loudly if so).
2. **Per-stage profile on Metal.** Reuse the 5B profiler (`galaxy-bench`) to produce the deposit /
   solve / grad / gather / integrate breakdown on Metal. Watch RV17 (gather dominated at 100M on
   CUDA) — note whether Metal shifts the bottleneck.
3. **Benchmark matrix.** steps/sec at N ∈ {1M, 10M, 30M, 100M} on Metal, with the per-stage
   breakdown and peak unified-memory footprint.

**Exit 6B:** Metal benchmark matrix produced; deposit characterized (winner chosen or Vulkan
escalation documented); profile captured.

### 6C — CUDA-vs-Metal write-up + headline run + doc sync ⬜

1. **Headline run (D24).** A `two_galaxy_4v` 100M / 256³ collision on the M5 Pro; capture
   steps/s + peak unified memory.
2. **CUDA-vs-Metal comparison.** Add a Metal section to `docs/performance.md` with the side-by-side
   matrix (3070 vs M5 Pro) and the per-stage profile; honest read on where each wins.
3. **Doc sync.** `AGENT.md` §7 (Stage 6 → ✅, mark RV15 done), §1 status, README status/roadmap,
   §11 review log; per the maintenance contract.

**Exit 6C (= Stage-6 exit gate):** 3-way parity green; Metal benchmark vs CUDA documented; Stage 3–4
suite green on Metal; docs synced.

## 4. Memory budget on the M5 Pro (64 GB unified)

Same per-particle footprint as Stage 5: **44 B/particle** (32 B SoA + 12 B accel) + ~0.57 GB
N-independent grid. Unified memory means the CUDA 8 GB ceiling does **not** apply.

| N | Particles (44 B) | + grid | Peak (approx.) | Headroom (64 GB) |
|---|---|---|---|---|
| 1M | 0.04 GB | ~0.6 GB | **~0.6 GB** | trivial |
| 10M | 0.44 GB | ~1.0 GB | **~1.0 GB** | trivial |
| 30M | 1.32 GB | ~1.9 GB | **~1.9 GB** | trivial |
| 100M | 4.40 GB | ~5.0 GB | **~5.0 GB** | ~59 GB free — **not the ceiling** |

Memory is a non-issue at 100M here; the open question is *throughput* (the deposit spike), not fit.
Unified memory may also change the cost model — no discrete host↔device PCIe copy — which is worth a
note in the write-up.

## 5. Risks & mitigations (Stage-6-specific)

- **Metal atomic-scatter underwhelms (R2)** → the Phase-0 spike (6B) is the mitigation;
  Vulkan-compute + VkFFT is the documented escalation if it's a hard wall.
- **fp32 reductions lose accuracy vs fp64 reference** → Kahan/compensated summation (D22), tolerance
  derived against the fp64 path and asserted in 6A; diagnostics are a *monitor*, not a hard gate at
  scale (consistent with §6 / Stage-4 scope).
- **Metal-specific kernel feature gaps** (unsupported ops, atomic types) → surfaced in 6A bring-up
  on small N before any scale run; one portable source (D2) keeps fixes localized.
- **fp32 trajectory divergence across backends** → test 6 asserts agreement *within an fp32
  tolerance*, not bitwise; the tolerance is derived and documented (as in the CPU↔CUDA test).

## 6. Out of scope (stays later)

Realtime GGUI + batch→movie + the 100M density-projection panel (Stage 7, incl. RV19); warm-disk ICs
/ research campaigns (Stage 8, D16); RV13 tracer-cadence decoupling and RV17 gather tuning (optional);
the Vulkan/VkFFT backend (only if the R2 spike forces it).
