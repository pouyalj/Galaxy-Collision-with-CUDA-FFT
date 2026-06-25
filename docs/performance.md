# Performance — Stage 5 / 5B (CUDA, RTX 3070)

Benchmark + per-stage profile of the PM force chain on the Stage-5 box (**NVIDIA RTX 3070,
8 GB**, driver 595.71.05 / CUDA 13.2, Taichi 1.7.4, 256³ grid, multigrid solver). Numbers from
`galaxy-bench` (`src/galaxy_collision/bench.py`): a clean throughput pass (one final `ti.sync`)
plus a per-stage pass (each stage `ti.sync`-bracketed). Disks equilibrated, warm-started — i.e.
representative of a real run, not a cold transient.

## Headline: 6.5–12.5× throughput

| N | before 5B (steps/s) | after 5B (steps/s) | speedup | ms/step | peak VRAM |
|---|---|---|---|---|---|
| 1M | 2.66 | **33.4** | 12.5× | 30.0 | 827 MiB |
| 10M | 2.58 | **24.8** | 9.6× | 40.3 | 1339 MiB |
| 30M | 2.43 | **15.8** | 6.5× | 63.3 | 2491 MiB |
| 100M | ~2.0 (est.) | **6.95** | ~3.5× | 143.9 | 6491 MiB |

"before 5B" = the committed 5A solver (fixed 20 V-cycles/step + the contended moment reduction
below). **100M fits the 8 GB card** with ~1.7 GB headroom — the Stage-5 exit-gate scale runs here.

## Per-stage breakdown (after 5B, ms/step)

| stage | 1M | 10M | 30M | 100M |
|---|---|---|---|---|
| deposit | 0.50 | 2.54 | 6.94 | 22.84 |
| **solve** | **27.5** | **27.8** | **27.7** | **27.6** |
| grad | 1.17 | 1.17 | 1.19 | 1.16 |
| gather | 0.63 | 5.88 | 17.70 | 57.89 |
| integrate | 0.40 | 3.18 | 9.90 | 33.95 |
| solve share | 91% | 68% | 44% | 19% |

Two regimes, exactly as the O(grid) vs O(N) cost model predicts:

- **The solve is N-independent (~27.6 ms).** It dominates at small N (91% at 1M) and is a minor
  cost at 100M (19%).
- **Particle stages scale linearly with N.** By 30M the step is particle-bound, and at 100M
  **gather is the single largest stage (40%)** — it reads three acceleration grids × 8 nodes per
  particle, so it moves more memory than deposit. Gather/integrate are the natural *next*
  optimization target (a later stage); deposit is not (see D19 below).

## What made the solve ~13× faster

The committed 5A solve was ~375 ms/step. Profiling decomposed it into two independent problems,
both fixed in 5B:

1. **Contended moment reduction (a 5A regression, ~175 ms → ~3.5 ms).** 5A moved the open-BC
   boundary multipole moments onto the device (RV6a) — correct for residency, but it reduced ten
   moments into an *indexed* field element (`out[c] += …`), which **defeats Taichi's thread-local
   reduction** and leaves 256³ cells contending on a few global atomics (~175 ms). Reducing into
   ten **0-D scalar** fields (`m[None] += …`) restores TLS: a measured **72×** on a single
   reduction (28 ms → 0.39 ms), and the full moment pass drops to ~3.5 ms. The same fix applies to
   the new residual-norm reduction used by (2). *This regression was invisible until 5B profiling —
   5A's tests checked correctness, not cost.*
2. **Over-iteration on warm steps (20 cycles → ~2).** The time loop warm-starts every solve from
   the previous step's Φ, which sits at the discretization floor; **one V-cycle already reaches
   residual 6.2e-5, identical to 20 cycles.** Adaptive cycling (`MultigridPoissonSolver(tol=…)`)
   stops when ‖residual‖ < `tol`·‖rhs‖ (device-side norm, no host copy): a cold step (residual
   ~2e-2 even at 20 cycles) still runs the full cap — preserving RV10 — while warm steps stop at
   ~1–2. Validated against a fixed-20-cycle reference over an 80-step 4v collision: energy drift,
   virial, half-mass-radius trajectory and tracer path all agree to ≤1e-4 (relative) / sub-cell.

Neither fix alone suffices: without (1), a 2-cycle solve would still cost ~175 ms (moments) and
the adaptive win would be masked; without (2), the TLS fix alone takes the solve only ~375→~210 ms.

## D19 — deposition tuning: closed by measurement

The plan (D19) scoped benchmarking three CIC-scatter strategies. The profile shows **deposit is
≤16% of the step even at 100M and <2% at ≤10M**, so the premise (deposit dominates at scale) is
false — solve, then gather, lead. A lightweight A/B confirms tuning deposit is not worth it:

| deposit @ 10M | ms |
|---|---|
| baseline atomics (as generated) | 2.57 |
| pre-sorted by cell index | 3.40 (**0.75× — slower**) |
| one host argsort (the sort cost) | 642 |

Sorting by cell **slowed** deposit here (the IC is already locally coherent; cell-contiguous order
raises intra-cell atomic contention) and the sort itself costs 250× the deposit. So the baseline
`ti.atomic_add` scatter is kept; per-block privatization and a device radix sort are **not
pursued** — they would optimize a non-bottleneck. If a future workload makes deposit dominant,
revisit; the harness is here to re-measure.

## 100M headline run (5C, D20)

The Stage-5 exit-gate run: a **100M-particle `two_galaxy_4v` collision at 256³**, 800 steps
(400 Myr), on the RTX 3070.

| metric | value |
|---|---|
| wall (run loop) | 207 s → **3.86 steps/s** |
| peak VRAM | **6491 MiB / 8192** (fits, ~1.7 GB headroom) |
| virial (final) | 0.91 |
| grid-energy drift (max) | 0.40 — cold-disk heating, *monitor-only* (not a gate; §6/RV11) |
| half-mass radius | 45.7 → **11.5 kpc** at first pericenter (~175 Myr) → 62.8 kpc (tidal re-expansion) |

The collision signature is physically sensible (infall → pericenter compression → tidal
spreading). The 3.86 steps/s here is below the **6.95** of the quiescent benchmark above: the
benchmark measures warm, quasi-static steps (~2 V-cycles), whereas a live collision drives ρ
hard near pericenter, so the *adaptive* solver correctly spends more cycles when the warm guess
is poorer — exactly the behavior it's designed for. No snapshots were written (a 100M snapshot is
~2.4 GB); the diagnostics history + tracer were recorded in memory. The full 400-Myr collision
completes in ~3.5 minutes of compute. Per D20 this run is a **performance/scale demonstration**
validated by diagnostics, not a new image (the science figures were 4C at 10M) — a 100M
density-projection panel would make a stronger headline and is a natural **Stage-7** viz item.

## FFT oracle on GPU (5C, D21)

The validation oracle (`solver/fft_oracle.py`) now runs its transform on **CuPy/cuFFT** when a
CUDA device + CuPy are present (auto-detected; `use_gpu=False` forces NumPy). The 512³ `rfftn`
that backs the 256³ open-BC solve is **~19× faster on the GPU** (≈0.11 s vs ≈2.05 s), so
multigrid can be validated against the spectrally-exact oracle at the *full production 256³ grid*
(`test_multigrid_matches_gpu_oracle_at_production_grid`) rather than only the 64³ grids CPU CI
affords. The math is identical — selected by the array module — so the NumPy path remains the CI
fallback (covered when CuPy is absent). `rho`/`phi` still bridge through the host (Taichi field ↔
NumPy); only the transform is offloaded, which is fine for a validation-only solver.

Install: `pip install -e ".[cuda]"` (pulls `cupy-cuda12x[ctk]` — the `[ctk]` extra is required, as
CuPy JITs its elementwise kernels and needs the CUDA headers; without it only precompiled cuFFT
works, not the Green's-function build).

## Method notes / caveats

- Per-stage times are individually `ti.sync`-bracketed (serialized), so they sum slightly above
  the overlapped throughput ms/step — read them as *relative* cost.
- "before 5B" 100M is an estimate (the pre-5B solver was never run at 100M); all "after" numbers
  are measured.
- f64 reductions (moments, residual norm) are CUDA/CPU only; the Metal Kahan-fp32 path is Stage 6
  (RV15).
- Reproduce: `galaxy-bench --n <N> --grid 256 --solver multigrid --backend cuda`.
