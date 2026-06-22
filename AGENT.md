# AGENT.md — Galaxy Collision (Milky Way × Andromeda), PM N-body

> **Purpose of this file.** This is the project bible and onboarding doc for any human or AI
> agent working in this repo. It does two things:
> 1. **Documents the existing code as-is** — what it is, how it works, every dependency, and
>    every known bug.
> 2. **Architects the path forward** — a modernized, cross-platform (CPU / NVIDIA CUDA /
>    Apple-Silicon GPU) version that is performant at 10–100M particles and research-grade in
>    physics.
>
> Decisions in §4 were scoped directly with the repo owner (Pouya, a co-author of the original
> 2020 paper). Treat §4 as binding unless the owner changes them.
>
> ### ⚠️ Maintenance contract (binding — applies to every agent, human or AI)
>
> **This file is the living source of truth and MUST be kept in sync with the code at all times.**
> Any agent that touches the project is required to update AGENT.md *in the same change* whenever it:
> - adds/moves/removes files or changes the repo layout (§3.1, §5.8),
> - changes dependencies, build, CI, or config schema (§3.5, §5.8),
> - completes or advances a stage, or alters its exit criteria (§7 — keep the status column current,
>   **and update the Status line in `README.md`** so the user-facing status never lags),
> - makes or revises a scoped decision (§4.2 Decision Record) or resolves an open question (§9),
> - discovers, fixes, or introduces a notable bug/limitation (§3.6).
>
> If a code change and AGENT.md disagree, that is a bug. Treat "update the docs" as part of "done,"
> never a follow-up. When a generated companion (`Galaxy_Collision_Modernization_Plan.docx`) exists,
> AGENT.md governs and should be updated first; regenerate the docx from it as needed.

---

## 1. Quick orientation (TL;DR)

- **What:** A Particle-Mesh (PM) gravitational N-body simulation of a Milky-Way–Andromeda-like
  galaxy collision. Particles are binned onto a 256³ density grid (1 kpc³ cells), Poisson's
  equation is solved on the grid (currently via cuFFT), the potential is differentiated into
  forces, and particles are advanced with a leapfrog-style integrator.
- **Origin:** University of Toronto Scarborough course project **PHYD57 (2020)**, paper
  *"Galaxy Collisions With CUDA and FFT"* — Adams, Lajevardi, LeBlanc, Movaseghi (Apr 15, 2020).
  The paper PDF is the authoritative statement of *intent*; the code is a partial, buggy
  realization of it.
- **State:** The original three near-duplicate CUDA C files (now archived under `legacy/`, see §3.1)
  compile against CUDA 9.2 + cuFFT + DISLIN on a specific 2020 lab machine. Several correctness and
  performance bugs (see §3.6) mean the headline configuration (M = 100,000,000 particles) will not
  realistically run as written. **Modernization status:** Stages 0–3 are **done** — see §7. The
  repo root hosts the `galaxy_collision` Python package with a scaffold, a tested (kpc, Myr, M☉)
  unit system, the SoA particle/grid data model, a reproducible two-galaxy IC generator
  (disk+bulge+central BH, 4v/2v), and — as of Stage 3 — a **correct CPU simulation**: CIC
  deposit/gather, both Poisson solvers (open-BC multigrid + the zero-padded FFT oracle, pulled
  forward from Stage 4), KDK leapfrog (grid-softened PM forces; explicit Plummer softening on the
  direct/diagnostic path), fp64 conservation diagnostics, and
  HDF5/npz snapshot I/O. It is validated by a stable Plummer sphere, a two-body Kepler orbit, and
  multigrid-vs-oracle agreement.
- **Goal (scoped):** Rebuild as **one portable source** (Taichi-style kernels compiling to
  CPU + CUDA + Metal), research-grade physics, **10–100M particles**, **Apple GPU as a
  first-class performance target**, with a **pluggable Poisson solver** (open-boundary multigrid
  as the portable default + a zero-padded FFT solver as an NVIDIA validation oracle). Outputs:
  real-time 3D viewer, batch→movie, data+analysis snapshots, and reproduction of the paper.

---

## 2. What the project is

### 2.1 The science (from the paper)

- Two spiral galaxies resembling the Milky Way and Andromeda, started a short distance apart and
  given bulk velocities toward each other, are evolved through a collision.
- Two runs in the paper: initial approach speed **4v** and **2v**, where v ≈ 125,000 mph is the
  measured approach speed (~402,000 km/h). Run 1 terminated at 40×10⁷ yr; run 2 at 45×10⁷ yr.
- A single "Sun-like" tracer particle (mass 840 M☉, initial galactocentric radius 8.32 kpc) is
  tracked through the collision. (Our uniform-mass reproduction tracks an *existing* disk particle
  nearest 8.32 kpc rather than an 840 M☉ insert — see §3.7.)
- **Known physical shortcomings acknowledged in the paper:** galaxies render as *squares* (crude
  mass deposition / no interpolation), density visibly *dissipates over time* (numerical heating
  + wrong boundary conditions), no gas, and dark matter / central black holes contribute *force
  but no mass* in the grid. These are exactly the issues the modernization should fix.

### 2.2 The numerical method (PM scheme)

Per timestep:

1. **Deposit** particle masses onto the 256³ density grid ρ.
2. **Solve Poisson** ∇²Φ = 4πGρ for the potential Φ. Currently done spectrally:
   FFT(ρ) → multiply by −1/k² → inverse FFT → Φ. (cuFFT, C2C.)
3. **Force** g = −∇Φ via finite differences of Φ.
4. **Integrate** particle positions/velocities (leapfrog, intended).
5. Repeat; periodically render a density image.

Time complexity is O(N) for particle work + O(G log G) for the FFT (G = grid cells), versus
O(N²) for direct summation — this is the whole point of the PM approach.

---

## 3. The current codebase (as-is)

### 3.1 Repository layout

As of **Stage 0**, the repo root hosts the modernized Python package and the original CUDA source
is archived under `legacy/`:

```
GalaxyCollision/
└── Galaxy-Collision-with-CUDA-FFT/     # the actual git repo (origin: github.com/pouyalj/...)
    ├── README.md  ·  AGENT.md  ·  LICENSE (Apache-2.0)
    ├── pyproject.toml                  # package + deps (taichi, numpy, h5py, matplotlib, pyyaml)
    ├── .github/workflows/ci.yml        # CI: ruff + pytest + CPU hello-sim smoke run
    ├── configs/
    │   └── smoke.yaml                  # Stage-0 one-step no-op smoke config
    ├── src/galaxy_collision/
    │   ├── __init__.py
    │   ├── config.py                   # SimConfig schema + YAML load/dump (Stage 0)
    │   ├── sim.py                      # orchestration + hello-sim/galaxy-sim CLIs (Stage 0/3)
    │   ├── units.py                    # (kpc, Myr, M_sun) system + derived G (Stage 1)
    │   ├── data.py                     # SoA particle/grid fields + memory estimator (Stage 1)
    │   ├── ic.py                       # two-galaxy + Plummer initial conditions (Stage 2/3)
    │   ├── deposit.py                  # CIC deposit + grad (−∇Φ) + force gather (Stage 3)
    │   ├── integrator.py              # KDK leapfrog + direct Plummer-softened force (Stage 3)
    │   ├── diagnostics.py              # energy/momentum/Lagrangian-radius, fp64 (Stage 3)
    │   ├── io.py                       # HDF5/npz snapshot I/O (Stage 3)
    │   ├── solver/                     # base + multigrid (open-BC) + fft_oracle (Stage 3)
    │   └── viz/                        # paper_repro static figures (Stage 4/4B); GGUI/movie = Stage 7
    ├── tests/                          # config/hello_sim/units/data/ic/deposit/solver_*/integrator/io/sim
    ├── docs/development.md             # contributor quickstart
    └── legacy/                         # original 2020 CUDA source, preserved for reference
        ├── README.md                   # build notes + bug pointers
        ├── CUDAfft2.0.cu               # 307 lines — FFT-Poisson TEST HARNESS
        ├── final_draft1.cu             # 632 lines — full simulation (primary)
        └── Multi-Parallel-CUDA-FFT.cu  # 619 lines — full simulation (OpenMP variant)
```

Git history shows earlier C-only prototypes (`Numerical_Differentiator`, `Potential_Array`,
`densA_particleA`, `Particle_Array`, test files) that were deleted — i.e. the project migrated
from a CPU/C prototype to the CUDA version, now preserved under `legacy/`.

### 3.2 The three files compared

| File | Role | Distinctive traits |
|---|---|---|
| `CUDAfft2.0.cu` | **Unit test** for the spectral Poisson solver. Places a few point masses in a 256³ grid, solves, prints potential values around them, renders a potential map. | No particles, no time loop. Static global grids. This is effectively the only existing "test." |
| `final_draft1.cu` | **Primary full simulation.** IC generation, density deposition, FFT-Poisson, integrator, time loop, image dumps. | `#pragma omp` directives present but mostly **inert** (no enclosing parallel region). |
| `Multi-Parallel-CUDA-FFT.cu` | **Later variant** of the full sim. ~90% identical to `final_draft1.cu`. | Defines a global `G`; `center_diff` uses a real `#pragma omp parallel for private(...)`; cleaner density deposit; extra image timesteps. Treat this as the most recent iteration. |

> The three files duplicate ~90% of their content. The modernization collapses them into a
> single library + entry points.

### 3.3 Function-by-function (full simulation files)

- **`FFT_poisson(den_array, grav_po)`** — Copies ρ into a flat array (folding in
  `4πG·m`), builds the wavenumber vector `k` (0…N/2 then −N/2…−1), uploads to device, runs
  `real2complex` → forward C2C FFT → `solve_poisson` (divide by −k²) → inverse C2C FFT →
  `complex2real_scaled` (×1/N³), downloads Φ. **Re-allocates device buffers and re-creates the
  cuFFT plan every call.**
- **`solve_poisson` / `real2complex` / `complex2real_scaled`** (CUDA `__global__` kernels) — the
  spectral solve in k-space and the real↔complex packing. The DC term (k=0) is special-cased to
  avoid divide-by-zero.
- **`densArray(particleArray, den_array)`** — Zeroes ρ then deposits each particle by
  **nearest-grid-point** (`ρ[round(x)][round(y)][round(z)] += 1`). Mass per particle is *not*
  applied here (folded into G later). Rebuilds the whole grid each step.
- **`center_diff(...)`** — The "integrator." Computes each galaxy's center of mass, then for every
  particle computes an acceleration from (a) central finite differences of Φ and (b) hard-coded
  **point-mass force terms** `4.9e-14 · 600 / R` toward each galaxy center (a kludge standing in
  for the central black holes, which carry no grid mass), and updates position/velocity.
- **`CM_finder(galaxy_ID, xyz, part_array)`** — Mean position of a galaxy's particle range.
- **`make_image(array, name, title)`** — Projects the central N/2³ sub-volume of the grid down
  the z-axis into a 2D image and renders a PNG via **DISLIN**.
- **`main()`** — Allocates the particle array, populates two galaxies (bulge + 10 radial rings
  each, with circular velocities), allocates the density grid, then loops `t = 0…500` (`dt = 1`,
  one step = 10⁴ yr) calling deposit → solve → integrate, dumping images at selected steps.

### 3.4 Key parameters & "magic numbers"

| Symbol / value | Meaning | Notes |
|---|---|---|
| `N = 256` | Grid side length; box is 256 kpc, cells are 1 kpc³ | Fixed at compile time |
| `M = 100,000,000` | Particle count | See §3.6 — unrunnable in current layout |
| 7 floats/particle | `[x, y, z, vx, vy, vz, galaxyID]` | galaxyID: 0.0 = MW, 1.0 = Andromeda |
| `dt = 1` | Timestep = 10⁴ yr | t=125 → 1.25 Myr in image captions |
| `block_size = 2×2×2` | CUDA threads/block = 8 | **Severe** under-utilization (warp = 32) |
| particle mass 840 M☉ | From paper | Not used in `densArray`; folded into G |
| `9.571212e-15` ≈ `1.139430e-17 × 840` | G × particle-mass, in grid units | "kpc³ / (M☉ · 10 kyr²)" per comments |
| `1190` | "fraction of central BH gravity"; sets circular speed `v = √(1190·G/R)` | Empirical |
| `4.9e-14 · 600 / R` | Hard-coded central point-mass force | Per-galaxy BH kludge |
| centers (96,96,128), (160,160,128) | Galaxy 1 / Galaxy 2 placement | See bug #8 (G2 uses wrong center for velocity) |

### 3.5 Build & dependencies (as-is)

From the source header comment:

```
nvcc final_draft1.cu -Xcompiler -fopenmp \
  -I/home/phyd57/N_Body1/9.2/include -L/home/phyd57/N_Body1/9.2/lib64 -lcufft \
  -o CUDAfftcu2.out -I/usr/local/dislin -ldislin
```

- **CUDA 9.2** (2018) + **cuFFT** — hard-coded lab paths.
- **DISLIN** — proprietary scientific plotting library at `/usr/local/dislin`; awkward to install,
  effectively unmaintained for this use, and the main portability blocker for visualization.
- **OpenMP** — used for CPU loop parallelism (mostly inert in `final_draft1.cu`).
- The *legacy* code had no Makefile/CMake, no tests, no config files, no snapshot I/O. (As of
  Stage 0 the repo as a whole now has `pyproject.toml`, `.gitignore`, CI, `pytest`, and YAML
  configs — see §3.1 — but none of that builds the legacy `.cu`; those remain reference-only.)

### 3.6 Known bugs & limitations

**Correctness — physics/logic:**

1. **Integrator is scrambled.** In `center_diff`, position/velocity components are written in the
   wrong order and overwritten (e.g. `particleArray[i][2]` is assigned twice; final `[0]=x`); the
   speed *magnitude* `√(vx²+vy²+vz²)` is mixed into the x-velocity update. This is not a clean
   kick-drift-kick leapfrog and is dimensionally inconsistent.
2. **Out-of-bounds grid access.** `grav_po[(int)x ± 1]…` has no bounds checks; particles drifting
   to/over the box edge index outside the array → garbage forces or segfault.
3. **`CM_finder` galaxy IDs are swapped.** Called with `0` and `1`, but the function branches on
   `==1` vs else, so the two galaxies' centers get computed from the wrong particle ranges.
4. **Galaxy 2 rotation center is wrong.** Velocity setup uses `X = pos − 96` for *both* galaxies;
   Galaxy 2 is centered at 160, so its initial rotation curve is built around the wrong center.
5. **Nearest-grid-point deposition** (not CIC). Causes the blocky/"square galaxy" artifact the
   paper noted. No interpolation on deposit *or* gather.
6. **Periodic boundary conditions** (inherent to a plain FFT solve). Physically wrong for two
   isolated galaxies in vacuum; contributes to the density "dissipating" over time. (See §5.4.)
7. **Finite-difference factor.** Central differences divided by 4 (rather than 2·Δx = 2),
   apparently to absorb an unexplained factor — units are ad hoc and inconsistent across files.

**Correctness — memory/C:**

8. **Particle layout is array-of-pointers** (`float**` + M separate `malloc`s). At M = 1e8 that's
   100M tiny allocations, ~2.8 GB of data plus ~0.8 GB of pointers, with catastrophic cache
   locality. Also `malloc(7 * sizeof(float*))` over-allocates (pointer size, not float size).
   **This single issue makes the headline config effectively unrunnable.**
9. **`cudaFree` on host pointers.** `FFT_poisson` calls `cudaFree(k_xyz/den/den_inital)` on
   `malloc`'d host memory (should be `free`), **frees `den` twice** (double free), and **never
   frees `den_d`** (leak).
10. **Broken zeroing loop in `densArray`:** `for(j=0; i<N; i++)` / `for(k=0; i<N; i++)` test and
    increment `i` instead of `j`/`k` — incorrect/again-buggy nested loop.
11. **`den_inital` copies only the first N·N elements** (should be N·N·N) and is then unused.

**Performance:**

12. **8 threads/block** (`2×2×2`) wastes 75%+ of every warp; kernels run far below capacity.
13. **No device-resident state.** Every step re-allocates device buffers, **re-creates the cuFFT
    plan**, and copies the full 67 MB grid host↔device twice. Enormous per-step overhead.
14. **`pow(x,2)` / `sqrt` everywhere** in hot loops instead of multiplies.
15. **OpenMP largely inert** in `final_draft1.cu` (`#pragma omp for` with no parallel region).

**Maintainability/portability:**

16. ~90% code duplication across three files; static globals; mixed `float***` and `float[N][N][N]`
    signatures for the same data; DISLIN + CUDA-9.2 hard-coded paths; no build system/tests/IO.

### 3.7 Code ↔ paper divergences (watch these when "reproducing")

- Paper describes a clean **leapfrog**; code's `center_diff` does not implement one (bug #1).
- Paper implies particle counts from physical masses (2.32×10¹¹ M☉ ÷ 840 M☉ ≈ 2.8×10⁸ per galaxy);
  code hard-codes M = 1×10⁸ total. Counts/masses are inconsistent — pick a documented convention.
- Paper's surface-density law (Gaia-derived, Eq. 1: Σ(r)=Σ₀ℓc/√((r−rb)²+ℓc²), rb=2 kpc, ℓc=2.5 kpc,
  Σ₀≈611 M☉/pc²) is **not** what the legacy code samples; *it* used uniform-in-disk radial rings.
  (The modernized IC generator *does* sample Eq. 1's shape — §5.5, D16 — so this is a legacy-only gap.)
- **Tracer particle (Stage 4 / 4B):** the paper tracks a special ~840 M☉ "Sun-like" particle at 8.32 kpc.
  Our run is uniform-mass, so `ic.select_tracer` *tracks an existing disk particle* nearest 8.32 kpc
  rather than inserting a distinct-mass one — the trajectory is what reproduces, not the tracer's mass.

---

## 4. Modernization goals & decisions (scoped with the owner)

### 4.1 Goals / non-goals

**Goals**
- **Research-grade physics**: correct units, CIC deposition + interpolation, energy-/momentum-aware
  kick-drift-kick leapfrog, force softening, open boundary conditions.
- **Central black holes** as massive, softened particles (one per galaxy) that carry real grid mass
  — fixing the 2020 force-only kludge.
- **One portable source** running on **CPU, NVIDIA CUDA, and Apple-Silicon GPU**, with **Apple GPU
  treated as a first-class performance target** (not a fallback).
- **Scale derived from physical galaxy mass**, landing in the **10–100M particle** working range
  (see §5.2); performance is non-negotiable.
- **All four output modes**: real-time 3D viewer, batch→movie, data+analysis snapshots, and faithful
  reproduction of the 2020 paper figures.

**Non-goals (for now)**
- **Dark matter halo** — deferred; it was not part of the 2020 project. Keep a clean hook for a
  static NFW halo, but it is out of scope for the initial modernization.
- TreePM / FMM / adaptive mesh (a larger rearchitecture — revisit after PM is solid).
- Hydrodynamics / gas (SPH or grid). Collisionless gravity only at first.
- Cosmological expansion / comoving coordinates. This is an isolated two-galaxy problem.

### 4.2 Decision record

| # | Decision | Choice | Rationale |
|---|---|---|---|
| D1 | Overall goal | Research-grade physics + visuals + full portability | Owner selected "all three" |
| D2 | Kernel strategy | **One portable source** (Taichi: Python-embedded kernels → native CPU/CUDA/Metal) | Owner wants single source at near-hand-tuned perf; lets us hand-write the CIC scatter |
| D3 | Scale target | **10–100M particles**, SoA layout mandatory | Owner; fits NVIDIA ≥12 GB and Apple unified ≥16 GB |
| D4 | Apple GPU priority | **First-class performance** | Owner |
| D5 | Poisson solver | **Pluggable**: open-BC **multigrid** = portable default; **zero-padded FFT** = NVIDIA validation oracle | FFT is the one piece Taichi can't do portably; multigrid is pure portable kernels and fixes the periodic-BC physics bug. FFT oracle gives spectral ground-truth for tests + paper repro |
| D6 | Precision | fp32 for compute on all backends; fp64 diagnostics on CPU/CUDA only | Apple Metal has **no hardware fp64** — a hard constraint |
| D7 | Dev hardware | **Apple M5 Pro Mac** is the primary dev box (CPU + Metal); a **CUDA NVIDIA workstation is not yet set up** — it is a Stage-5 prerequisite to configure together (§7). | Updated 2026-06-21; Apple hardware now on hand (resolves R1) |
| D8 | Particle count | **Derive N from physical galaxy mass**: N = M_galaxy / m_particle. Particle mass is the resolution knob, bounded by memory. Full-res (840 M☉/particle) ≈ 5.5×10⁸ total; coarsen to land in the 10–100M working range (§5.2). | Owner |
| D9 | Central black holes | **In scope** — one massive, softened particle per galaxy that carries real grid mass. | Owner |
| D10 | Dark matter | **Deferred** — not part of the 2020 project; keep a hook, revisit later. | Owner |
| D11 | Repo placement | **In-place on `master`** — modernized package at repo root, legacy `.cu` moved to `legacy/`, Apache-2.0 + git history kept | Owner (2026-06-16); single source of truth, preserves history (relicensed GPL-3.0 → Apache-2.0 in 0bbd01e) |
| D12 | Language / runtime | **Python ≥ 3.11** + Taichi kernels; `hatchling` build, `src/` layout | Owner (2026-06-16) |
| D13 | CI & lint | **GitHub Actions** (origin is GitHub): install → `ruff` lint → `pytest` → CPU `hello-sim` smoke run | Owner (2026-06-16) |
| D14 | Config format | **YAML** run configs (`pyyaml`), validated by the `SimConfig` schema | Owner (2026-06-16) |
| D15 | Unit system | **(kpc, Myr, M☉)**; single `G ≈ 4.498×10⁻¹² kpc³ M☉⁻¹ Myr⁻²` *derived* from the standard `4.30091×10⁻³ pc M☉⁻¹ (km/s)²` and unit-tested | Owner (2026-06-16); resolves §9 Q2. Aligns dt=0.01 Myr with the 2020 code's 10⁴ yr step |
| D16 | Stage-2 IC modeling | **Hernquist bulge (~10% of stellar mass)** + Gaia-profile disk; **same SMBH ≈ 1×10⁸ M☉ in both galaxies**; **cold circular disk** (v_c from softened enclosed mass); disk profile *shape* from paper Eq. 1 **normalized to the mass-derived target** | Owner (2026-06-16). Bulge dispersion + thin-disk rotation are virial/spherical approximations to be validated at Stage 3; per-galaxy distinct masses & warm disk are Stage-8 refinements |
| D17 | Apple Metal target | **Optimize Metal for the M5 *Pro* tier** (the owner's dev Mac). | Owner (2026-06-21); resolves §9 Q1. This is the in-hand validation/perf-target hardware |
| D18 | Stage-4 reproduction scope | **Higher-fidelity, CPU-only repro:** 256³ grid, **10–30M particles**, both 4v & 2v over ~400–450 Myr; **headline figures only** (collision density-projection sequence + Sun-like tracer trajectory). Production solver = **multigrid** (portable default, D5); FFT oracle as a spot-check. Done in **three review checkpoints** (4A validation hardening → 4B repro machinery → 4C production runs + write-up). | Owner (2026-06-21). Measured CPU cost (M5 Pro): ~0.6 s/step (multigrid) / ~1.3 s/step (FFT) at 256³; a ~900-step 30M-particle run is ~18 min, so the full repro is CPU-affordable — **no CUDA (Stage 5) dependency**. Static matplotlib figures only; GGUI/movie viz stays Stage 7 |

---

## 5. Target architecture

### 5.1 Overview & module map

A thin, well-separated pipeline with swappable backends and solvers:

```
            ┌───────────────┐
config ───► │  IC Generator │ (disk+bulge sampling, velocities, bulk approach)
            └───────┬───────┘
                    ▼
            ┌───────────────┐   per step (all device-resident):
            │   Simulator   │   kick½ → drift → deposit(CIC) → PoissonSolver → grad → kick½
            └───┬───────┬───┘
                │       │
        ┌───────▼──┐ ┌──▼─────────────────┐
        │ Deposit/ │ │  PoissonSolver IF  │── Multigrid (portable, default)
        │ Gather   │ │  (pluggable)       │── ZeroPaddedFFT (oracle: cuFFT/pocketfft)
        │ (CIC)    │ └────────────────────┘
        └──────────┘
                    │
            ┌───────▼────────┐      ┌──────────────┐
            │  Diagnostics   │      │   Snapshot   │ (HDF5/npz: pos,vel,ρ,Φ,energy)
            │ (E, p, L)      │      │     I/O      │
            └────────────────┘      └──────┬───────┘
                                           ▼
                          ┌────────────────────────────────┐
                          │ Visualization                  │
                          │  • realtime viewer (Taichi GGUI)│
                          │  • offline frames → ffmpeg movie│
                          │  • paper-repro plots (matplotlib)│
                          └────────────────────────────────┘
```

Everything between IC generation and snapshot output stays **resident on the device** across all
timesteps — host↔device transfers happen only for snapshots/visualization. (This alone fixes
bug #13.)

### 5.2 Data model & memory budget

- **Structure-of-Arrays (SoA)**, flat contiguous fields, never array-of-pointers:
  `pos_x[N], pos_y[N], pos_z[N], vel_x[N], vel_y[N], vel_z[N], mass[N], gid[N]`.
  **As implemented in Stage 1, `mass[N]` is included** (a per-particle f32) so the central black
  holes carry real mass (D9) — making the layout **32 B/particle**, not the 28 B of an early
  mass-less sketch. (Drop `mass` for a uniform-mass run and it would be 28 B; see §11 RV2.)
- Grid fields: `rho[256³]`, `phi[256³]`, the acceleration field `g=−∇Φ` (3×256³), plus the
  multigrid hierarchy (phi/rhs/res across levels). As implemented in Stage 3 this is
  **≈ 8.43 × 256³ f32 ≈ 0.57 GB** (see `data.py` `_GRID_FIELDS`; §11 RV4b).

**Memory at fp32 (production = multigrid, no FFT padding; particles at 32 B):**

| N particles | Particles (6×f32 + f32 mass + i32) | Grid (ρ, Φ, g, MG hierarchy) | Total (approx.) |
|---|---|---|---|
| 10M | ~0.32 GB | ~0.57 GB | **~0.9 GB** |
| 30M | ~0.96 GB | ~0.57 GB | **~1.5 GB** |
| 100M | ~3.2 GB | ~0.57 GB | **~3.8 GB** |

The **zero-padded FFT oracle** adds a 512³ complex buffer (~1.07 GB) — used only on NVIDIA for
validation, so it never constrains the Apple/production path.

**Particle count is derived from physical mass (D8).** Each galaxy's mass (paper value
M ≈ 2.32×10¹¹ M☉ within 25 kpc) divided by the per-particle mass sets N; the per-particle mass is
the resolution knob, bounded by memory:

| Per-particle mass | N (both galaxies) | Particle memory (32 B/particle) | Fits on |
|---|---|---|---|
| 840 M☉ (paper, full resolution) | ~5.5×10⁸ | ~17.6 GB | 24 GB+ NVIDIA / ≥32 GB Apple unified |
| ~4,640 M☉ | 100M | ~3.2 GB | ≥8–12 GB GPU |
| ~15,500 M☉ | 30M | ~0.96 GB | most GPUs |
| ~46,400 M☉ | 10M | ~0.32 GB | laptop / entry GPU |

So full physical resolution (the paper's own 840 M☉ particle) is ~5.5×10⁸ particles (~17.6 GB) — a
stretch target needing a 24 GB+ card; the 10–100M working range is reached by coarsening the
per-particle mass. (The paper assumed both galaxies have *equal* mass; real Andromeda is more massive, so using
true per-galaxy masses is a Stage-8 refinement.)

### 5.3 Backend strategy & precision policy

- **Taichi** as the kernel layer: write each kernel once with `@ti.kernel`; select
  `ti.init(arch=ti.cpu | ti.cuda | ti.metal)` at runtime. Compiles to native code per backend.
- **Hardware reporting (`sim.init_backend`):** every run prints the **detected device it is actually
  running on** (e.g. `Apple M5 Pro — Apple GPU (Metal)`, `NVIDIA <name>`, or `… (CPU)`). Taichi
  *silently* falls back to CPU when a requested GPU backend is unavailable (e.g. `cuda` on a Mac → the
  native `arm64` arch); `init_backend` detects that mismatch and prints a loud **WARNING** so a run
  never *looks* like it used a GPU it didn't. The resolved arch + device are returned in the run
  summary (`backend_resolved`, `device`).
- **Precision (D6):** all force/position/velocity math in **fp32** on every backend (positions live
  in grid units 0–256, so fp32 is plenty for forces). **fp64 only for diagnostics** (total
  energy/momentum drift), computed on CPU or CUDA; on Apple, diagnostics either run in fp32 with
  compensated (Kahan) summation or are offloaded to a CPU reduction. **Never** require fp64 inside a
  Metal kernel.
- **Determinism/portability test:** identical IC + seed must produce matching trajectories across
  CPU/CUDA/Metal within an fp32 tolerance (see §6).

### 5.4 Poisson solver (the §3 deep-dive, settled as D5)

The solver maps ρ → Φ via ∇²Φ = 4πGρ. The boundary condition matters more than the algorithm:

- **Periodic (a plain FFT):** the box tiles space → ghost copies of both galaxies pull on the real
  ones, and forces wrap across edges. This is the artifact source in the 2020 results. **Rejected.**
- **Open / isolated (target):** Φ → 0 far from the mass; each galaxy feels only itself + its partner.

Two implementations behind one `PoissonSolver` interface:

1. **Multigrid (portable default).** Real-space V-cycles (e.g. red-black Gauss-Seidel or Jacobi
   smoother, full-weighting restriction, trilinear prolongation) solving ∇²Φ = 4πGρ on the 256³
   grid. **Open BCs** are imposed by setting the box-face potential from a multipole expansion
   (monopole + low-order moments) of the enclosed mass, then relaxing the interior. Pure Taichi
   stencil kernels — one source, full speed on CPU/CUDA/Metal, **no FFT dependency**. This is what
   runs in production on all backends.
2. **Zero-padded FFT (NVIDIA validation oracle).** Embed ρ in a 512³ zero-padded box and convolve
   with the *isolated* Green's function (Hockney–Eastwood) to get spectrally exact open-BC
   potential. Backed by cuFFT on CUDA (and pocketfft/NumPy on CPU). Used to (a) unit-test multigrid
   accuracy and (b) reproduce the paper's spectral results exactly.

**Why multigrid is the right default here:** the solve is grid-bound (~constant cost regardless of
N), so at 100M particles it's only ~10–20% of a step — the O(N) deposition/gather dominate. That
means we should choose the solver for *correctness + portability*, and multigrid is both fully
portable (the one thing a plain FFT is not under Taichi) and physically correct (open BCs).

### 5.5 Physics correctness plan

- **Units:** adopt a single documented system — recommend **(kpc, Myr, M☉)**. In these units
  **G ≈ 4.498×10⁻¹² kpc³ M☉⁻¹ Myr⁻²** (equivalently the familiar **4.498×10⁻⁶ kpc³ M☉⁻¹ Gyr⁻²**,
  or **4.30091×10⁻³ pc M☉⁻¹ (km/s)²**); pin this down in a unit test. Note the original code's
  timestep `dt = 1` corresponds to 10⁴ yr (0.01 Myr). Replace all magic numbers with named,
  unit-checked constants.
- **Deposition + gather:** **Cloud-In-Cell (CIC)** — trilinear weights on both deposit and force
  interpolation. Fixes the "square galaxy" artifact and conserves momentum (symmetric weights).
- **Integrator:** **kick-drift-kick (KDK) leapfrog**, symplectic, second-order; this is what the
  paper intended. Global fixed `dt` to start; adaptive dt is a later option.
- **Softening:** the PM force is **grid-softened** at the cell scale by CIC + the finite grid (no
  explicit ε on the deposit→solve→grad→gather path). Explicit **Plummer softening** (ε) is used only
  by the *direct* O(N²) force (Kepler test, central BHs) and by the IC circular-velocity setup
  (`GalaxyModel.softening = 0.3 kpc`). *Stage-3 status:* reconciling that 0.3 kpc IC ε with the ~1 kpc
  effective grid softening (so disk ICs are launched on the same potential the grid feels) is a
  Stage-4 refinement.
- **Initial conditions (research-grade):** derive the particle count from physical mass —
  N = M_galaxy / m_particle (§5.2) — then sample the Gaia-derived surface density (paper Eq. 1) for a
  disk + bulge, set circular velocities from the enclosed-mass rotation curve, add a vertical scale
  height, and add the bulk approach velocities (the 4v / 2v runs). Place a **central black hole as a
  massive, softened particle** at each galaxy center so it carries real grid mass (fixing the 2020
  force-only kludge). A **dark-matter halo is deferred** (not in 2020 scope); leave a clean,
  documented hook to add a static NFW halo later.
- **Diagnostics:** track total energy, linear & angular momentum, and the tracer-particle path;
  energy drift is the headline correctness metric.

### 5.6 Performance plan

- **Device-resident state** for the whole run; transfers only for snapshots (fixes bug #13).
- **CIC deposition** is the hot kernel (atomic scatter, 8 adds/particle). Options to benchmark:
  plain `ti.atomic_add`; per-cell **privatization**; or **sort particles by cell index** (radix
  sort) to coalesce and cut atomic contention. This is where we earn the "performance
  non-negotiable" requirement.
- **Gather + KDK** are bandwidth-bound and embarrassingly parallel — straightforward.
- Persistent solver state (multigrid hierarchy buffers; cuFFT plan for the oracle) allocated once.
- Avoid `pow`/redundant `sqrt`; precompute reciprocals; fuse kernels where Taichi allows.
- **Benchmark matrix:** steps/sec at N ∈ {1M, 10M, 30M, 100M} on each backend, plus a per-stage
  breakdown (deposit / solve / gather / integrate).

### 5.7 Visualization & output

- **Real-time 3D viewer:** Taichi **GGUI** particle renderer (in-process, GPU-resident points) —
  rotate, scrub, tweak parameters live. (If a standalone app is later wanted, export to a viewer.)
- **Batch → movie:** headless run dumps frames (offscreen render or projected-density PNGs) →
  stitch with `ffmpeg` to MP4/GIF. **Replaces DISLIN entirely.**
- **Data + analysis:** periodic snapshots to **HDF5** (or `.npz`) — positions, velocities, ρ, Φ,
  and diagnostics — for offline plotting/analysis.
- **Reproduce the paper:** a script that recreates the density-projection panels and the
  tracked-particle trajectory using the FFT oracle and the original-style ICs.

### 5.8 Proposed repo layout

> **Stages 0–3 realized this layout** (see §3.1 for what exists on disk today): `pyproject.toml`,
> `configs/`, `src/galaxy_collision/{config,sim,units,data,ic,deposit,integrator,diagnostics,io}.py`,
> the `solver/` package (`base` + `multigrid` + `fft_oracle`), `tests/`, `docs/`, and `legacy/`.
> Stage 4 / 4B added `viz/paper_repro.py` (static paper-reproduction figures + the `paper-repro`
> CLI); the remaining `viz/` placeholder is the Stage-7 realtime GGUI viewer + ffmpeg movie path.

```
galaxy_collision/
├── pyproject.toml                  # deps: taichi, numpy, h5py, matplotlib, (cupy optional)
├── README.md  ·  AGENT.md  ·  LICENSE (Apache-2.0)
├── configs/                        # YAML/TOML run configs (N, dt, ICs, solver, backend)
├── src/galaxy_collision/
│   ├── units.py                    # unit system + G, with tests
│   ├── data.py                     # SoA particle + grid fields
│   ├── ic.py                       # initial conditions (disk/bulge/halo, 4v/2v runs)
│   ├── deposit.py                  # CIC deposit + gather kernels (Taichi)
│   ├── solver/
│   │   ├── base.py                 # PoissonSolver interface
│   │   ├── multigrid.py            # open-BC multigrid (portable default)
│   │   └── fft_oracle.py           # zero-padded isolated FFT (validation)
│   ├── integrator.py               # KDK leapfrog + softening
│   ├── diagnostics.py              # energy/momentum/tracer (fp64 on CPU/CUDA)
│   ├── io.py                       # HDF5/npz snapshots
│   ├── viz/ (realtime GGUI, offline frames, paper_repro)
│   └── sim.py                      # orchestration + CLI entry
├── tests/                          # see §6
└── legacy/                         # the original three .cu files, preserved for reference
```

---

## 6. Verification & testing strategy

Correctness is the gate for "research-grade." Tests, roughly in order of authority:

1. **Solver vs analytic point mass.** A single mass → Φ(r) ≈ −GM/r outside the cell; check radial
   profile. (The legacy `CUDAfft2.0.cu` already did a hand-version of this — formalize it.)
2. **Multigrid vs FFT oracle.** Same ρ through both solvers must agree within tolerance.
3. **Plummer sphere equilibrium.** A sampled Plummer/iso sphere should stay in steady state (no
   spurious expansion/collapse) — directly tests the "density dissipates" bug.
4. **Two-body Kepler orbit.** Period and energy match the analytic solution → validates integrator
   + softening.
5. **Conservation.** Total energy drift < threshold over a long run; momentum conserved to fp
   tolerance (CIC symmetry).
6. **Cross-backend determinism.** CPU vs CUDA vs Metal agree within fp32 tolerance on a fixed IC.
7. **Performance regression.** Benchmark matrix (§5.6) tracked over time.

For high-stakes correctness checks, run the comparison as an isolated verification pass (e.g. a
dedicated test job / subagent) rather than eyeballing plots.

**On "energy drift < threshold" (test 5):** the grid PE ½ΣρΦ·dV is a *PM-tolerance* metric, not
the exact integrated Hamiltonian — it carries a dt-independent grid-discretization floor *plus* a
fluctuating CIC self-energy term that makes its drift realization-fragile (it crosses 1% on ~2/5
seeds). So the **exit-gate energy check uses the direct softened pair-sum PE** (softening = 1 cell,
matching the grid's effective softening), which omits the self-energy noise and is robust: across 10
seeds drift is **~0.1–0.9%** (FFT ~0.3%, multigrid worst ~0.9%) (gate < 2%, a named constant in
`tests/test_sim.py`, giving ~2× margin plus cross-platform fp32 headroom).
The grid PE is still *recorded* in the run history for production monitoring of large runs (where the
O(N²) direct PE is infeasible). The trajectory is symplectic regardless — the Kepler test proves that
on the direct path, where energy and force share one Hamiltonian.

---

## 7. Staged implementation plan

The full, shareable version of this plan — with per-stage objectives, tasks, deliverables, and exit
criteria — is in **`Galaxy_Collision_Modernization_Plan.docx`** (a generated artifact, gitignored;
this section is the tracked source of truth).

> **Current status (2026-06-21):** Stages 0–3 ✅ done. The FFT oracle was pulled forward into
> Stage 3 (so multigrid is validated against it now). **Stage 4 in progress** (scope = D18) —
> **4A ✅ done, 4B ✅ done** (machinery); 4C (production runs) next — structured as three review
> checkpoints:
> - **4A — Validation hardening ✅:** principled multigrid≈oracle tolerance + analytic point-mass
>   multigrid test (RV9); reconcile the 0.3 kpc IC velocity softening with the ~1 kpc grid
>   softening so disks launch in equilibrium (§5.5); a big-run grid-energy diagnostic that omits
>   the CIC self-energy offset (RV11); two-galaxy integration smoke test.
> - **4B — Repro machinery ✅:** Sun-like tracer particle (~8.32 kpc galactocentric) selection +
>   path recording (`run_simulation(tracer_indices=…)`); matplotlib density-projection sequence +
>   tracer-trajectory figures (`viz/paper_repro.py`, the DISLIN `make_image` replacement, static
>   only); the `paper-repro` CLI driver + `configs/paper_{4v,2v}.yaml`.
> - **4C — Production runs + write-up:** 4v & 2v at 256³ / 10–30M particles, headline figures,
>   qualitative comparison to the 2020 paper; flip this row to ✅. **Validation bar (explicit):**
>   at 10–30M the O(N²) direct PE is infeasible, so 4C has **no hard energy-conservation gate** —
>   conservation rests on grid-energy *drift* (a monitor, not pass/fail; RV7/RV11) + the virial
>   trend. 4C validation is therefore **qualitative** (morphology figures + tracer path + virial
>   consistency), not a numeric energy gate. (A cold thin disk also heats/spreads over many t_dyn
>   in a PM code — mitigated by the chosen resolution; a warm disk is the Stage-8 fix, D16.)

| Stage | Outcome | Backend | Exit gate | Status |
|---|---|---|---|---|
| **0 — Scaffold & guardrails** | Repo, build, CI, config schema; legacy `.cu` archived | — | CI green; a trivial `hello-sim` runs | ✅ **Done** (2026-06-16) |
| **1 — Units & data model** | Unit system (G) + SoA particle/grid fields | CPU | G unit test passes; config round-trips | ✅ **Done** (2026-06-16) |
| **2 — Initial conditions** | Mass-derived N, disk+bulge sampling, central BH, 4v/2v setups | CPU | ICs match target mass/profile; reproducible from seed | ✅ **Done** (2026-06-16) |
| **3 — CPU reference (anchor)** | A *correct* sim: CIC + open-BC multigrid + KDK + softening + diagnostics + I/O | CPU | Plummer stays stable; two-body Kepler matches; energy drift < threshold | ✅ **Done** (2026-06-20) |
| **4 — Validation & FFT oracle** | Zero-padded isolated FFT (✅ done in Stage 3) + full test suite + paper reproduction | CPU/CUDA | Multigrid ≈ FFT within tolerance (✅ test 2 passing); paper figures reproduced (⬜) | 🟡 partial |
| **5 — CUDA & scale-up** | Device-resident state, deposition tuning, 100M+ runs | CUDA | 100M-particle run; benchmark + per-stage profile | ⬜ **Prereq:** a CUDA NVIDIA workstation + Taichi CUDA runtime is **not yet set up** (D7) — configure it together as the first task when Stage 5 begins, before any CUDA work |
| **6 — Apple / Metal** | First-class Apple GPU, fp32 compute policy | Metal | Cross-backend parity (test 6); perf benchmark vs CUDA | ⬜ |
| **7 — Visualization & output** | Realtime GGUI, batch→movie, paper figures, tracer particle | all | All four output modes working | ⬜ |
| **8 — Research campaigns** | 4v/2v studies, central-BH experiments, longer runs | all | Reproduce + extend 2020 results (future hooks: DM halo, TreePM) | ⬜ |

Stage 3 is the linchpin: get the physics provably correct on CPU first, *then* chase speed and
portability on top of a trusted reference. Stages 5 and 6 each must **re-pass the Stage 3–4 test
suite** to count as done.

---

## 8. Risks & mitigations

- **R1 — ~~No Apple hardware on hand~~ → RESOLVED (2026-06-21).** The owner's primary dev box is now
  an **Apple M5 Pro** Mac, which is exactly the Metal target tier (D17) — Metal validation/benchmarking
  runs directly on it (Taichi `arch=metal` confirmed working). The flipped risk is now **CUDA access**:
  the NVIDIA workstation is not yet set up and is a Stage-5 prerequisite (D7, §7).
- **R2 — Taichi Metal feature/perf gaps** (e.g. large atomic-scatter throughput, kernel features).
  *Mitigation:* a Phase-0/1 spike that benchmarks CIC deposition on Metal; if it underperforms,
  fall back plan is **Vulkan-compute + VkFFT** (one GPU source covering NVIDIA + Apple via
  MoltenVK/native-Metal, peak FFT everywhere) — keep the `PoissonSolver`/kernel interfaces
  abstract enough to swap.
- **R3 — fp32 accuracy over long integrations.** *Mitigation:* energy-drift monitoring,
  compensated summation for diagnostics, optional fp64 reference run on CUDA.
- **R4 — Multigrid open-BC accuracy** (boundary multipole order, convergence). *Mitigation:* it's
  continuously checked against the FFT oracle (test 2); tune V-cycles/tolerance to hit target.

---

## 9. Open questions for the owner

**Resolved 2026-06-16:** particle count is derived from physical galaxy mass (D8); a massive central
black hole is in scope (D9); dark matter is deferred (D10). Stage-0 toolchain settled: in-place on
`master` (D11), Python ≥ 3.11 + Taichi (D12), GitHub Actions CI (D13), YAML configs (D14) — this
resolves former Q4 (repo placement). Unit system standardized on **(kpc, Myr, M☉)** with a single
derived, unit-tested G (D15) — this resolves former Q2.

**Resolved 2026-06-21:** former Q1 (Apple hardware) — the owner's dev Mac is an **Apple M5 Pro**;
Metal is validated/optimized directly on it at the **M5 Pro tier** (D17), resolving Risk R1. The
counterpart is that the CUDA NVIDIA workstation is not yet provisioned and is a Stage-5 prerequisite
(D7, §7 Stage-5 row).

Still open:

1. **Real-time viewer:** is an in-process Taichi GGUI window sufficient, or do you want a standalone
   app / web viewer?

---

## 10. References

- **Paper:** Adams, Lajevardi, LeBlanc, Movaseghi, *Galaxy Collisions With CUDA and FFT*, PHYD57,
  Apr 15 2020. (Authoritative statement of intent; see §2.1, §3.7.)
- Hockney & Eastwood, *Computer Simulation Using Particles* — PM method, CIC, zero-padded
  (isolated) FFT Green's functions.
- Binney & Tremaine, *Galactic Dynamics* (2nd ed.) — units, equilibria, ICs (cited in the paper).
- Springel, GADGET / TreePM lecture notes — the gold-standard next step beyond PM.
- VkFFT (Tolmachev) — cross-platform FFT (CUDA/Vulkan/Metal/OpenCL), MIT — the R2 fallback FFT.
- Taichi documentation — `ti.cpu / ti.cuda / ti.metal` backends, GGUI renderer.

---

## 11. Review log & deferred decisions

Items surfaced by the **Stage 0–1 code review (2026-06-16)**. None of these block the Stage 0/1
exit criteria — which are met (units/config tests pass, `ruff` clean, G independently verified, the
scaffold + CI are sound). Each is a small fix or a doc-sync left for the owner to schedule.
`Status: Pending` = not yet decided/applied. (These use `RV…` IDs to mark them as review findings,
distinct from the `D…` decision record in §4.2.) RV1–RV3 and RV4(a) were resolved 2026-06-16.

| ID | Area | Finding | Suggested action | Status |
|---|---|---|---|---|
| RV1 | Test hygiene | `tests/test_data.py` put `pytest.importorskip("taichi")` at **module level, below** its pure-Python memory-estimator tests, so in any environment without Taichi the **whole file skipped**. Those tests have no Taichi dependency and they guard the §5.2 memory math. | Move the `importorskip` into the two Taichi-only test functions. | ✅ **Done** (2026-06-16) — `importorskip` is now per-test; the estimator tests always run. |
| RV2 | Doc ↔ code drift | `data.py` uses **32 B/particle** (it includes a per-particle `mass` field for the central black holes, D9), but §5.2 previously stated **28 B** / "100M → 2.8 GB". True footprint is ~3.2 GB at 100M. | Keep the mass field and update the §5.2 tables to 32 B (option **a**). | ✅ **Done** (2026-06-16) — §5.2 tables + prose now state 32 B (100M → ~3.2 GB; full-res → ~17.6 GB). |
| RV3 | Packaging / compat | `pyproject.toml` set `requires-python = ">=3.11"`, but Taichi wheels lag the newest CPython (no 3.13 wheels yet), so `pip install` could fail on 3.13. | Cap to `>=3.11,<3.13`. | ✅ **Done** (2026-06-16) — cap applied in `pyproject.toml`. |
| RV4 | Nits | (a) `estimate_memory` accepted `n_particles == 0` while `ParticleState`/`SimConfig` require positive — a minor inconsistency. (b) `_GRID_FIELDS = 2` is a Stage-1 placeholder; it must grow when the multigrid hierarchy (Stage 3) and the zero-padded FFT buffer (Stage 4) land — the docstring already flags this. | (a) Align the zero-handling. (b) Revisit grid-memory accounting at Stage 3/4. | (a) ✅ **Done** (2026-06-16). (b) ✅ **Done** (2026-06-20) — `_GRID_FIELDS` now counts ρ+Φ (2) + accel (3) + MG hierarchy (≈3.43) ≈ 8.43 G³; §5.2 grid column updated to ~0.57 GB. FFT-oracle pad stays excluded (validation-only). |

Items surfaced by the **Stage 3 adversarial verification (2026-06-20)** — deferred, non-blocking
(the Stage-3 exit criteria are met: Plummer stable, two-body Kepler matches, energy drift small,
multigrid ≈ FFT oracle):

| ID | Area | Finding | Suggested action | Status |
|---|---|---|---|---|
| RV5 | Perf (Stage 5) | The open-BC multigrid converges at a healthy but non-optimal ~0.5–0.7/cycle (measured on small/64³ blobs). The power-of-two vertex-centered coarsening puts the coarse far-face a half-cell inside the fine one, which slows (does not break) convergence; correctness comes from the fine-grid operator + multipole BCs, so the converged Φ is right regardless. | Tune at Stage 5 (e.g. proper 2^L+1 nesting or W-cycles/Galerkin coarsening); fine for a CPU reference. | ⬜ Deferred to Stage 5 |
| RV6 | Perf (Stage 5) | Device-residency gaps for the CPU reference: (a) the multigrid boundary multipole moments (M, CM, quadrupole) are computed host-side in NumPy each solve (a host↔device hop); (b) diagnostics pull full fields to host via `to_numpy()` each sample; (c) grid state is split across owners — `GridState` holds only ρ/Φ, the acceleration grids are allocated ad-hoc in `run_simulation`, and the MG hierarchy lives inside the solver. Clear and correct for the reference, but legacy bug #13 territory. | Move the moment reduction into a Taichi kernel; keep diagnostics device-side (Kahan fp32 on Metal); consolidate grid ownership when chasing device-resident performance. | ⬜ Deferred to Stage 5 |

Items from the **exit-gate hardening review (2026-06-20)** — the FFT energy gate was found
realization-fragile; fixed, plus an oracle memory fix and two deferrals logged:

| ID | Area | Finding | Suggested action | Status |
|---|---|---|---|---|
| RV7 | Test robustness | The Plummer exit-gate energy drift was measured with the grid PE (½ΣρΦ·dV), whose fluctuating CIC self-energy term made it realization-fragile — it crossed the 1% gate on ~2/5 seeds (verified). | Gate on the **direct softened pair-sum PE** (softening = 1 cell) instead; ~0.1–0.9% drift across 10 seeds/both solvers, gate < 2%. Grid PE still recorded for production monitoring. | ✅ **Done** (2026-06-20) — `run_simulation(direct_pe_softening=…)`; gate in `tests/test_sim.py`. |
| RV8 | Perf / memory | `fft_oracle._build_greens_ft` used `np.meshgrid(d,d,d)`, materializing three (2N)³ f64 arrays (~3.2 GB transient at 256³, exceeding the "~1 GB" the §5.2 oracle note implies). | Build the squared-distance grid by broadcasting (one (2N)³ array). | ✅ **Done** (2026-06-20) — broadcasting in `_build_greens_ft`. |
| RV9 | Test margin (Stage 4) | The multigrid≈FFT-oracle agreement (test 2) sits at ~1.7% vs the 3% cap on a centered Gaussian fixture; fixture-sensitive. (Stage 3 added an off-center case, also ~1.7%.) | Re-derive a principled tolerance and/or add a compact/point-mass multigrid-vs-analytic case at Stage 4. | ✅ **Done** (2026-06-21, 4A) — added `test_multigrid_point_mass_far_field_monopole` (deviation obeys the *derived* 0.26→0.4/r² discrete-Green's law, monotone) + a scale-free `_resid_ratio` convergence gate; the fixture-sensitive MG≈oracle cap is demoted to a documented loose sanity check (the tight checks are convergence + the analytic anchor). |
| RV10 | Coverage (Stage 5) | Production-scale multigrid (n_cycles=20 @ 256³) has no residual-convergence assertion — only small grids are tested for convergence. | Add a 256³ residual-norm regression once the CUDA path makes it cheap. | ⬜ Deferred to Stage 5 |
| RV11 | Diagnostics (Stage 4) | The grid PE's self-energy offset (RV7) follows into Stage-4 paper-reproduction energy diagnostics on large runs, where the O(N²) direct PE is infeasible. | Use a self-energy-subtracted grid energy estimate (or accept the offset and report drift, not absolute E) for big-run diagnostics. | ✅ **Done** (2026-06-21, 4A) — took the *report-drift* path (evidence-based): probes showed the CIC cloud self-energy is sub-percent while the grid-PE↔direct-PE gap is a ~2–4% *discretization* offset (locked by `tests/test_diagnostics.py`), so a self-energy subtraction would change absolute E by ~0.1% and not shrink the drift. Big-run conservation is the grid-energy **drift** (a *monitoring* metric, not a gate, per RV7) plus the offset-insensitive **virial** ratio (`virial_ratio`/`total_energy_grid` helpers, virial recorded in the run history). |

*Verified good in the same review:* the (kpc, Myr, M☉) unit system and derived G (independently
re-derived to 4.4985×10⁻¹²; circular-velocity sanity check 231.9 km/s for 10¹¹ M☉ at 8 kpc),
config validation + YAML round-trip, `ruff` lint, and the overall scaffold/CI.

Items from the **Stage 4 / 4B review (2026-06-22)** — paper-repro machinery; non-blocking:

| ID | Area | Finding | Suggested action | Status |
|---|---|---|---|---|
| RV12 | Viz units (4B) | The density figures labeled the colorbar `Σ [M☉/kpc²]` but `viz/paper_repro._hist2d` returned mass *per bin*, not per kpc² (off by a constant ~1.4× at the default 300 bins over a 256 kpc box). The image is unaffected (LogNorm), but a paper-comparison figure should plot the labeled quantity. | Divide the histogram by the bin area so it is a true surface density. | ✅ **Done** (2026-06-22) — `_hist2d` divides by the (exact, from edges) bin area; label is now correct, image unchanged. |
| RV13 | Tracer sampling (4C) | The Sun-like tracer path is recorded at `history_cadence`, so a smooth orbit needs a fine history cadence even when only a few *density* snapshots are wanted — the two samplings are coupled. | If 4C wants a smooth tracer path with sparse density snapshots, decouple by adding a separate `tracer_cadence` to `run_simulation`. Controllable today via `history_cadence`; only worth it if the coupling bites. | ⬜ Deferred (optional, 4C) |

---

*Document status: architecture + staged plan, scoped with the owner on 2026-06-16 (decisions
D1–D15). Companion deliverable: `Galaxy_Collision_Modernization_Plan.docx` (gitignored generated
artifact; this file is the tracked source of truth). Update the Decision Record (§4.2), Open
Questions (§9), and the Review Log (§11) as choices are finalized during implementation.*
