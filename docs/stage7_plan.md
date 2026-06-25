# Stage 7 — Visualization & output (DRAFT plan)

> **Status:** draft for owner review (2026-06-25), not yet folded into `AGENT.md` §7.
> Builds on Stages 0–6 (CPU reference → paper repro → CUDA → Metal; all three backends live).

## 0. Goal & exit gate (from AGENT.md §5.7, §7)

> **Outcome:** the four output modes all working — (1) realtime 3D viewer, (2) batch→movie,
> (3) data+analysis snapshots, (4) paper reproduction.
> **Exit gate:** all four output modes working.

**Two of the four already exist** and just get re-confirmed: **(3) snapshots** (HDF5/npz, `io.py`,
Stage 3) and **(4) paper reproduction** (static figures, `viz/paper_repro.py`, Stage 4). So Stage 7's
*new* work is **(1) the realtime GGUI viewer** and **(2) the batch→movie pipeline**, plus the RV19
100M density panel.

Three owner decisions scoped this draft (2026-06-25):

- **D25 (proposed) — Realtime viewer form (resolves §9 Q1):** an **in-process Taichi GGUI window**
  (`ti.ui.Window`) — live, GPU-resident, interactive; no new deps. (A standalone/web viewer was
  rejected as out-of-scope for the Taichi-native path.)
- **D26 (proposed) — Viewer content:** **both** a subsampled **3D particle cloud** (orbit camera)
  **and** a full-N **2D density-projection** mode, toggled with a key.
- **D27 (proposed) — Batch→movie:** **density-projection frames → MP4/GIF via `imageio-ffmpeg`**
  (bundles ffmpeg — none is installed on the CUDA box, so no system dependency). Deliver 4v & 2v
  collision movies + the RV19 100M density panel.

Stage 7 closes **RV19** (100M figure) and optionally **RV13** (tracer-cadence decoupling).

## 1. The one idea that makes it scale: device-side 2D projection

Rendering a density frame at 100M by pulling all positions to the host is **~2.4 GB/frame** — a
non-starter for a movie. Instead, a small Taichi kernel **bins particles into a 2D image on the
device** (a deposit-style scatter into an `(H, W)` float field, one axis projected out), and only
that image (~1 MB at 512²) crosses to the host. This single kernel powers **both** new modes:

- the **movie**: per-frame device projection → PNG (matplotlib `LogNorm`, reusing `paper_repro`'s
  colormap/units) → `imageio-ffmpeg`;
- the **viewer's 2D mode**: device projection → `canvas.set_image` (no host round-trip needed for
  display on GPU backends).

It also retires the snapshot-bloat problem: movies never write the 2.4 GB snapshots, just frames.

## 2. Checkpoints (mirrors 5A/5B/5C, 6A/6B/6C)

### 7A — Batch→movie (headless; fully testable on the CUDA box)
*The robust, display-free half — and it builds the device projection kernel 7B reuses.*

1. **Device 2D projection kernel** (`viz/project.py`): scatter particle mass into an `(H, W)`
   field along a chosen projection (xy/xz/yz), CIC or NGP; returns the small image. Unit-tested vs
   a NumPy `histogram2d` reference (like the deposit tests).
2. **Frame renderer + encoder** (`viz/movie.py`): project → matplotlib density PNG (shared
   colormap/`LogNorm`/surface-density units with `paper_repro._hist2d`) → encode to MP4/GIF via
   `imageio-ffmpeg`. Consistent global color scale across frames.
3. **`galaxy-movie` CLI**: run a sim (or read snapshots) and emit a movie at a chosen cadence;
   optional Sun-like-tracer overlay. If the smooth-tracer/sparse-frame coupling bites, add
   `tracer_cadence` to `run_simulation` (RV13).
4. **Deliverables:** a 4v and a 2v collision movie (10–30M), and the **RV19 100M
   density-projection panel** (a strong headline still).

**Exit 7A:** `galaxy-movie` produces an MP4/GIF headless on CPU and CUDA; projection kernel tested;
100M panel rendered.

### 7B — Realtime GGUI viewer ✅ **Done (2026-06-25)**
*The interactive half. Code written/linted/offscreen-tested here; the live window is exercised on the M5 Pro.*

> **Feasibility spike first (the gating unknown).** Taichi GGUI was confirmed working on this
> **M5 Pro on `arch=metal`** (and on `vulkan`) — offscreen `show_window=False` renders a real
> particle image — so the sim *and* the renderer both run on the Apple GPU (no sim-on-metal /
> render-on-vulkan split), and the offscreen path is CI-smoke-testable.
>
> **Delivered.** (1) `viz/viewer.py` + `galaxy-view` CLI: an in-process `ti.ui.Window` driving the
> live force chain (deposit→solve→grad→gather→KDK), runnable on cpu/cuda/**metal**. (2) **3D mode** —
> `scene.particles` over a fixed ceil-capped subsample (`--max-points`), colored by galaxy id
> (MW blue / Andromeda orange), orbit camera. (3) **2D mode** — the §1 device projection
> (`viz.project`, the kernel the movie reuses) → LogNorm colormap → `canvas.set_image`; **M** toggles
> 3D↔2D, **P** cycles the plane (xy/xz/yz). (4) **Controls** — SPACE pause · N single-step · R restart
> to t=0 (reloads the stashed IC, no realloc) · `[`/`]` speed · `-`/`=` point size · ESC/Q quit · an
> on-screen text overlay (step/time/mode/speed). (5) **Offscreen mode** (`--offscreen --frames`) for a
> headless quick-check / CI. (6) `tests/test_viewer.py` — offscreen render smoke test (skips where
> there's no graphics device, like the CUDA tests); ruff clean; CPU+Metal suite green.
> *Interactive paths (keys, 2D `set_image`, text overlay) are validated on the Mac; the offscreen 3D
> render + the 2D projection/colormap path are verified here.*

1. **`viz/viewer.py` + `galaxy-view` CLI**: an in-process `ti.ui.Window` driving the live sim loop.
2. **3D mode:** `scene.particles` over a **fixed subsample** (~0.5–2M indices chosen once),
   colored by galaxy id / speed; orbit camera. Interactive at any N via the subsample.
3. **2D mode:** the §1 device projection → `canvas.set_image`; toggle 3D↔2D with a key.
4. **Controls:** pause / single-step / restart, a speed (steps-per-frame) control, and a couple of
   live knobs (e.g. point size, projection axis). Shows step/time/diagnostics as on-screen text.

**Exit 7B:** the viewer runs interactively on the Mac (Metal) at an interactive frame rate for a
10–30M sim; both modes + controls work. (Headless CI only imports/constructs it — no window.)

### 7C — Integration, headline artifacts, write-up
1. **All four modes confirmed** end-to-end (viewer, movie, snapshots, paper figures).
2. **Headline artifacts:** the 100M density panel (RV19) and a polished 4v collision movie; drop
   them in `docs/` and reference from the README.
3. **Docs/sync:** a short `docs/visualization.md` (how to drive viewer/movie); update `AGENT.md`
   §5.7/§7 (Stage 7 → ✅) + §9 (Q1 resolved) + §11; README "Highlights"/Usage; `pyproject` `viz`
   extra. Suite green on CPU/CUDA/Metal.

**Exit 7C (= Stage-7 exit gate):** all four output modes working + documented; suite green.

## 3. Dependencies & packaging

- New optional extra `viz = ["imageio>=2.31", "imageio-ffmpeg>=0.4"]` (movie encoding; ffmpeg
  bundled). `matplotlib` is already a core dep (figures). GGUI ships with Taichi (no new dep).
- The viewer needs a display → it's a Mac/desktop deliverable; the movie path is headless and runs
  on the CUDA box and in CI.

## 4. Risks & mitigations

- **Headless offscreen 3D render** (a 3D-particle *movie*) is unreliable on a display-less Linux box
  (needs a GL/Vulkan context) → **out of scope**; the movie uses the device 2D projection (D27), the
  3D view is interactive-only (Mac). (Matches the rejected "3D-particle render movie" option.)
- **GGUI interactivity at scale** → the 3D mode renders a *subsample*, not all N (D26); 2D mode is
  full-N but image-based (cheap).
- **Per-frame cost at 100M** → the device projection kernel keeps host transfer at ~1 MB/frame, not
  ~2.4 GB (§1).
- **`imageio-ffmpeg` wheel** carries a bundled ffmpeg; if a user declines the `viz` extra,
  `galaxy-movie` emits frames + the ffmpeg command (graceful degrade).

## 5. Out of scope (Stage 8+)

Science campaigns (4v/2v parameter studies, central-BH experiments, DM halo, TreePM) — Stage 8.
A standalone/web viewer (D25 rejected). On-the-fly volumetric ray-marching of the density field.
