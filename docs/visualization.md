# Visualization & output

The four ways to get pictures and data out of a run (AGENT.md §5.7, Stage 7). Two were always
here (snapshots, paper figures); Stage 7 added the **movie** and the **realtime viewer**.

| Mode | Command | Headless? | Best run on |
|---|---|---|---|
| Realtime 3D/2D viewer | `galaxy-view` | No — needs a display | **Mac/desktop** (Metal) |
| **Live web view (browser link)** | `galaxy-serve` | **Yes** — streams over HTTP | **headless box** (e.g. CUDA) |
| Batch → movie (MP4/GIF) | `galaxy-movie` | Yes | CUDA (faster at scale) or Mac |
| Data snapshots (HDF5/npz) | `galaxy-sim` | Yes | anywhere |
| Paper-figure reproduction | `paper-repro` | Yes | anywhere |

All four share one rule: **the heavy data never leaves the GPU unnecessarily.** The viewer's 2D mode
and the movie both render from a *device-side* 2D density projection (`viz/project.py`) — a
deposit-style scatter that bins particle mass into a small `(bins, bins)` image on the GPU, so only
~1 MB/frame crosses to the host instead of the ~2.4 GB it would cost to pull 100M positions.

---

## 1. Realtime viewer — `galaxy-view`

Opens a window and animates the collision **as it computes**, in-process on the same GPU the physics
runs on (verified on Apple **Metal**; also CUDA/Vulkan). No extra install — it uses Taichi's built-in
GGUI renderer — but it **needs a display**, so it's a desktop/Mac tool, not a headless-CI one.

```bash
galaxy-view --config configs/paper_4v.yaml --backend metal
# lighter/smoother while exploring:
galaxy-view --config configs/paper_4v.yaml --backend metal --max-points 200000 --steps-per-frame 2
```

**Two modes** (toggle with **M**):

- **3D** — particles as GPU points, colored by galaxy (Milky Way blue, Andromeda orange), orbit
  camera. Only a fixed **subsample** of `--max-points` is drawn (drawing 10–100M points is neither
  feasible nor legible); the physics still runs at the config's full N.
- **2D** — the full-N device density projection, log-scaled through an on-GPU colormap. Cycle the
  projection plane (xy/xz/yz) with **P**.

**Controls:** RMB-drag orbit · WASD move · **SPACE** pause · **N** single-step · **R** restart to t=0 ·
**[ / ]** slower/faster (steps per frame) · **- / =** smaller/larger points (3D) · **M** 2D↔3D ·
**P** projection plane · **ESC/Q** quit. An on-screen overlay shows step / time / mode / speed.

**Interactive frame rate** tracks the sim's steps/s at the config's N (see
[`performance.md`](performance.md)): ~1M particles is fluid, 10–30M is watchable, 100M is a slideshow
(use `--steps-per-frame` to trade smoothness for evolution speed, or a smaller N to explore). The
*render* cost is bounded by `--max-points`, so it's the *physics* per step that sets the pace.

**Flags:** `--config` `--backend` `--max-points` `--steps-per-frame` `--bins` (2D resolution)
`--radius` (3D point size; 0 = auto). Headless quick-check: `--offscreen --frames N --out DIR` dumps
PNGs with no window (the CI smoke path).

---

## 1b. Live web view — `galaxy-serve` (headless, via a link)

The desktop `galaxy-view` needs a display. When the sim runs on a **headless** box (e.g. the CUDA
workstation), `galaxy-serve` streams the **live 2D density projection** over HTTP instead, so you
watch the collision in a **browser** — no display, no X server, no Vulkan context (it's pure compute
+ an on-device colormap, the same projection the movie uses). Needs the `viz` extra (for JPEG).

```bash
# on the headless box:
galaxy-serve --config configs/paper_4v.yaml --backend cuda            # serves 127.0.0.1:8080
# from your laptop, tunnel the port and open the link in any browser:
ssh -N -L 8080:localhost:8080 user@cuda-box                           # then: http://localhost:8080
```

The page shows a self-updating image (MJPEG `multipart/x-mixed-replace`) plus buttons to
**pause/resume**, change **speed** (steps/frame), and cycle the **projection plane** (xy/xz/yz). The
sim runs at full speed on the box; the browser just receives frames.

**Flags:** `--config` `--backend` `--host` (default `127.0.0.1`) `--port` (8080) `--bins` (projection
resolution) `--steps-per-frame` `--fps` (frame cap).

**Access & safety:** the **SSH tunnel above is the recommended way** to get a link — it's encrypted
and needs no open ports. `--host 0.0.0.0` instead binds the stream to the network directly (anyone
who can reach the box can watch — **no authentication**), so use it only on a trusted network. This
is a *view* of the 2D projection, not the interactive 3D window (that's `galaxy-view`, on a display).

---

## 2. Batch → movie — `galaxy-movie`

Runs a collision headless and encodes a density-projection **movie**. Needs the `viz` extra:

```bash
pip install -e ".[viz]"          # imageio + imageio-ffmpeg (bundles ffmpeg — no system install)
galaxy-movie --config configs/paper_4v.yaml --backend cuda --frame-cadence 5 --out collision.mp4
galaxy-movie --config configs/paper_4v.yaml --frame-cadence 5 --out collision.gif --panel final.png
```

Each frame is the device 2D projection → matplotlib `LogNorm` (one global color scale across the
whole movie, so brightness is comparable frame-to-frame, matching `paper_repro`'s surface-density
convention) → encoded to **MP4 or GIF** (by the `--out` extension). `--panel PATH` also saves the
final frame as a standalone PNG.

**Flags:** `--config` `--backend` `--frame-cadence` (steps between frames) `--bins` (frame resolution)
`--axes {xy,xz,yz}` `--fps` `--out` (`.mp4`/`.gif`) `--panel`. Without the `viz` extra it degrades
gracefully — writes the PNG frames + prints the `ffmpeg` command to stitch them.

**Where to render:** the movie path is headless, so the CUDA workstation renders frames fastest at
large N (Metal's CIC deposit is its throughput ceiling at scale — see `performance.md`/RV20). At
≤10M either backend is fine; for a fresh 100M movie, prefer CUDA.

---

## 3. Data snapshots — `galaxy-sim`

The analysis path: periodic **HDF5** (or `.npz`) snapshots of positions, velocities, ρ, Φ, and the
diagnostics history, for offline plotting. Driven from the YAML config:

```bash
galaxy-sim --config configs/plummer.yaml          # writes to outputs/ per output_cadence
```

Set `output_cadence` (steps between snapshots; `0` disables) and `output_dir` in the config. A 100M
snapshot is ~2.4 GB, which is exactly why the movie path projects on-device instead of dumping these.

---

## 4. Paper-figure reproduction — `paper-repro`

Recreates the 2020 paper's signature static figures — the collision density-projection sequence and
the Sun-like tracer-particle trajectory — with matplotlib (headless Agg):

```bash
paper-repro --config configs/paper_4v.yaml --out figures/
```

Flags: `--config` `--backend` `--tracer-radius` (galactocentric kpc of the tracer; default 8.32, the
Sun's distance) `--out`. See [`paper_reproduction.md`](paper_reproduction.md) for the science write-up.

---

## Headline artifacts

The 100M density-projection panels (the RV19 scale headline, rendered via the device projection at
the full particle count) live in [`figures/`](figures/):
`density_100M_t175_pericenter.png` (first pericenter) and `density_100M_t375_tidal.png` (tidal
re-expansion).

## Which machine?

- **Interactive viewer** (`galaxy-view`) → the Mac (or any box with a display); the Metal-native
  deliverable with full 3D + controls.
- **Watching a headless run via a link** (`galaxy-serve`) → run it on the headless box (CUDA), open
  the SSH-tunnelled URL in your browser. Live 2D view, no display needed on the box.
- **Movie / figures / snapshots** → headless, run anywhere; the CUDA box is fastest for large-N
  movie frame rendering, the Mac is fine for everything at ≤10M.
