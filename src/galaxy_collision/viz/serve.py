"""Headless web-stream viewer (Stage 7 follow-on) — watch a live run in a browser via a link.

The desktop ``galaxy-view`` needs a display; this serves the live collision over **HTTP** instead,
so a run on a **headless** box (e.g. the CUDA workstation) is viewable in any browser. It streams
the live **2D density projection** as **MJPEG** (`multipart/x-mixed-replace`) — the browser shows it
as a self-updating image. The projection is pure compute + an on-device colormap (``viz.project`` +
``viz.viewer`` kernels), so **no GGUI/display/Vulkan context is needed** — it works anywhere the sim
runs (CPU/CUDA/Metal).

Usage::

    galaxy-serve --config configs/paper_4v.yaml --backend cuda            # serves on 127.0.0.1:8080
    # from your laptop, tunnel the headless box's port and open the link:
    ssh -N -L 8080:localhost:8080 user@cuda-box      # then browse http://localhost:8080

``--host 0.0.0.0`` exposes it on the network directly (no auth — only on a trusted network; the
SSH-tunnel above is the safe default). Browser buttons pause/resume, change speed, and cycle the
projection plane. Needs the ``viz`` extra (``pip install -e ".[viz]"``) for JPEG encoding.

This module deliberately omits ``from __future__ import annotations`` (consistency with the kernels
it imports; see ``deposit.py``).
"""

import io
import math
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import numpy as np

from galaxy_collision.viz.project import scatter_density
from galaxy_collision.viz.viewer import (
    _LOG_DECADES,
    _LUT_N,
    _PLANE_NAMES,
    _PROJ_PLANES,
    _colormap_log,
    _inferno_lut,
    _max_field,
    _setup_chain,
)

_PAGE = """<!doctype html><html><head><meta charset=utf-8><title>Galaxy Collision (live)</title>
<style>body{margin:0;background:#000;color:#ccc;font:14px system-ui;text-align:center}
img{max-width:100vw;max-height:88vh;object-fit:contain}
button{background:#222;color:#ccc;border:1px solid #444;padding:6px 12px;margin:3px;
border-radius:5px;cursor:pointer}
button:hover{background:#333}#bar{padding:6px}</style></head><body>
<img src="/stream" alt="live density projection">
<div id=bar>
<button onclick="c('pause')">&#9199; pause/resume</button>
<button onclick="c('slower')">&laquo; slower</button>
<button onclick="c('faster')">faster &raquo;</button>
<button onclick="c('plane')">&#8635; plane (xy/xz/yz)</button>
<span id=s></span></div>
<script>function c(a){fetch('/control?action='+a).then(r=>r.text())
.then(t=>document.getElementById('s').textContent=' '+t)}</script>
</body></html>"""


def _frame_jpeg(rgb_xy: np.ndarray) -> bytes:
    """Encode a colormapped ``[x, y]`` device image (floats 0–1) as JPEG bytes, oriented for a
    browser (row 0 = top, y increasing upward — matching the density figures' origin-lower view)."""
    import imageio.v3 as iio

    img = (np.clip(rgb_xy, 0.0, 1.0) * 255.0).astype(np.uint8)
    img = np.transpose(img, (1, 0, 2))[::-1]  # [x,y]→[y,x], flip so +y is up
    buf = io.BytesIO()
    iio.imwrite(buf, img, extension=".jpg")
    return buf.getvalue()


def serve_viewer(
    config,
    *,
    host: str = "127.0.0.1",
    port: int = 8080,
    bins: int = 512,
    steps_per_frame: int = 1,
    fps_cap: float = 30.0,
    _frames_limit: int = 0,  # test hook: stop the sim loop after N frames (0 = run forever)
):
    """Run the sim headless and stream its live 2D density projection as MJPEG over HTTP.

    Blocks (serving) until interrupted. Returns a small summary dict (mainly for the test hook)."""
    import taichi as ti

    from galaxy_collision.deposit import deposit_density, gather_acceleration, potential_to_accel
    from galaxy_collision.integrator import kdk_step
    from galaxy_collision.sim import init_backend

    backend = init_backend(config)
    icr, parts, grid, acc, solver, n, gs, dx = _setup_chain(config)
    acc_x, acc_y, acc_z = acc

    proj_img = ti.field(ti.f32, shape=(bins, bins))
    rgb_img = ti.Vector.field(3, ti.f32, shape=(bins, bins))
    lut = ti.field(ti.f32, shape=(_LUT_N, 3))
    lut.from_numpy(_inferno_lut())
    maxbuf = ti.field(ti.f32, shape=())
    inv_logspan = 1.0 / (_LOG_DECADES * math.log(10.0))
    extent_full = (0.0, float(gs), 0.0, float(gs))

    warm = {"on": False}

    def accel_fn():
        deposit_density(parts, grid.rho, dx)
        solver.solve(grid.rho, grid.phi, warm_start=warm["on"])
        warm["on"] = True
        potential_to_accel(grid.phi, grid.ax, grid.ay, grid.az, dx)
        gather_acceleration(parts, grid.ax, grid.ay, grid.az, acc_x, acc_y, acc_z, dx)

    accel_fn()  # prime acc_* for the opening half-kick

    # Shared state: the latest JPEG frame + the controls the browser drives. Plain dict reads/writes
    # are atomic enough under the GIL for these scalars; all Taichi calls stay on the main thread.
    state = {"frame": None, "paused": False, "spf": max(1, steps_per_frame), "plane": 0,
             "step": 0, "fps_cap": fps_cap}

    def render_frame():
        scatter_density(parts, proj_img, extent_full, axes=_PROJ_PLANES[state["plane"]])
        _max_field(proj_img, maxbuf)
        vmax = float(maxbuf[None]) or 1.0
        log_vmin = math.log(vmax) - _LOG_DECADES * math.log(10.0)
        _colormap_log(proj_img, lut, log_vmin, inv_logspan, rgb_img)
        state["frame"] = _frame_jpeg(rgb_img.to_numpy())

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *_):  # keep the console quiet
            pass

        def _send(self, code, ctype, body=b""):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            if body:
                self.wfile.write(body)

        def do_GET(self):
            path = urlparse(self.path).path
            if path in ("/", "/index.html"):
                self._send(200, "text/html; charset=utf-8", _PAGE.encode())
            elif path == "/control":
                q = parse_qs(urlparse(self.path).query)
                action = (q.get("action") or [""])[0]
                if action == "pause":
                    state["paused"] = not state["paused"]
                elif action == "slower":
                    state["spf"] = max(1, state["spf"] - 1)
                elif action == "faster":
                    state["spf"] += 1
                elif action == "plane":
                    state["plane"] = (state["plane"] + 1) % len(_PROJ_PLANES)
                msg = (f"step {state['step']} | {_PLANE_NAMES[_PROJ_PLANES[state['plane']]]} | "
                       f"{state['spf']} steps/frame | {'paused' if state['paused'] else 'running'}")
                self._send(200, "text/plain; charset=utf-8", msg.encode())
            elif path == "/stream":
                self.send_response(200)
                self.send_header("Content-Type", "multipart/x-mixed-replace; boundary=frame")
                self.send_header("Cache-Control", "no-store")
                self.end_headers()
                try:
                    while True:
                        frame = state["frame"]
                        if frame is not None:
                            self.wfile.write(b"--frame\r\nContent-Type: image/jpeg\r\n"
                                             + f"Content-Length: {len(frame)}\r\n\r\n".encode()
                                             + frame + b"\r\n")
                        time.sleep(1.0 / max(1.0, state["fps_cap"]))
                except (BrokenPipeError, ConnectionResetError):
                    pass
            else:
                self._send(404, "text/plain", b"not found")

    server = ThreadingHTTPServer((host, port), Handler)
    server.daemon_threads = True
    threading.Thread(target=server.serve_forever, daemon=True).start()
    shown = host if host != "0.0.0.0" else "<this-host>"
    print(f"galaxy-serve: live view on {backend['device']} at http://{shown}:{port}/  "
          f"(Ctrl-C to stop)")
    if host == "0.0.0.0":
        print("  WARNING: bound to 0.0.0.0 — the stream is reachable on the network with no auth.")

    render_frame()  # first frame before any client connects
    try:
        while True:
            if not state["paused"]:
                for _ in range(state["spf"]):
                    kdk_step(parts, acc_x, acc_y, acc_z, accel_fn, config.dt)
                    state["step"] += 1
            render_frame()
            if _frames_limit and state["step"] >= _frames_limit:
                break
            time.sleep(1.0 / max(1.0, state["fps_cap"]))
    except KeyboardInterrupt:
        print("\ngalaxy-serve: stopped.")
    finally:
        server.shutdown()

    return {"device": backend["device"], "backend_resolved": backend["resolved"],
            "n_particles": n, "steps_run": state["step"], "host": host, "port": port}


def serve_cli(argv: list[str] | None = None) -> int:
    """CLI: ``galaxy-serve --config C [--backend B] [--host H] [--port P] [--bins N]
    [--steps-per-frame K] [--fps F]``."""
    import argparse

    from galaxy_collision.config import load_config

    p = argparse.ArgumentParser(
        prog="galaxy-serve",
        description="Stream a live PM run to a browser over HTTP (headless; no display needed).",
    )
    p.add_argument("--config", type=__import__("pathlib").Path, required=True)
    p.add_argument("--backend", default=None, help="cpu | cuda | metal (overrides the config).")
    p.add_argument("--host", default="127.0.0.1",
                   help="bind address (0.0.0.0 to expose on the network; default localhost).")
    p.add_argument("--port", type=int, default=8080)
    p.add_argument("--bins", type=int, default=512, help="density-projection resolution.")
    p.add_argument("--steps-per-frame", type=int, default=1)
    p.add_argument("--fps", type=float, default=30.0,
                   help="frame cap (the sim still steps at full speed).")
    args = p.parse_args(argv)

    config = load_config(args.config)
    if args.backend is not None:
        config.backend = args.backend
        config.validate()

    serve_viewer(config, host=args.host, port=args.port, bins=args.bins,
                 steps_per_frame=args.steps_per_frame, fps_cap=args.fps)
    return 0


if __name__ == "__main__":
    raise SystemExit(serve_cli())
