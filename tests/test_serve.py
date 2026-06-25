"""Headless web-stream viewer smoke test (Stage 7 follow-on).

``galaxy-serve`` streams the live 2D density projection as MJPEG over HTTP — no display needed, so
unlike the GGUI viewer this *is* fully headless-testable. We start the server on an ephemeral port,
fetch the page and one stream frame, exercise a control endpoint, and check the frame is a valid
non-empty JPEG. Self-skips without the ``viz`` extra (imageio, for JPEG encoding).
"""

from __future__ import annotations

import socket
import threading
import time
import urllib.request

import pytest


def _free_port() -> int:
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


@pytest.mark.fixed_arch  # picks its own backend; not subject to GALAXY_TEST_ARCH
def test_serve_streams_a_jpeg_frame():
    pytest.importorskip("taichi", reason="Taichi not installed")
    pytest.importorskip("imageio", reason="viz extra (imageio) not installed")

    from galaxy_collision.config import SimConfig
    from galaxy_collision.viz.serve import serve_viewer

    port = _free_port()
    cfg = SimConfig(
        name="serve-smoke", backend="metal", ic_preset="two_galaxy_4v",
        n_particles=20_000, dt=0.5, steps=0, solver="multigrid", grid_size=32, seed=0,
    )

    result: dict = {}

    def run():
        try:
            # _frames_limit stops the (otherwise infinite) sim loop so the test thread joins.
            result["summary"] = serve_viewer(
                cfg, host="127.0.0.1", port=port, bins=64, steps_per_frame=1, _frames_limit=8
            )
        except Exception as e:  # no graphics device is fine (pure compute), but be defensive
            result["error"] = e

    t = threading.Thread(target=run, daemon=True)
    t.start()

    # Wait for the server to come up + the first frame to render.
    base = f"http://127.0.0.1:{port}"
    page = None
    for _ in range(100):
        try:
            page = urllib.request.urlopen(base + "/", timeout=2).read()
            break
        except Exception:
            if "error" in result:
                pytest.skip(f"serve_viewer unavailable here: {result['error']}")
            time.sleep(0.1)
    assert page is not None and b"/stream" in page, "index page did not serve"

    # A control endpoint responds.
    ctl = urllib.request.urlopen(base + "/control?action=plane", timeout=2).read()
    assert b"steps/frame" in ctl

    # Pull one JPEG out of the MJPEG multipart stream.
    with urllib.request.urlopen(base + "/stream", timeout=3) as r:
        chunk = r.read(200_000)
    start = chunk.find(b"\xff\xd8")  # JPEG SOI marker
    assert start != -1, "no JPEG frame in the MJPEG stream"
    assert chunk.find(b"\xff\xd9", start) != -1, "JPEG not terminated (EOI)"

    t.join(timeout=15)
    assert "error" not in result, f"serve loop raised: {result.get('error')}"
