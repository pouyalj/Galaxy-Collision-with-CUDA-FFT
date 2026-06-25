"""The terminal progress bar (Stage 7 follow-on): renders to 100%, bounded redraws, TTY-gated."""

from __future__ import annotations

import io

from galaxy_collision.progress import ProgressBar


def test_progress_bar_renders_to_full():
    buf = io.StringIO()
    bar = ProgressBar(40, label="run", enabled=True, stream=buf, width=10)
    for s in range(1, 41):
        bar.update(s)
    bar.finish()
    out = buf.getvalue()
    assert "100%" in out
    assert "█" in out  # a filled bar was drawn
    assert "40/40" in out
    assert out.endswith("\n")  # finish() terminates the line


def test_progress_bar_redraws_bounded_by_percent():
    """Redraws only when the integer percent changes → ≤101 over any run (cheap)."""
    buf = io.StringIO()
    bar = ProgressBar(10_000, enabled=True, stream=buf)
    for s in range(1, 10_001):
        bar.update(s)
    assert buf.getvalue().count("\r") <= 101


def test_progress_bar_auto_disabled_off_tty():
    """enabled=None auto-detects: a non-TTY stream (StringIO) stays silent."""
    buf = io.StringIO()
    bar = ProgressBar(50, enabled=None, stream=buf)
    assert not bar.enabled
    for s in range(1, 51):
        bar.update(s)
    bar.finish()
    assert buf.getvalue() == ""  # nothing written when disabled


def test_progress_bar_disabled_when_total_zero():
    buf = io.StringIO()
    bar = ProgressBar(0, enabled=True, stream=buf)
    assert not bar.enabled
    bar.update(0)
    assert buf.getvalue() == ""
