"""A tiny terminal progress bar for headless runs (Stage 7 follow-on).

Fills 0→100% over a known number of steps, drawn in place with a carriage return on stderr.
Deliberately dependency-free and near-zero overhead: it redraws **only when the integer percent
changes** (≤101 redraws over a whole run), so the per-step cost is one integer divide + compare.

It is **auto-enabled only when the stream is a TTY**, so piped/redirected output, logs, and the test
suite never see it — a library call stays silent unless a human is watching a terminal.
"""

from __future__ import annotations

import sys
import time


class ProgressBar:
    """Step-based progress bar. ``update(step)`` each iteration; ``finish()`` ends the line.

    ``enabled=None`` (default) auto-detects a TTY on ``stream``; pass ``True``/``False`` to force.
    Disabled (a no-op) when not a TTY or when ``total <= 0``.
    """

    def __init__(self, total, *, label="", enabled=None, width=28, stream=None):
        self.total = int(total)
        self.label = f"{label} " if label else ""
        self.width = width
        self.stream = stream if stream is not None else sys.stderr
        auto = self.stream.isatty() if hasattr(self.stream, "isatty") else False
        self.enabled = (auto if enabled is None else bool(enabled)) and self.total > 0
        self._last_pct = -1
        self.t0 = time.perf_counter()

    def update(self, step: int) -> None:
        if not self.enabled:
            return
        step = min(int(step), self.total)
        pct = step * 100 // self.total
        if pct == self._last_pct and step != self.total:
            return  # nothing visibly changed — skip the redraw
        self._last_pct = pct
        filled = self.width * step // self.total
        bar = "█" * filled + "░" * (self.width - filled)
        elapsed = time.perf_counter() - self.t0
        rate = step / elapsed if elapsed > 0 else 0.0
        tail = ""
        if rate > 0:
            tail = f"  {rate:5.1f} it/s  ETA {self._fmt((self.total - step) / rate)}"
        self.stream.write(f"\r{self.label}[{bar}] {pct:3d}%  {step}/{self.total}{tail}")
        self.stream.flush()

    def finish(self) -> None:
        """Terminate the bar's line with a newline (so following output starts fresh)."""
        if self.enabled:
            self.stream.write("\n")
            self.stream.flush()

    @staticmethod
    def _fmt(seconds: float) -> str:
        s = int(seconds)
        if s >= 3600:
            return f"{s // 3600}h{(s % 3600) // 60:02d}m"
        if s >= 60:
            return f"{s // 60}m{s % 60:02d}s"
        return f"{s}s"


__all__ = ["ProgressBar"]
