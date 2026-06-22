"""Visualization & output (per AGENT.md §5.7).

``paper_repro`` (Stage 4 / 4B) provides the **static** paper-reproduction figures — the
projected-density collision sequence and the Sun-like tracer path (the DISLIN ``make_image``
replacement), plus the ``paper-repro`` CLI driver. The realtime Taichi GGUI viewer and the
offline frames -> ffmpeg movie pipeline remain Stage 7. Replaces the legacy DISLIN dependency
entirely.
"""
