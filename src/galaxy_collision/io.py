"""Snapshot I/O: periodic dumps of the simulation state (Stage 3).

AGENT.md §5.7: periodic snapshots to **HDF5** (or ``.npz``) holding positions,
velocities, mass, galaxy id, optionally the ρ/Φ grids, plus the run-time diagnostics —
the offline record for analysis, movies, and paper reproduction. This replaces the
legacy code's DISLIN PNG dumps as the *data* output (rendering is Stage 7).

This module is deliberately **Taichi-free**: it operates on plain NumPy arrays bundled in
a :class:`Snapshot`, so it imports and tests without a GPU runtime. The simulator converts
its device-resident fields to NumPy (``field.to_numpy()``) via :func:`snapshot_from_states`
at the configured ``output_cadence`` and hands the result here.

Format is chosen by file extension: ``.h5``/``.hdf5`` → HDF5 (h5py), ``.npz`` → NumPy.
Both round-trip a :class:`Snapshot` losslessly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Particle datasets always present; grids optional.
_PARTICLE_KEYS = ("pos", "vel", "mass", "gid")
_GRID_KEYS = ("rho", "phi")
# Scalar metadata carried as attributes.
_META_KEYS = ("step", "time", "dx")


@dataclass
class Snapshot:
    """One simulation snapshot as NumPy arrays (internal units: kpc, Myr, M_sun)."""

    step: int
    time: float  # Myr
    pos: np.ndarray  # (N, 3) kpc
    vel: np.ndarray  # (N, 3) kpc/Myr
    mass: np.ndarray  # (N,) M_sun
    gid: np.ndarray  # (N,) int32
    dx: float = 1.0  # cell size, kpc
    rho: np.ndarray | None = None  # (G, G, G) M_sun/kpc^3, optional
    phi: np.ndarray | None = None  # (G, G, G) potential, optional
    diagnostics: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return self.pos.shape[0]


def snapshot_from_states(step, time, parts, grid=None, diagnostics=None, dx=1.0) -> Snapshot:
    """Build a :class:`Snapshot` by pulling a Taichi ``ParticleState`` (+ optional grid) to host.

    Imports no Taichi itself — it just calls ``.to_numpy()`` on the passed fields.
    """
    pos = np.stack(
        [parts.pos_x.to_numpy(), parts.pos_y.to_numpy(), parts.pos_z.to_numpy()], axis=1
    )
    vel = np.stack(
        [parts.vel_x.to_numpy(), parts.vel_y.to_numpy(), parts.vel_z.to_numpy()], axis=1
    )
    rho = grid.rho.to_numpy() if grid is not None else None
    phi = grid.phi.to_numpy() if grid is not None else None
    diag = {k: np.asarray(v) for k, v in (diagnostics or {}).items()}
    return Snapshot(
        step=int(step),
        time=float(time),
        pos=pos,
        vel=vel,
        mass=parts.mass.to_numpy(),
        gid=parts.gid.to_numpy(),
        dx=float(dx),
        rho=rho,
        phi=phi,
        diagnostics=diag,
    )


def write_snapshot(snap: Snapshot, path: str | Path) -> Path:
    """Write a snapshot to ``path``; HDF5 or npz chosen by suffix. Returns the path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix in (".h5", ".hdf5"):
        _write_hdf5(snap, path)
    elif suffix == ".npz":
        _write_npz(snap, path)
    else:
        raise ValueError(f"unsupported snapshot extension {path.suffix!r} (use .h5/.hdf5/.npz)")
    return path


def read_snapshot(path: str | Path) -> Snapshot:
    """Read a snapshot from ``path``; HDF5 or npz chosen by suffix."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"snapshot not found: {path}")
    suffix = path.suffix.lower()
    if suffix in (".h5", ".hdf5"):
        return _read_hdf5(path)
    if suffix == ".npz":
        return _read_npz(path)
    raise ValueError(f"unsupported snapshot extension {path.suffix!r} (use .h5/.hdf5/.npz)")


# --- HDF5 backend ---------------------------------------------------------------


def _write_hdf5(snap: Snapshot, path: Path) -> None:
    import h5py

    with h5py.File(path, "w") as f:
        for key in _META_KEYS:
            f.attrs[key] = getattr(snap, key)
        for key in _PARTICLE_KEYS:
            f.create_dataset(key, data=getattr(snap, key), compression="gzip")
        for key in _GRID_KEYS:
            arr = getattr(snap, key)
            if arr is not None:
                f.create_dataset(key, data=arr, compression="gzip")
        if snap.diagnostics:
            grp = f.create_group("diagnostics")
            for key, val in snap.diagnostics.items():
                grp.create_dataset(key, data=np.asarray(val))


def _read_hdf5(path: Path) -> Snapshot:
    import h5py

    with h5py.File(path, "r") as f:
        kwargs = {key: f.attrs[key] for key in _META_KEYS}
        for key in _PARTICLE_KEYS:
            kwargs[key] = f[key][()]
        for key in _GRID_KEYS:
            kwargs[key] = f[key][()] if key in f else None
        diag = {}
        if "diagnostics" in f:
            # np.asarray so scalars come back as 0-d arrays — matching the npz backend
            # (h5py would otherwise return a bare numpy scalar). Lossless *and* consistent.
            diag = {k: np.asarray(f["diagnostics"][k][()]) for k in f["diagnostics"]}
    return Snapshot(
        step=int(kwargs["step"]),
        time=float(kwargs["time"]),
        pos=kwargs["pos"],
        vel=kwargs["vel"],
        mass=kwargs["mass"],
        gid=kwargs["gid"],
        dx=float(kwargs["dx"]),
        rho=kwargs["rho"],
        phi=kwargs["phi"],
        diagnostics=diag,
    )


# --- npz backend ----------------------------------------------------------------

_DIAG_PREFIX = "diag__"


def _write_npz(snap: Snapshot, path: Path) -> None:
    arrays = {
        "step": np.asarray(snap.step),
        "time": np.asarray(snap.time),
        "dx": np.asarray(snap.dx),
    }
    for key in _PARTICLE_KEYS:
        arrays[key] = getattr(snap, key)
    for key in _GRID_KEYS:
        arr = getattr(snap, key)
        if arr is not None:
            arrays[key] = arr
    for key, val in snap.diagnostics.items():
        arrays[_DIAG_PREFIX + key] = np.asarray(val)
    np.savez_compressed(path, **arrays)


def _read_npz(path: Path) -> Snapshot:
    with np.load(path) as data:
        diag = {
            k[len(_DIAG_PREFIX) :]: np.asarray(data[k])
            for k in data.files
            if k.startswith(_DIAG_PREFIX)
        }
        return Snapshot(
            step=int(data["step"]),
            time=float(data["time"]),
            pos=data["pos"],
            vel=data["vel"],
            mass=data["mass"],
            gid=data["gid"],
            dx=float(data["dx"]),
            rho=data["rho"] if "rho" in data.files else None,
            phi=data["phi"] if "phi" in data.files else None,
            diagnostics=diag,
        )


__all__ = ["Snapshot", "snapshot_from_states", "write_snapshot", "read_snapshot"]
