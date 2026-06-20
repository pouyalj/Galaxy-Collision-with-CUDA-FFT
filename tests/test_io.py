"""Snapshot I/O round-trips (Stage 3, Checkpoint D).

HDF5 and npz must each round-trip a Snapshot losslessly — particles, optional ρ/Φ grids,
diagnostics, and scalar metadata. The pure-NumPy path needs no Taichi; snapshot_from_states
is tested separately and self-skips without Taichi.
"""

from __future__ import annotations

import numpy as np
import pytest

from galaxy_collision.io import Snapshot, read_snapshot, write_snapshot


def _make_snapshot(rng, n=200, g=8, with_grid=True, with_diag=True):
    return Snapshot(
        step=12,
        time=0.12,
        pos=rng.uniform(0, g, size=(n, 3)).astype(np.float32),
        vel=rng.normal(size=(n, 3)).astype(np.float32),
        mass=rng.uniform(1.0, 2.0, size=n).astype(np.float32),
        gid=rng.integers(0, 2, size=n).astype(np.int32),
        dx=1.5,
        rho=(rng.random((g, g, g)).astype(np.float32) if with_grid else None),
        phi=(rng.random((g, g, g)).astype(np.float32) if with_grid else None),
        diagnostics=(
            {
                "energy_total": np.float64(-1.23e10),
                "momentum": np.array([0.1, -0.2, 0.3]),
                "ang_momentum": np.array([1.0, 2.0, 3.0]),
            }
            if with_diag
            else {}
        ),
    )


def _assert_equal(a: Snapshot, b: Snapshot):
    assert a.step == b.step
    assert a.time == pytest.approx(b.time)
    assert a.dx == pytest.approx(b.dx)
    # Values AND dtypes must survive the round-trip (assert_array_equal ignores dtype,
    # so a silent f32->f64 upcast would otherwise pass while doubling file size).
    for key in ("pos", "vel", "mass", "gid"):
        va, vb = getattr(a, key), getattr(b, key)
        np.testing.assert_array_equal(va, vb)
        assert va.dtype == vb.dtype, f"{key} dtype drift: {va.dtype} -> {vb.dtype}"
    for key in ("rho", "phi"):
        ga, gb = getattr(a, key), getattr(b, key)
        assert (ga is None) == (gb is None)
        if ga is not None:
            np.testing.assert_array_equal(ga, gb)
            assert ga.dtype == gb.dtype
    assert set(a.diagnostics) == set(b.diagnostics)
    for key in a.diagnostics:
        # Always an ndarray on read, regardless of backend (consistency guarantee).
        assert isinstance(b.diagnostics[key], np.ndarray)
        np.testing.assert_allclose(a.diagnostics[key], b.diagnostics[key])
        assert np.asarray(a.diagnostics[key]).dtype == b.diagnostics[key].dtype


@pytest.mark.parametrize("ext", [".h5", ".hdf5", ".npz"])
def test_round_trip_full(tmp_path, ext):
    rng = np.random.default_rng(0)
    snap = _make_snapshot(rng)
    path = write_snapshot(snap, tmp_path / f"snap{ext}")
    assert path.is_file()
    _assert_equal(snap, read_snapshot(path))


@pytest.mark.parametrize("ext", [".h5", ".npz"])
def test_round_trip_minimal_no_grid_no_diag(tmp_path, ext):
    rng = np.random.default_rng(1)
    snap = _make_snapshot(rng, with_grid=False, with_diag=False)
    path = write_snapshot(snap, tmp_path / f"min{ext}")
    out = read_snapshot(path)
    _assert_equal(snap, out)
    assert out.rho is None and out.phi is None and out.diagnostics == {}


def test_write_creates_parent_dirs(tmp_path):
    rng = np.random.default_rng(2)
    snap = _make_snapshot(rng, with_grid=False, with_diag=False)
    path = write_snapshot(snap, tmp_path / "nested" / "deeper" / "snap.npz")
    assert path.is_file()


def test_unsupported_extension(tmp_path):
    rng = np.random.default_rng(3)
    snap = _make_snapshot(rng, with_grid=False, with_diag=False)
    with pytest.raises(ValueError, match="unsupported snapshot extension"):
        write_snapshot(snap, tmp_path / "snap.txt")


def test_read_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_snapshot(tmp_path / "absent.h5")


def test_snapshot_from_states(tmp_path):
    pytest.importorskip("taichi", reason="Taichi not installed")
    import taichi as ti

    from galaxy_collision.data import GridState, ParticleState
    from galaxy_collision.io import snapshot_from_states

    ti.init(arch=ti.cpu)
    n, g = 50, 8
    parts = ParticleState(n)
    rng = np.random.default_rng(4)
    pos = rng.uniform(0, g, size=(n, 3)).astype(np.float32)
    parts.pos_x.from_numpy(np.ascontiguousarray(pos[:, 0]))
    parts.pos_y.from_numpy(np.ascontiguousarray(pos[:, 1]))
    parts.pos_z.from_numpy(np.ascontiguousarray(pos[:, 2]))
    parts.mass.from_numpy(rng.uniform(1.0, 2.0, size=n).astype(np.float32))
    parts.gid.from_numpy(rng.integers(0, 2, size=n).astype(np.int32))
    grid = GridState(g)
    grid.rho.from_numpy(rng.random((g, g, g)).astype(np.float32))

    diag = {"energy_total": np.float64(-5.0)}
    snap = snapshot_from_states(7, 0.07, parts, grid=grid, diagnostics=diag, dx=2.0)
    assert snap.n == n
    np.testing.assert_array_equal(snap.pos[:, 0], pos[:, 0])
    assert snap.rho is not None and snap.phi is not None

    out = read_snapshot(write_snapshot(snap, tmp_path / "fromstate.h5"))
    np.testing.assert_array_equal(out.pos, snap.pos)
    assert out.diagnostics["energy_total"] == pytest.approx(-5.0)
