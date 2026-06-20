"""Kick-drift-kick (KDK) leapfrog integrator + Plummer-softened forces (Stage 3).

AGENT.md §5.5: a symplectic, second-order **KDK leapfrog** is the integrator the paper
intended (the legacy ``center_diff`` was a scrambled, dimensionally-inconsistent mess —
bug #1). One full step at fixed ``dt``:

    kick(½dt) → drift(dt) → [recompute acceleration] → kick(½dt)

The acceleration between the kicks is recomputed by the PM force chain (deposit → solve
→ ``potential_to_accel`` → gather; see ``deposit.py`` and ``solver/``). For validating
the *integrator itself* — independent of the grid — this module also provides a direct
**Plummer-softened** pairwise force (O(N²), for tiny systems like the two-body Kepler
test and the softened central black holes):

    a_i = G Σ_{j≠i} m_j (x_j − x_i) / (|x_i − x_j|² + ε²)^{3/2}

The softening ε removes the 1/r² singularity at close approach (Plummer kernel).

This module deliberately omits ``from __future__ import annotations`` — PEP 563 would
stringize the ``ti.template()`` kernel annotations (see ``deposit.py``).
"""

import taichi as ti


@ti.kernel
def kick(
    vel_x: ti.template(),
    vel_y: ti.template(),
    vel_z: ti.template(),
    acc_x: ti.template(),
    acc_y: ti.template(),
    acc_z: ti.template(),
    n: ti.i32,
    dt: ti.f32,
):
    """Velocity update v += a·dt (call with ½dt for each leapfrog half-kick)."""
    for i in range(n):
        vel_x[i] += acc_x[i] * dt
        vel_y[i] += acc_y[i] * dt
        vel_z[i] += acc_z[i] * dt


@ti.kernel
def drift(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    vel_x: ti.template(),
    vel_y: ti.template(),
    vel_z: ti.template(),
    n: ti.i32,
    dt: ti.f32,
):
    """Position update x += v·dt (the full-step drift)."""
    for i in range(n):
        pos_x[i] += vel_x[i] * dt
        pos_y[i] += vel_y[i] * dt
        pos_z[i] += vel_z[i] * dt


@ti.kernel
def direct_accel(
    pos_x: ti.template(),
    pos_y: ti.template(),
    pos_z: ti.template(),
    mass: ti.template(),
    acc_x: ti.template(),
    acc_y: ti.template(),
    acc_z: ti.template(),
    n: ti.i32,
    grav: ti.f32,
    eps2: ti.f32,
):
    """Direct O(N²) Plummer-softened gravitational acceleration (for small N / tests)."""
    for i in range(n):
        ax = 0.0
        ay = 0.0
        az = 0.0
        xi = pos_x[i]
        yi = pos_y[i]
        zi = pos_z[i]
        for j in range(n):
            if j != i:
                rx = pos_x[j] - xi
                ry = pos_y[j] - yi
                rz = pos_z[j] - zi
                r2 = rx * rx + ry * ry + rz * rz + eps2
                inv_r3 = 1.0 / (r2 * ti.sqrt(r2))
                f = grav * mass[j] * inv_r3
                ax += f * rx
                ay += f * ry
                az += f * rz
        acc_x[i] = ax
        acc_y[i] = ay
        acc_z[i] = az


def kdk_step(parts, acc_x, acc_y, acc_z, accel_fn, dt: float) -> None:
    """Advance one KDK leapfrog step with a single acceleration evaluation per step.

    The symplectic ordering is ``kick(½dt) → drift(dt) → recompute a → kick(½dt)``. The
    trick that makes it one force eval per step (not two) is reusing the acceleration left
    over from the *previous* step for the opening half-kick:

    - **Precondition:** ``acc_{x,y,z}`` already hold a(x) for ``parts``' current positions.
      Prime this once before the loop by calling ``accel_fn()`` (which must fill the same
      ``acc_*`` fields from the current positions).
    - **Postcondition:** ``acc_{x,y,z}`` hold a(x_new) at exit — exactly what the next
      step's opening half-kick needs.

    ``accel_fn`` is the only thing that recomputes acceleration, and it is called **once**
    here (after the drift). This locks the ordering in one place so callers (``sim.py``)
    can't reintroduce a DKD/scrambled-order regression (legacy bug #1).

    ⚠️ **Contract — the cached accel must stay in sync with positions.** Because the opening
    half-kick trusts ``acc_*`` from the previous step without recomputing, *any* mutation of
    positions outside this function (teleport, boundary wrap, particle add/remove, a manual
    nudge) silently invalidates ``acc_*``: the next step would kick with a stale, wrong
    acceleration — no exception, no NaN, just a corrupted, non-symplectic trajectory. The
    rule: **after any out-of-band position change, call ``accel_fn()`` again to re-prime
    before the next ``kdk_step``.** The cleanest way to honor this is to let the main loop in
    ``sim.py`` own priming as the single point where ``acc_*`` and positions are reconciled.
    """
    n = parts.n
    half = 0.5 * dt
    kick(parts.vel_x, parts.vel_y, parts.vel_z, acc_x, acc_y, acc_z, n, half)
    drift(parts.pos_x, parts.pos_y, parts.pos_z, parts.vel_x, parts.vel_y, parts.vel_z, n, dt)
    accel_fn()
    kick(parts.vel_x, parts.vel_y, parts.vel_z, acc_x, acc_y, acc_z, n, half)


__all__ = ["kick", "drift", "direct_accel", "kdk_step"]
