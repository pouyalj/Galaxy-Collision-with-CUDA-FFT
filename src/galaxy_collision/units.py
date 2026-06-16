"""Unit system for the simulator (AGENT.md §5.5, decision §9 Q2 → resolved 2026-06-16).

**Internal units are (kpc, Myr, M_sun).** Everything in the simulation — positions,
velocities, the gravitational constant, forces — is expressed in these units. The one
documented system; no ad-hoc magic numbers (a primary lesson from the 2020 code).

The gravitational constant is *derived* here from the IAU/astronomy standard

    G = 4.30091e-3  pc (km/s)^2 / M_sun

via the conversion factors below, rather than hard-coded, so the derivation is
itself testable (see ``tests/test_units.py``). In internal units this lands at

    G ≈ 4.498e-12  kpc^3 M_sun^-1 Myr^-2.

The original code's ``dt = 1`` corresponds to 1e4 yr = 0.01 Myr in this system.
"""

from __future__ import annotations

# --- Base conversion factors ----------------------------------------------------
# Lengths. IAU 2015: 1 pc = 3.0856775814913673e16 m, exactly via the au definition.
M_PER_PC = 3.0856775814913673e16
PC_PER_KPC = 1.0e3
KPC_PER_PC = 1.0 / PC_PER_KPC
M_PER_KPC = M_PER_PC * PC_PER_KPC  # 1 kpc in metres (3.0857e19 m)
KM_PER_KPC = M_PER_KPC * 1.0e-3  # 1 kpc in km (3.0857e16 km)

# Time. Julian year = 365.25 d × 86400 s; 1 Myr = 1e6 yr.
S_PER_YR = 365.25 * 86400.0  # 3.15576e7
S_PER_MYR = S_PER_YR * 1.0e6  # 3.15576e13
YR_PER_MYR = 1.0e6

# Velocity. 1 kpc/Myr expressed in km/s.
KMS_PER_KPC_PER_MYR = KM_PER_KPC / S_PER_MYR  # ≈ 977.79 km/s

# --- The gravitational constant -------------------------------------------------
# Astronomy-standard value (Binney & Tremaine; consistent with the paper's citations).
G_STANDARD_PC_KMS = 4.30091e-3  # pc (km/s)^2 / M_sun

# Derived: convert pc → kpc and (km/s)^2 → (kpc/Myr)^2.
G = G_STANDARD_PC_KMS * KPC_PER_PC / KMS_PER_KPC_PER_MYR**2
"""Gravitational constant in internal units: kpc^3 M_sun^-1 Myr^-2 (≈ 4.498e-12)."""

# Useful reference value sometimes quoted in Gyr (G ≈ 4.498e-6 kpc^3 M_sun^-1 Gyr^-2).
# 1 Myr^-2 = (1e-3 Gyr)^-2 = 1e6 Gyr^-2, so G_[Gyr] = G_[Myr] × 1e6.
G_KPC_GYR = G * 1.0e6


def kms_to_kpc_per_myr(v_kms: float) -> float:
    """Convert a speed from km/s to internal velocity units (kpc/Myr)."""
    return v_kms / KMS_PER_KPC_PER_MYR


def kpc_per_myr_to_kms(v_internal: float) -> float:
    """Convert a speed from internal units (kpc/Myr) to km/s."""
    return v_internal * KMS_PER_KPC_PER_MYR


def myr_to_yr(t_myr: float) -> float:
    """Convert a time from Myr to years."""
    return t_myr * YR_PER_MYR


def yr_to_myr(t_yr: float) -> float:
    """Convert a time from years to Myr."""
    return t_yr / YR_PER_MYR


__all__ = [
    "G",
    "G_STANDARD_PC_KMS",
    "G_KPC_GYR",
    "M_PER_PC",
    "M_PER_KPC",
    "KM_PER_KPC",
    "PC_PER_KPC",
    "KPC_PER_PC",
    "S_PER_YR",
    "S_PER_MYR",
    "YR_PER_MYR",
    "KMS_PER_KPC_PER_MYR",
    "kms_to_kpc_per_myr",
    "kpc_per_myr_to_kms",
    "myr_to_yr",
    "yr_to_myr",
]
