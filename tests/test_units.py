"""Tests for the (kpc, Myr, M_sun) unit system (Stage 1).

The headline check: the derived gravitational constant matches the astronomy
standard once converted into internal units.
"""

from __future__ import annotations

import math

from galaxy_collision import units


def test_G_matches_expected_internal_value():
    # AGENT.md §5.5: G ≈ 4.498e-12 kpc^3 M_sun^-1 Myr^-2.
    assert math.isclose(units.G, 4.498e-12, rel_tol=1e-3)


def test_G_round_trips_to_standard():
    # Convert internal G back to pc (km/s)^2 / M_sun and recover the input constant.
    back = units.G / units.KPC_PER_PC * units.KMS_PER_KPC_PER_MYR**2
    assert math.isclose(back, units.G_STANDARD_PC_KMS, rel_tol=1e-12)


def test_G_gyr_reference():
    # G in (kpc, Gyr, M_sun) ≈ 4.498e-6.
    assert math.isclose(units.G_KPC_GYR, 4.498e-6, rel_tol=1e-3)


def test_kpc_per_myr_in_kms():
    # 1 kpc/Myr ≈ 977.79 km/s (a standard textbook figure).
    assert math.isclose(units.KMS_PER_KPC_PER_MYR, 977.79, rel_tol=1e-3)


def test_velocity_conversions_round_trip():
    for v_kms in (0.0, 1.0, 125.0, 402_000.0 / 3600.0):  # incl. the paper's ~v
        internal = units.kms_to_kpc_per_myr(v_kms)
        assert math.isclose(units.kpc_per_myr_to_kms(internal), v_kms, rel_tol=1e-12)


def test_paper_approach_speed_sanity():
    # Paper: v ≈ 125,000 mph ≈ 402,000 km/h ≈ 111.7 km/s → a small internal speed.
    v_kms = 402_000.0 / 3600.0
    v_internal = units.kms_to_kpc_per_myr(v_kms)
    # ~0.11 kpc/Myr: well below 1, sanity-checks the conversion scale.
    assert 0.05 < v_internal < 0.2


def test_time_conversions():
    assert math.isclose(units.myr_to_yr(0.01), 1e4)  # the 2020 code's dt=1 step
    assert math.isclose(units.yr_to_myr(1e6), 1.0)


def test_length_factors_consistent():
    # 1 kpc = 1000 pc, expressed in metres.
    assert math.isclose(units.M_PER_KPC, units.M_PER_PC * 1000.0, rel_tol=1e-12)
