"""Tests for ge.design.earth_pressure and ge.design.walls."""

import numpy as np
import pytest

from geoeq.design import (
    K0, Ka, Kp, earth_pressure, tension_crack_depth,
    wall_overturning, wall_sliding, wall_bearing, sheet_pile,
)


# --- K-factors -----------------------------------------------------

def test_K0_jaky_phi_30():
    # 1 - sin(30) = 0.5
    assert K0(30) == pytest.approx(0.5, abs=0.01)


def test_K0_mayne_OC():
    # OCR=1 should equal Jaky.
    assert K0(30, OCR=1, method="mayne") == pytest.approx(K0(30, method="jaky"))
    # OCR=4 should be larger.
    assert K0(30, OCR=4, method="mayne") > K0(30)


def test_Ka_rankine_phi_30():
    # tan^2(45 - 15) = tan^2(30) = 1/3
    assert Ka(30) == pytest.approx(1 / 3, abs=0.001)


def test_Kp_rankine_phi_30():
    # tan^2(45+15) = tan^2(60) = 3
    assert Kp(30) == pytest.approx(3.0, abs=0.001)


def test_Ka_Kp_reciprocal_smooth():
    # Smooth vertical walls: Ka * Kp = 1.
    assert Ka(30) * Kp(30) == pytest.approx(1.0, rel=0.001)


def test_Ka_coulomb_smooth_matches_rankine():
    # When delta=alpha=beta=0, Coulomb = Rankine.
    assert Ka(30, delta=0, alpha=0, beta=0, method="coulomb") == \
        pytest.approx(Ka(30, method="rankine"), abs=0.001)


def test_Ka_sloping_backfill():
    # Beta > 0 -> Ka > horizontal-backfill case.
    Ka_flat = Ka(30, beta=0)
    Ka_slope = Ka(30, beta=10)
    assert Ka_slope > Ka_flat


# --- distribution / resultant -------------------------------------

def test_earth_pressure_active_resultant():
    # gamma=18, H=5, phi=30: Pa = 0.5 * Ka * gamma * H^2 = 0.5 * 1/3 * 18 * 25 = 75
    res = earth_pressure(gamma=18, H=5, phi=30, kind="active")
    assert res["P_total"] == pytest.approx(75.0, rel=0.02)
    assert res["K"] == pytest.approx(1 / 3, abs=0.001)


def test_earth_pressure_passive_larger():
    a = earth_pressure(gamma=18, H=5, phi=30, kind="active")["P_total"]
    p = earth_pressure(gamma=18, H=5, phi=30, kind="passive")["P_total"]
    assert p > 5 * a


def test_earth_pressure_with_water():
    # Add water -- total should include hydrostatic.
    res = earth_pressure(gamma=18, H=5, phi=30, kind="active",
                          water_table=2)
    # P_water = 0.5 * 9.81 * 9 = 44.1 kN/m
    assert res["P_water"] == pytest.approx(44.1, rel=0.02)


def test_tension_crack_depth():
    # z_c = 2c / (gamma sqrt(Ka)); c=20, gamma=18, phi=30 -> Ka=1/3
    # z_c = 2*20 / (18 * sqrt(0.333)) = 40 / 10.39 = 3.85 m
    z_c = tension_crack_depth(c=20, gamma=18, phi=30)
    assert z_c == pytest.approx(3.85, rel=0.02)


# --- wall stability ------------------------------------------------

def test_overturning_basic():
    res = wall_overturning(resisting_moments=[300, 150], driving_moments=[100])
    assert res["FS"] == 4.5
    assert res["M_resisting"] == 450


def test_overturning_zero_driving_raises():
    with pytest.raises(ValueError):
        wall_overturning([100], [0])


def test_sliding_friction_only():
    # V=200, mu=0.6, H=80 -> FS = 200*0.6/80 = 1.5
    res = wall_sliding(horizontal_forces=[80],
                       vertical_forces=[200], mu=0.6)
    assert res["FS"] == pytest.approx(1.5)


def test_sliding_with_passive():
    res = wall_sliding(horizontal_forces=[80], vertical_forces=[200],
                       mu=0.6, Pp=20)
    assert res["FS"] > 1.5


def test_bearing_within_kern():
    # B=4, V=400, M=200 -> e = 0.5 < B/6 = 0.667
    res = wall_bearing(V=400, M_net=200, B=4)
    assert res["within_kern"]
    assert res["q_min"] > 0
    assert res["q_max"] > res["q_min"]


def test_bearing_outside_kern():
    # B=2, V=400, M=600 -> e = 1.5 > B/6 = 0.333
    res = wall_bearing(V=400, M_net=600, B=2)
    assert not res["within_kern"]
    assert res["q_min"] == 0


def test_bearing_FS():
    res = wall_bearing(V=400, M_net=100, B=4, q_ult=600)
    assert "FS" in res
    assert res["FS"] > 0


# --- sheet pile ----------------------------------------------------

def test_sheet_pile_basic():
    res = sheet_pile(gamma=18, H=4, phi=30, kind="cantilever")
    assert res["D_theory"] > 0
    assert res["D_design"] > res["D_theory"]
    assert res["Ka"] == pytest.approx(1 / 3, abs=0.001)
    assert res["Kp"] == pytest.approx(3.0, abs=0.001)
