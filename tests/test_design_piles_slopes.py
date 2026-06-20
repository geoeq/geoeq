"""Tests for ge.design.piles and ge.design.slopes."""

import numpy as np
import pytest

from geoeq.design import (
    pile_end_bearing, pile_skin_friction, pile_capacity,
    pile_settlement, pile_group_efficiency,
    infinite_slope, culmann, taylor_stability, bishop,
)


# --- pile end bearing ---------------------------------------------

def test_end_bearing_clay():
    res = pile_end_bearing(Su=100)
    assert res["q_p"] == pytest.approx(900.0)
    assert res["method"] == "skempton_clay"


def test_end_bearing_sand_meyerhof():
    # phi=30, sigma'_v=100 -> N_q* ~ 18.4 -> q_p ~ 1840
    res = pile_end_bearing(phi=30, sigma_v_eff=100, method="meyerhof")
    assert 1500 < res["q_p"] < 2500


def test_end_bearing_validation():
    with pytest.raises(ValueError):
        pile_end_bearing(phi=30, sigma_v_eff=100, method="bogus")


# --- skin friction ------------------------------------------------

def test_alpha_method_soft_clay():
    # Su=20 (very soft) -> alpha=1.0 -> f_s = 20
    res = pile_skin_friction(Su=20, method="alpha")
    assert res["f_s"] == pytest.approx(20.0)
    assert res["alpha"] == pytest.approx(1.0)


def test_alpha_method_stiff_clay():
    # Su=200 (stiff) -> alpha ~ 0.35
    res = pile_skin_friction(Su=200, method="alpha")
    assert res["alpha"] < 0.5
    assert res["f_s"] < 100


def test_beta_method_sand():
    res = pile_skin_friction(sigma_v_eff=100, method="beta", phi=30)
    # K=1-sin(30)=0.5, delta=30, beta = 0.5 * tan(30) = 0.289
    assert res["beta"] == pytest.approx(0.289, rel=0.02)
    assert res["f_s"] == pytest.approx(28.9, rel=0.02)


def test_lambda_method():
    res = pile_skin_friction(Su=50, sigma_v_eff=100, method="lambda")
    # f_s = 0.3 * (100 + 100) = 60
    assert res["f_s"] == pytest.approx(60.0)


# --- total capacity -----------------------------------------------

def test_pile_capacity_total():
    # D=0.3, L=10, q_p=2000, f_s=30
    res = pile_capacity(D=0.3, L=10, q_p=2000, f_s=30)
    # A_base = pi * 0.09 / 4 = 0.0707; perim = 0.942
    # Q_p = 141.4; Q_s = 30*0.942*10 = 282.7; Q_ult = 424
    assert res["Q_p"] == pytest.approx(141.4, rel=0.02)
    assert res["Q_s"] == pytest.approx(282.7, rel=0.02)
    assert res["Q_ult"] == pytest.approx(424.1, rel=0.02)
    assert res["Q_all"] == pytest.approx(141.4, rel=0.02)


def test_pile_capacity_custom_area():
    res = pile_capacity(D=0.3, L=10, q_p=2000, f_s=30,
                        area_base=0.1, perimeter=1.0)
    assert res["Q_p"] == 200.0
    assert res["Q_s"] == 300.0


# --- group efficiency --------------------------------------------

def test_group_efficiency_converse_labarre():
    # 3x3 group, D=0.3, s=0.9: theta = atan(1/3) ~ 18.43 deg
    # eta = 1 - 18.43 * (6+6) / (90*9) = 1 - 221.2/810 = 0.727
    res = pile_group_efficiency(3, 3, D=0.3, s=0.9)
    assert res["eta"] == pytest.approx(0.727, rel=0.05)


def test_group_efficiency_wide_spacing():
    # Very wide spacing -> eta close to 1.
    res = pile_group_efficiency(3, 3, D=0.3, s=3.0)
    assert res["eta"] > 0.9


def test_group_efficiency_feld():
    res = pile_group_efficiency(3, 3, D=0.3, s=1, method="feld")
    assert 0 < res["eta"] < 1


# --- pile settlement ---------------------------------------------

def test_pile_settlement_components():
    res = pile_settlement(
        Q_w=400, Q_p=200, Q_s=200,
        D=0.3, L=10, Es=20000, Ep=25e6)
    assert res["s_total"] > 0
    for k in ("s1_elastic", "s2_tip", "s3_shaft"):
        assert res[k] >= 0


# --- infinite slope ----------------------------------------------

def test_infinite_slope_dry_cohesionless():
    # FS = tan(30)/tan(20) = 0.577/0.364 = 1.587
    res = infinite_slope(phi=30, beta=20)
    assert res["FS"] == pytest.approx(1.587, rel=0.01)


def test_infinite_slope_with_seepage_reduces_FS():
    fs_dry = infinite_slope(phi=30, beta=20)["FS"]
    fs_seep = infinite_slope(phi=30, beta=20,
                              gamma=20, gamma_sat=20,
                              seepage=True)["FS"]
    assert fs_seep < fs_dry


def test_infinite_slope_cohesion_helps():
    fs0 = infinite_slope(phi=30, beta=20, c=0, H=3)["FS"]
    fs1 = infinite_slope(phi=30, beta=20, c=10, H=3, gamma=18)["FS"]
    assert fs1 > fs0


# --- Culmann -----------------------------------------------------

def test_culmann_basic():
    res = culmann(c=20, phi=20, gamma=18, beta=60)
    assert res["H_cr"] > 0


def test_culmann_requires_beta_gt_phi():
    with pytest.raises(ValueError):
        culmann(c=20, phi=30, gamma=18, beta=20)


# --- Taylor ------------------------------------------------------

def test_taylor_returns_positive_FS():
    res = taylor_stability(phi=0, c=25, gamma=18, H=5, beta=45)
    assert res["FS"] > 0
    assert 0 < res["m"] < 0.4


# --- Bishop ------------------------------------------------------

def _make_simple_slope_slices():
    """Toy 5-slice slope, all same gamma=18, phi=20, c=10."""
    slices = []
    for alpha in (-10, 5, 20, 35, 50):
        slices.append({
            "b": 2.0, "h": 5.0, "alpha": alpha,
            "c": 10, "phi": 20, "gamma": 18, "u": 0,
        })
    return slices


def test_bishop_converges():
    slices = _make_simple_slope_slices()
    res = bishop(slices)
    assert res["converged"]
    assert res["FS"] > 0


def test_bishop_high_cohesion_stable():
    slices = _make_simple_slope_slices()
    for s in slices:
        s["c"] = 200  # very high
    res = bishop(slices)
    assert res["FS"] > 2


def test_bishop_empty_slices_raises():
    with pytest.raises(ValueError):
        bishop([])
