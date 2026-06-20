"""Tests for ge.design.bearing_capacity."""

import numpy as np
import pytest

from geoeq.design import (
    bearing_factors, bearing_capacity, bearing_allowable,
    bearing_shape_factors, bearing_depth_factors,
    bearing_inclination_factors,
)


# --- N-factors -- cross-check against Das Table 4.1-4.3 -------------

def test_Nc_at_phi_zero():
    f = bearing_factors(0, method="meyerhof")
    assert f["Nc"] == pytest.approx(5.14, abs=0.05)
    assert f["Nq"] == pytest.approx(1.0)
    assert f["Ngamma"] == pytest.approx(0.0)


def test_meyerhof_phi_30():
    # Das Table 4.2: phi=30 -> Nc=30.14, Nq=18.40, Ng=15.67
    f = bearing_factors(30, method="meyerhof")
    assert f["Nc"] == pytest.approx(30.14, rel=0.02)
    assert f["Nq"] == pytest.approx(18.40, rel=0.02)
    assert f["Ngamma"] == pytest.approx(15.67, rel=0.05)


def test_hansen_phi_30():
    # Das Table 4.2: hansen Ng at phi=30 = 14.39
    f = bearing_factors(30, method="hansen")
    assert f["Ngamma"] == pytest.approx(14.39, rel=0.05)


def test_vesic_phi_30():
    # Vesic Ng at phi=30 = 22.40 (Das Table 4.3)
    f = bearing_factors(30, method="vesic")
    assert f["Ngamma"] == pytest.approx(22.40, rel=0.05)


def test_method_validation():
    with pytest.raises(ValueError):
        bearing_factors(30, method="bogus")


# --- Shape / depth / inclination factors -----------------------------

def test_shape_factors_strip_neutral():
    # B/L -> 0 should give sc, sq, sg -> 1.
    s = bearing_shape_factors(B=1, L=1e6, phi=30, method="meyerhof")
    assert s["sc"] == pytest.approx(1.0, abs=0.001)
    assert s["sq"] == pytest.approx(1.0, abs=0.001)


def test_depth_factors_no_embedment():
    d = bearing_depth_factors(Df=0, B=2, phi=30, method="meyerhof")
    assert d["dc"] == pytest.approx(1.0)
    assert d["dq"] == pytest.approx(1.0)


def test_inclination_factors_no_inclination():
    i = bearing_inclination_factors(beta=0, phi=30)
    assert i["ic"] == 1.0
    assert i["iq"] == 1.0
    assert i["i_gamma"] == 1.0


# --- Full bearing-capacity workflow ----------------------------------

def test_bearing_strip_meyerhof_textbook():
    # Das Example 4.1 (close): c=0, phi=30, gamma=18, Df=1, B=2
    # q_u (strip) = q*Nq + 0.5 * gamma * B * Ng
    #             = 18*1*18.4 + 0.5*18*2*15.67
    #             = 331.2 + 282.06 ~ 613 kPa
    res = bearing_capacity(c=0, gamma=18, Df=1, B=2, phi=30,
                           L=None, method="meyerhof",
                           shape=False, depth=False)
    assert res["q_u"] == pytest.approx(613, rel=0.05)


def test_bearing_square_increases_capacity():
    # Square should have higher capacity than strip via shape factors.
    res_strip = bearing_capacity(c=10, gamma=18, Df=1, B=2, phi=30,
                                 L=None, method="meyerhof")
    res_square = bearing_capacity(c=10, gamma=18, Df=1, B=2, phi=30,
                                  L=2, method="meyerhof")
    assert res_square["q_u"] > res_strip["q_u"]


def test_bearing_phi0_clay():
    # Undrained (Skempton-like) for c=Su, phi=0, square footing.
    # q_u = 5.14 * Su * sc + q*1 (no friction contribution)
    res = bearing_capacity(c=50, gamma=18, Df=1, B=2, phi=0,
                           L=2, method="meyerhof")
    # Should be around 5.14 * 50 * 1.2 + 18 = ~325 kPa
    assert 280 < res["q_u"] < 380


def test_bearing_allowable_FS3():
    assert bearing_allowable(900, FS=3) == 300.0


def test_terzaghi_strip():
    # Das Ex. 4.3 (Terzaghi general shear): c=0, phi=30, gamma=18, Df=1.5, B=1.5
    res = bearing_capacity(c=0, gamma=18, Df=1.5, B=1.5, phi=30,
                           method="terzaghi")
    assert res["q_u"] > 0
    assert res["method"] == "terzaghi"
