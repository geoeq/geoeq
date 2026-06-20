"""Tests for ge.design.settlement."""

import numpy as np
import pytest

from geoeq.design import (
    settlement_immediate, settlement_primary, settlement_secondary,
    settlement_schmertmann, time_factor, consolidation_degree,
    consolidation_time,
)


# --- immediate / elastic --------------------------------------------

def test_immediate_simple_square():
    # q=200 kPa, B=2 m, Es=10 MPa = 10000 kPa, mu=0.3, square rigid
    # S = 200 * 2 * (1 - 0.09) / 10000 * 0.88 = 0.032 m = 32 mm
    s = settlement_immediate(q=200, B=2, Es=10000, mu=0.3,
                             shape="rigid_square")
    assert s == pytest.approx(0.032, rel=0.05)


def test_immediate_known_Ip():
    s = settlement_immediate(q=100, B=1, Es=50000, mu=0.3, Ip=1.0)
    expected = 100 * 1 * (1 - 0.09) / 50000 * 1.0
    assert s == pytest.approx(expected)


def test_immediate_bad_shape_raises():
    with pytest.raises(ValueError):
        settlement_immediate(100, 1, 50000, shape="weird")


# --- primary consolidation ------------------------------------------

def test_primary_NC():
    # Das Ex 5.1: Cc=0.27, e0=0.92, H=3.66, sigma0=80, dsigma=15
    # S = (0.27 * 3.66 / 1.92) * log10(95/80) = 0.515 * 0.0744 = 0.0383 m
    s = settlement_primary(Cc=0.27, e0=0.92, H=3.66,
                            sigma0=80, delta_sigma=15)
    assert s == pytest.approx(0.0383, rel=0.05)


def test_primary_OC_below_pc():
    # Final stress < pc -> uses Cr (smaller settlement).
    s = settlement_primary(Cc=0.27, e0=0.92, H=3, sigma0=80,
                            delta_sigma=10, Cr=0.05, sigma_pc=100)
    # Should be much smaller than NC case.
    assert 0 < s < 0.01


def test_primary_OC_crossing_pc():
    # Final stress crosses pc -> two-segment formula.
    s = settlement_primary(Cc=0.27, e0=0.92, H=3, sigma0=80,
                            delta_sigma=80, Cr=0.05, sigma_pc=100)
    # Should be between Cr-only and Cc-only solutions.
    assert s > 0.01


def test_primary_OC_requires_cr_and_pc():
    with pytest.raises(ValueError):
        settlement_primary(Cc=0.27, e0=0.92, H=3, sigma0=80,
                           delta_sigma=10, kind="OC")


# --- secondary compression ------------------------------------------

def test_secondary_log_time():
    # C_alpha = 0.005 (strain form), H=2, t1=100 d, t2=1000 d
    # S = 0.005 * 2 * log10(10) = 0.01 m
    s = settlement_secondary(C_alpha=0.005, H=2, t1=100, t2=1000)
    assert s == pytest.approx(0.01)


def test_secondary_void_ratio_form():
    # C_alpha = 0.02, e_p = 1.0 -> C_alpha_eps = 0.01
    s = settlement_secondary(C_alpha=0.02, H=2, t1=100, t2=1000, e_p=1.0)
    assert s == pytest.approx(0.02)


def test_secondary_t2_must_exceed_t1():
    with pytest.raises(ValueError):
        settlement_secondary(0.005, 2, t1=1000, t2=100)


# --- Schmertmann ---------------------------------------------------

def test_schmertmann_square():
    # q_net=100, B=2, Es=20000, Df=1, gamma=18 -> S ~ a few mm.
    res = settlement_schmertmann(
        q_net=100, B=2, Es=20000, Df=1, shape="square")
    assert res["S"] > 0
    assert res["Iz_peak"] > 0.5
    assert res["z_max"] == 4.0  # 2B


def test_schmertmann_strip_deeper_influence():
    res = settlement_schmertmann(
        q_net=100, B=2, Es=20000, Df=1, shape="strip")
    assert res["z_max"] == 8.0  # 4B for strip


def test_schmertmann_creep_amplifies():
    res1 = settlement_schmertmann(q_net=100, B=2, Es=20000, Df=1, t_years=1)
    res10 = settlement_schmertmann(q_net=100, B=2, Es=20000, Df=1, t_years=10)
    assert res10["S"] > res1["S"]


def test_schmertmann_multilayer():
    res = settlement_schmertmann(
        q_net=100, B=2, Es=[10000, 30000],
        layers=[(0, 2), (2, 4)], Df=1, shape="square")
    assert "S" in res
    assert len(res["layers"]) == 2


# --- time factor / U% ---------------------------------------------

def test_time_factor_U_small():
    # Tv at U=0.5 = pi/4 * 0.25 = 0.196
    assert time_factor(0.5) == pytest.approx(0.196, rel=0.01)


def test_time_factor_U_large():
    # Tv at U=0.9: 1.781 - 0.933 * log10(10) = 0.848
    assert time_factor(0.9) == pytest.approx(0.848, rel=0.01)


def test_consolidation_degree_inverse():
    # Should be close to identity round-trip.
    for U in (0.3, 0.5, 0.7, 0.9, 0.95):
        Tv = time_factor(U)
        U_recovered = consolidation_degree(Tv)
        assert U_recovered == pytest.approx(U, abs=0.005)


def test_consolidation_time():
    # U=0.5, Hdr=1 m, cv=1e-7 m^2/s -> t = 0.196 / 1e-7 = 1.96e6 s
    t = consolidation_time(U=0.5, Hdr=1, cv=1e-7)
    assert t == pytest.approx(1.96e6, rel=0.01)


def test_time_factor_U_out_of_range():
    with pytest.raises(ValueError):
        time_factor(-0.1)
    with pytest.raises(ValueError):
        time_factor(1.0)
