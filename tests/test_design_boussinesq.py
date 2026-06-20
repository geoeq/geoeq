"""Tests for ge.design.boussinesq -- stress distribution."""

import numpy as np
import pytest

from geoeq.design import (
    boussinesq_point, boussinesq_line, boussinesq_strip,
    boussinesq_circular, boussinesq_rect, westergaard_point,
    newmark_influence, stress_2to1, stress_bulb,
)


# --- point load ------------------------------------------------------

def test_point_directly_below():
    # Below P=100 at z=2: dsigma = 3*100/(2*pi*4) * 1 = 11.94 kPa
    expected = 3 * 100 / (2 * np.pi * 4)
    assert boussinesq_point(100, 2, 0) == pytest.approx(expected)


def test_point_off_axis_smaller():
    # Off-axis stress always less than on-axis at same depth.
    assert boussinesq_point(100, 2, 3) < boussinesq_point(100, 2, 0)


def test_point_depth_decay():
    # Stress decays as ~1/z^2 directly below.
    s1 = boussinesq_point(100, 1)
    s2 = boussinesq_point(100, 2)
    assert s2 == pytest.approx(s1 / 4, rel=0.01)


# --- line load -------------------------------------------------------

def test_line_load_on_axis():
    # 2q/(pi z) at x=0
    expected = 2 * 50 / (np.pi * 3)
    assert boussinesq_line(50, 3, 0) == pytest.approx(expected)


# --- strip load ------------------------------------------------------

def test_strip_at_great_depth_approaches_zero():
    # Stress decays at depth.
    s_shallow = boussinesq_strip(100, B=2, z=1)
    s_deep = boussinesq_strip(100, B=2, z=20)
    assert s_deep < s_shallow / 5


def test_strip_at_surface_equals_q():
    # Right at the surface centred under the strip, sigma -> q (within model).
    s = boussinesq_strip(100, B=10, z=0.01)
    assert s > 90  # close to q (=100)


# --- circular load ---------------------------------------------------

def test_circular_textbook_das_table():
    # Das Table 6.6: q=100, z/R = 1 -> dsigma/q = 0.6465.
    # Centre below a circle of radius R=2 at z=2.
    s = boussinesq_circular(q=100, R=2, z=2)
    assert s == pytest.approx(64.65, rel=0.01)


def test_circular_at_surface():
    # z very small -> stress approaches q.
    s = boussinesq_circular(q=100, R=5, z=0.001)
    assert s == pytest.approx(100.0, abs=0.01)


# --- rectangular / Fadum ---------------------------------------------

def test_fadum_textbook_value():
    # Das Table 6.5: m=1, n=1 -> I = 0.1752.
    I = newmark_influence(m=1.0, n=1.0)
    assert I == pytest.approx(0.1752, rel=0.005)


def test_fadum_m05_n05():
    # Das Table 6.5: m=0.5, n=0.5 -> I = 0.0840
    I = newmark_influence(m=0.5, n=0.5)
    assert I == pytest.approx(0.0840, rel=0.01)


def test_rect_corner_vs_centre():
    # Centre stress = 4 * corner stress for sub-rectangles of half size.
    q, B, L, z = 100, 4, 6, 2
    s_corner = boussinesq_rect(q, B, L, z, position="corner")
    s_centre = boussinesq_rect(q, B, L, z, position="centre")
    # Centre will be > corner.
    assert s_centre > s_corner


# --- Westergaard -----------------------------------------------------

def test_westergaard_textbook_value_mu0():
    # Per Das (2010) Eq. 6.13, mu=0, r=0: Δσz = 2P/(π·z²).
    P, z = 100, 2
    s = westergaard_point(P, z, r=0, mu=0.0)
    expected = 2 * P / (np.pi * z ** 2)
    assert s == pytest.approx(expected, rel=0.005)


def test_westergaard_decays_with_r():
    # Stress decreases as r increases.
    P, z = 100, 2
    s_on = westergaard_point(P, z, r=0, mu=0.0)
    s_off = westergaard_point(P, z, r=3, mu=0.0)
    assert s_off < s_on


def test_westergaard_positive():
    # Sanity: positive, decays with depth.
    s_near = westergaard_point(100, 1, r=0, mu=0.0)
    s_far = westergaard_point(100, 5, r=0, mu=0.0)
    assert s_far < s_near


# --- 2:1 approximation ----------------------------------------------

def test_2to1_at_surface():
    # z=0 -> stress = q.
    assert stress_2to1(100, B=2, L=4, z=1e-6) == pytest.approx(100, abs=0.01)


def test_2to1_decay():
    s = stress_2to1(100, B=2, L=4, z=5)
    # = 100 * 2 * 4 / (7 * 9) = 800/63 = 12.7
    assert s == pytest.approx(12.7, rel=0.01)


# --- bulb helper ----------------------------------------------------

def test_stress_bulb_shape():
    res = stress_bulb(P=100, n=10)
    assert res["sigma"].shape == (10, 10)
    assert np.all(res["sigma"] > 0)
