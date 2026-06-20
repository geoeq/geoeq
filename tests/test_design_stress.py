"""Tests for ge.design.stress and ge.design.seepage."""

import numpy as np
import pytest

from geoeq.design import (
    total_stress, pore_pressure, effective_stress, capillary_rise,
    darcy_flow, hydraulic_gradient, critical_gradient,
    equivalent_k, flow_net,
)


# --- total_stress / pore_pressure / effective_stress -----------------

def test_total_stress_scalar():
    assert total_stress(18.0, 5.0) == pytest.approx(90.0)


def test_total_stress_layers():
    assert total_stress([17, 19, 20], [2, 3, 5]) == pytest.approx(
        17 * 2 + 19 * 3 + 20 * 5)


def test_total_stress_mismatched_lengths():
    with pytest.raises(ValueError):
        total_stress([17, 19], [2, 3, 5])


def test_pore_pressure_above_wt_is_zero():
    assert pore_pressure(2.0, z_w=5.0) == 0.0


def test_pore_pressure_below_wt():
    assert pore_pressure(5.0, z_w=2.0) == pytest.approx(9.81 * 3, rel=1e-3)


def test_pore_pressure_array():
    z = np.array([0, 2, 5, 10])
    u = pore_pressure(z, z_w=2)
    assert u[0] == 0
    assert u[1] == 0
    assert u[2] == pytest.approx(9.81 * 3, rel=1e-3)
    assert u[3] == pytest.approx(9.81 * 8, rel=1e-3)


def test_effective_stress_scalar():
    assert effective_stress(154.0, 58.86) == pytest.approx(95.14)


def test_effective_stress_array():
    res = effective_stress(np.array([100, 200]), np.array([30, 80]))
    np.testing.assert_allclose(res, [70, 120])


# --- capillary rise --------------------------------------------------

def test_capillary_rise_positive():
    # D10 = 0.05 mm, e = 0.7 -> h_c = 0.2 / (0.7 * 0.005 cm) ~ 57 cm
    h = capillary_rise(D10=0.05, e=0.7)
    assert 0.3 < h < 1.0  # metres


def test_capillary_rise_requires_positive():
    with pytest.raises(ValueError):
        capillary_rise(D10=-0.1, e=0.7)


# --- Darcy's law -----------------------------------------------------

def test_darcy_flow():
    # k=1e-4 m/s, i=0.1, A=1 m^2 -> Q = 1e-5 m^3/s
    assert darcy_flow(1e-4, 0.1, 1.0) == pytest.approx(1e-5)


def test_hydraulic_gradient():
    assert hydraulic_gradient(2.0, 10.0) == pytest.approx(0.2)


def test_critical_gradient():
    # Gs=2.65, e=0.65: i_cr = 1.65/1.65 = 1.0
    assert critical_gradient(2.65, 0.65) == pytest.approx(1.0)


def test_critical_gradient_typical_value():
    # For Gs=2.7, e=0.7: i_cr = 1.7/1.7 = 1.0
    assert critical_gradient(2.7, 0.7) == pytest.approx(1.0)


def test_critical_gradient_invalid_Gs():
    with pytest.raises(ValueError):
        critical_gradient(0.5, 0.7)


# --- Equivalent k ----------------------------------------------------

def test_equivalent_k_horizontal():
    # Two equal-thickness layers, k1=1e-3, k2=1e-5: average = 5.05e-4
    k_eq = equivalent_k([1e-3, 1e-5], [5, 5], direction="horizontal")
    assert k_eq == pytest.approx(5.05e-4)


def test_equivalent_k_vertical_lowest_dominates():
    # Vertical (series) is harmonic-like: dominated by the smaller k.
    k_eq = equivalent_k([1e-3, 1e-5], [5, 5], direction="vertical")
    # 10 / (5/1e-3 + 5/1e-5) = 10 / (5000 + 500000) ~ 1.98e-5
    assert k_eq == pytest.approx(1.98e-5, rel=0.05)


def test_equivalent_k_direction_check():
    with pytest.raises(ValueError):
        equivalent_k([1e-3], [1], direction="diagonal")


# --- Flow net --------------------------------------------------------

def test_flow_net_textbook():
    # Das example: Nf=4, Nd=8, k=1e-5 m/s, dh=10 m -> Q = 5e-5 m^3/s per m
    assert flow_net(Nf=4, Nd=8, k=1e-5, dh=10) == pytest.approx(5e-5)
