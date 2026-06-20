"""Tests for geoeq.site.field_perm — Field permeability tests."""

import numpy as np
import pytest
from geoeq.site.field_perm import (
    slug_test,
    pumping_test_confined,
    pumping_test_unconfined,
    lefranc_test,
)


class TestSlugTest:
    def test_hvorslev(self):
        # k = r²*ln(Le/R) / (2*Le*T0)
        r, R, Le, T0 = 0.025, 0.05, 1.5, 300
        expected = r**2 * np.log(Le / R) / (2 * Le * T0)
        result = slug_test(r, R, Le, T0)
        assert result == pytest.approx(expected, rel=1e-6)

    def test_order_of_magnitude(self):
        # Typical clay: k ~ 1e-8 to 1e-6 m/s
        k = slug_test(0.025, 0.05, 1.5, 3000)
        assert 1e-9 < k < 1e-5

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            slug_test(0.025, 0.05, 1.5, 300, method="bouwer_rice")


class TestPumpingTestConfined:
    def test_thiem_equation(self):
        # k = Q*ln(r2/r1) / (2*pi*H*(h2-h1))
        Q, h1, h2, r1, r2, H = 0.004, 18.0, 19.5, 10, 50, 20
        expected = Q * np.log(r2 / r1) / (2 * np.pi * H * (h2 - h1))
        res = pumping_test_confined(Q, h1, h2, r1, r2, H)
        assert res["k"] == pytest.approx(expected, rel=1e-6)
        assert res["T"] == pytest.approx(expected * H, rel=1e-6)

    def test_r2_less_than_r1_raises(self):
        with pytest.raises(ValueError):
            pumping_test_confined(0.004, 18, 19.5, 50, 10, 20)

    def test_h2_less_than_h1_raises(self):
        with pytest.raises(ValueError):
            pumping_test_confined(0.004, 19.5, 18.0, 10, 50, 20)


class TestPumpingTestUnconfined:
    def test_thiem_equation(self):
        # k = Q*ln(r2/r1) / (pi*(h2²-h1²))
        Q, h1, h2, r1, r2 = 0.003, 15.0, 18.0, 10, 50
        expected = Q * np.log(r2 / r1) / (np.pi * (h2**2 - h1**2))
        res = pumping_test_unconfined(Q, h1, h2, r1, r2)
        assert res["k"] == pytest.approx(expected, rel=1e-6)

    def test_order_of_magnitude(self):
        res = pumping_test_unconfined(0.003, 15, 18, 10, 50)
        assert 1e-6 < res["k"] < 1e-3


class TestLefrancTest:
    def test_basic(self):
        k = lefranc_test(Q=1e-5, H=2.0, D=0.076)
        assert k > 0
        assert 1e-8 < k < 1e-3

    def test_higher_flow_rate_higher_k(self):
        k1 = lefranc_test(1e-5, 2.0, 0.076)
        k2 = lefranc_test(2e-5, 2.0, 0.076)
        assert k2 > k1
