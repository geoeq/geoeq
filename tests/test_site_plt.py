"""Tests for geoeq.site.plt_test — Plate Load Test."""

import numpy as np
import pytest
from geoeq.site.plt_test import (
    plt_bearing,
    plt_subgrade_modulus,
    plt_settlement_correction,
    plt_elastic_modulus,
    plt_plot,
)


class TestPltBearing:
    def test_settlement_criterion(self):
        p = [0, 50, 100, 150, 200, 250, 300, 350, 400]
        s = [0, 1, 3, 6, 12, 20, 32, 48, 65]
        res = plt_bearing(p, s, failure_criterion="settlement")
        assert res["settlement_at_failure"] == 30.0
        assert 200 < res["qu"] < 400

    def test_log_log_criterion(self):
        p = [50, 100, 150, 200, 250, 300, 350, 400]
        s = [0.5, 1.2, 2.5, 5.0, 10.0, 20.0, 35.0, 55.0]
        res = plt_bearing(p, s, failure_criterion="log_log")
        assert res["qu"] > 0

    def test_tangent_criterion(self):
        p = [0, 50, 100, 150, 200, 250, 300, 350, 400]
        s = [0, 0.5, 1.2, 2.5, 5.0, 10.0, 20.0, 35.0, 55.0]
        res = plt_bearing(p, s, failure_criterion="tangent")
        assert res["qu"] > 0

    def test_invalid_criterion(self):
        with pytest.raises(ValueError):
            plt_bearing([0, 100], [0, 5], failure_criterion="bogus")

    def test_insufficient_data(self):
        with pytest.raises(ValueError):
            plt_bearing([100, 200], [1, 2])


class TestPltSubgradeModulus:
    def test_scalar(self):
        # ks = 200 / 0.005 = 40000 kN/m³
        assert plt_subgrade_modulus(200, 5.0) == pytest.approx(40000.0, rel=1e-6)

    def test_array(self):
        p = np.array([100, 200, 300])
        s = np.array([2.5, 5.0, 7.5])
        res = plt_subgrade_modulus(p, s)
        assert "mean_ks" in res
        assert res["mean_ks"] == pytest.approx(40000.0, rel=1e-6)


class TestPltSettlementCorrection:
    def test_sand_terzaghi_peck(self):
        # s_f/s_p = (Bf*(Bp+0.3) / (Bp*(Bf+0.3)))²
        Bp, Bf = 0.3, 2.0
        ratio = (Bf * (Bp + 0.3) / (Bp * (Bf + 0.3)))**2
        expected = 5.0 * ratio
        result = plt_settlement_correction(5.0, 0.3, 2.0, "sand")
        assert result == pytest.approx(expected, rel=1e-3)

    def test_clay_linear(self):
        # s_f/s_p = Bf/Bp = 2.0/0.3
        expected = 5.0 * (2.0 / 0.3)
        result = plt_settlement_correction(5.0, 0.3, 2.0, "clay")
        assert result == pytest.approx(expected, rel=1e-3)

    def test_same_size(self):
        result = plt_settlement_correction(5.0, 0.3, 0.3, "sand")
        assert result == pytest.approx(5.0, rel=1e-3)

    def test_invalid_soil(self):
        with pytest.raises(ValueError):
            plt_settlement_correction(5.0, 0.3, 2.0, "rock")


class TestPltElasticModulus:
    def test_basic(self):
        # E = q*B*(1-nu²)*Is / delta
        # E = 200 * 0.3 * (1 - 0.09) * 0.79 / 0.005
        E = 200 * 0.3 * 0.91 * 0.79 / 0.005
        result = plt_elastic_modulus(200, 5.0)
        assert result == pytest.approx(E, rel=1e-3)

    def test_increases_with_pressure(self):
        E1 = plt_elastic_modulus(100, 5.0)
        E2 = plt_elastic_modulus(200, 5.0)
        assert E2 > E1


class TestPltPlot:
    def test_smoke(self):
        import matplotlib
        matplotlib.use("Agg")
        res = plt_plot([0, 100, 200, 300], [0, 2, 5, 12])
        assert "ax" in res
