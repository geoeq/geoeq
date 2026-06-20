"""Tests for geoeq.site.pile_load — Pile Load Test interpretation."""

import numpy as np
import pytest
from geoeq.site.pile_load import (
    davisson,
    chin,
    de_beer,
    hansen_80,
    fhwa_5_percent,
    case_method,
    hiley,
    danish_formula,
    enr,
    pile_impedance,
    beta_integrity,
)


class TestDavisson:
    def test_basic(self):
        Q = [0, 200, 400, 600, 800, 1000, 1200, 1400]
        s = [0, 1.0, 2.5, 5.0, 8.5, 14.0, 22.0, 35.0]
        res = davisson(Q, s, diameter=300, length=12000, area=70686,
                       elastic_modulus=30)
        assert not np.isnan(res["Qu_davisson"])
        assert res["Qu_davisson"] > 0
        assert res["settlement_at_failure"] > 0

    def test_offset_formula(self):
        res = davisson([0, 500, 1000, 1500, 2000], [0, 2, 5, 12, 30],
                       diameter=300, length=12000, area=70686, elastic_modulus=30)
        # Offset = (0.15 + D_in/120) * 25.4, D_in = 300/25.4
        D_in = 300 / 25.4
        expected = (0.15 + D_in / 120.0) * 25.4
        assert res["elastic_compression_offset"] == pytest.approx(expected, rel=1e-2)

    def test_short_data_raises(self):
        with pytest.raises(ValueError):
            davisson([0, 100], [0, 1], diameter=300, length=12000,
                     area=70686, elastic_modulus=30)


class TestChin:
    def test_basic(self):
        Q = [100, 200, 400, 600, 800, 1000, 1200]
        s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
        res = chin(Q, s)
        assert res["Qu_chin"] > 0
        assert res["r_squared"] > 0.8

    def test_overestimates(self):
        # Chin typically overestimates actual capacity
        Q = [100, 200, 400, 600, 800, 1000, 1200]
        s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
        res = chin(Q, s)
        assert res["Qu_chin"] > max(Q)


class TestDeBeer:
    def test_basic(self):
        Q = [50, 100, 200, 400, 600, 800, 1000, 1200]
        s = [0.3, 0.7, 1.5, 3.5, 7.0, 14.0, 25.0, 42.0]
        res = de_beer(Q, s)
        assert res["Qu_de_beer"] > 0

    def test_negative_values_raise(self):
        with pytest.raises(ValueError):
            de_beer([0, 100, 200], [0, 1, 2])


class TestHansen80:
    def test_basic(self):
        Q = [100, 200, 400, 600, 800, 1000, 1200]
        s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
        res = hansen_80(Q, s)
        assert res["Qu_hansen"] > 0
        assert res["settlement_at_failure"] > 0
        assert res["r_squared"] > 0.8


class TestFhwa5Percent:
    def test_basic(self):
        Q = [0, 500, 1000, 1500, 2000, 2500, 3000]
        s = [0, 2.0, 5.0, 12.0, 25.0, 45.0, 75.0]
        res = fhwa_5_percent(Q, s, diameter=600)
        assert res["settlement_criterion"] == 30.0
        assert res["Qu_fhwa"] > 0

    def test_insufficient_settlement(self):
        Q = [0, 100, 200]
        s = [0, 1, 2]
        res = fhwa_5_percent(Q, s, diameter=600)
        assert np.isnan(res["Qu_fhwa"])


class TestCaseMethod:
    def test_basic(self):
        res = case_method(F1=2500, F2=800, v1=3.0, v2=0.5,
                          impedance=4000, Jc=0.4)
        assert res["R_total"] > 0
        assert res["R_static"] > 0
        assert res["R_static"] < res["R_total"]

    def test_sand_low_Jc(self):
        res = case_method(F1=2000, F2=600, v1=2.5, v2=0.3,
                          impedance=3500, Jc=0.15)
        assert res["R_static"] > 0


class TestHiley:
    def test_basic(self):
        res = hiley(50, 1.5, 0.85, 0.010, 0.005)
        # Qu = 50*1.5*0.85 / (0.010 + 0.005/2) = 63.75 / 0.0125 = 5100
        assert res["Qu_hiley"] == pytest.approx(5100.0, rel=1e-3)
        assert res["energy_input"] == pytest.approx(63.75, rel=1e-3)

    def test_zero_set_raises(self):
        with pytest.raises(ValueError):
            hiley(50, 1.5, 0.85, 0.0, 0.005)


class TestDanishFormula:
    def test_basic(self):
        res = danish_formula(50, 1.5, 0.85, 0.010, 15, 0.07, 25e6)
        assert res["Qu_danish"] > 0
        assert res["elastic_compression_term"] > 0


class TestEnr:
    def test_basic(self):
        res = enr(30, 1.0, 0.010)
        # Qa = 30*1.0 / (6 * (0.010 + 0.0254)) = 30 / 0.2124 = 141.2
        Qa = 30 / (6 * (0.010 + 0.0254))
        assert res["Qa_enr"] == pytest.approx(Qa, rel=1e-2)
        assert res["Qu_enr"] == pytest.approx(Qa * 6, rel=1e-2)


class TestPileImpedance:
    def test_basic(self):
        # Z = 30e6 * 0.07 / 3800 = 552.6
        res = pile_impedance(30e6, 0.07, 3800)
        assert res["impedance"] == pytest.approx(552.6, rel=1e-2)


class TestBetaIntegrity:
    def test_intact(self):
        res = beta_integrity(14.8, 15.0)
        assert 0.95 <= res["beta"] <= 1.05
        assert "Intact" in res["interpretation"]

    def test_defect(self):
        res = beta_integrity(12.0, 15.0)
        assert res["beta"] < 0.95
        assert "defect" in res["interpretation"].lower()
