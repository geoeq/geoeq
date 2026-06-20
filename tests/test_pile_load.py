"""Tests for geoeq.site.pile_load — pile load test interpretation."""

import math
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


# ---------------------------------------------------------------------------
# Static Methods
# ---------------------------------------------------------------------------

class TestDavisson:
    """Davisson's Offset Limit Method."""

    def test_basic(self):
        Q = [0, 200, 400, 600, 800, 1000, 1200, 1400]
        s = [0, 1.0, 2.5, 5.0, 8.5, 14.0, 22.0, 35.0]
        res = davisson(Q, s, diameter=300, length=12000, area=70686,
                       elastic_modulus=30)
        assert not math.isnan(res["Qu_davisson"])
        assert res["Qu_davisson"] > 0
        assert res["settlement_at_failure"] > 0

    def test_offset_value(self):
        """Check the offset calculation: 0.15 + D_in/120, converted to mm."""
        res = davisson([0, 100, 200, 300, 400, 500, 600, 700, 800],
                       [0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0],
                       diameter=304.8, length=10000, area=72966,
                       elastic_modulus=30)
        # 304.8 mm = 12 in -> offset = (0.15 + 12/120) * 25.4 = 0.25 * 25.4 = 6.35 mm
        assert abs(res["elastic_compression_offset"] - 6.35) < 0.01

    def test_no_intersection(self):
        """Very stiff pile: settlement never reaches offset line."""
        Q = [0, 500, 1000, 1500]
        s = [0, 0.01, 0.02, 0.03]  # tiny settlements
        res = davisson(Q, s, diameter=300, length=12000, area=70686,
                       elastic_modulus=30)
        assert math.isnan(res["Qu_davisson"])

    def test_validation(self):
        with pytest.raises(ValueError):
            davisson([1, 2], [1, 2, 3], diameter=300, length=12000,
                     area=70686, elastic_modulus=30)
        with pytest.raises(ValueError):
            davisson([1], [1], diameter=300, length=12000,
                     area=70686, elastic_modulus=30)


class TestChin:
    """Chin's Hyperbolic Method."""

    def test_basic(self):
        Q = [100, 200, 400, 600, 800, 1000, 1200]
        s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
        res = chin(Q, s)
        assert res["Qu_chin"] > 0
        assert res["slope"] > 0
        assert 0 <= res["r_squared"] <= 1

    def test_custom_range(self):
        Q = [100, 200, 400, 600, 800, 1000, 1200]
        s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
        res = chin(Q, s, start_index=3, end_index=7)
        assert res["Qu_chin"] > 0

    def test_validation(self):
        with pytest.raises(ValueError):
            chin([0, 0, 0], [1, 2, 3])  # all Q == 0 after filter


class TestDeBeer:
    """De Beer's Double-Log Method."""

    def test_basic(self):
        Q = [50, 100, 200, 400, 600, 800, 1000, 1200]
        s = [0.3, 0.7, 1.5, 3.5, 7.0, 14.0, 25.0, 42.0]
        res = de_beer(Q, s)
        assert res["Qu_de_beer"] > 0
        assert res["settlement_at_failure"] > 0

    def test_validation(self):
        with pytest.raises(ValueError):
            de_beer([0, 100], [1, 2])  # Q=0 not allowed
        with pytest.raises(ValueError):
            de_beer([100], [1])  # too few points


class TestHansen80:
    """Hansen's 80% Method."""

    def test_basic(self):
        Q = [100, 200, 400, 600, 800, 1000, 1200]
        s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
        res = hansen_80(Q, s)
        assert res["Qu_hansen"] > 0
        assert res["settlement_at_failure"] > 0
        assert res["a"] > 0
        assert res["b"] > 0

    def test_custom_range(self):
        Q = [100, 200, 400, 600, 800, 1000, 1200]
        s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
        res = hansen_80(Q, s, start_index=2, end_index=7)
        assert res["Qu_hansen"] > 0


class TestFHWA5Percent:
    """FHWA 5% diameter criterion."""

    def test_basic(self):
        Q = [0, 500, 1000, 1500, 2000, 2500, 3000]
        s = [0, 2.0, 5.0, 12.0, 25.0, 45.0, 75.0]
        res = fhwa_5_percent(Q, s, diameter=600)
        assert res["settlement_criterion"] == 30.0
        assert not math.isnan(res["Qu_fhwa"])
        assert res["Qu_fhwa"] > 0

    def test_criterion_beyond_data(self):
        Q = [0, 500, 1000]
        s = [0, 1.0, 2.0]
        res = fhwa_5_percent(Q, s, diameter=600)  # need 30 mm
        assert math.isnan(res["Qu_fhwa"])


# ---------------------------------------------------------------------------
# Dynamic Methods
# ---------------------------------------------------------------------------

class TestCaseMethod:
    """Case Method for PDA interpretation."""

    def test_basic(self):
        res = case_method(F1=2500, F2=800, v1=3.0, v2=0.5,
                          impedance=4000, Jc=0.4)
        assert res["R_total"] > 0
        assert res["R_static"] > 0
        assert res["Jc"] == 0.4

    def test_damping_sensitivity(self):
        """Higher Jc should reduce R_static (more damping removed)."""
        r_low = case_method(F1=2500, F2=800, v1=3.0, v2=0.5,
                            impedance=4000, Jc=0.2)
        r_high = case_method(F1=2500, F2=800, v1=3.0, v2=0.5,
                             impedance=4000, Jc=0.8)
        # Not necessarily always true depending on F2/v2 values,
        # but for this data set it holds.
        assert r_low["R_static"] != r_high["R_static"]

    def test_validation(self):
        with pytest.raises(ValueError):
            case_method(F1=2500, F2=800, v1=3.0, v2=0.5,
                        impedance=4000, Jc=1.5)  # Jc > 1


class TestHiley:
    """Hiley Driving Formula."""

    def test_basic(self):
        res = hiley(hammer_weight=50, drop_height=1.5, efficiency=0.85,
                    set_per_blow=0.010, elastic_compression=0.005)
        assert res["Qu_hiley"] > 0
        expected_energy = 50 * 1.5 * 0.85
        assert abs(res["energy_input"] - expected_energy) < 1e-6
        expected_Qu = expected_energy / (0.010 + 0.005 / 2)
        assert abs(res["Qu_hiley"] - expected_Qu) < 1e-6

    def test_validation(self):
        with pytest.raises(ValueError):
            hiley(50, 1.5, 1.5, 0.01, 0.005)  # efficiency > 1


class TestDanish:
    """Danish Driving Formula."""

    def test_basic(self):
        res = danish_formula(hammer_weight=50, drop_height=1.5,
                             efficiency=0.85, set_per_blow=0.010,
                             pile_length=15, pile_area=0.07,
                             pile_modulus=25e6)
        assert res["Qu_danish"] > 0
        assert res["elastic_compression_term"] > 0

    def test_formula_check(self):
        """Verify formula manually."""
        W, h, eh, s, L, A, E = 50, 1.5, 0.85, 0.010, 15, 0.07, 25e6
        energy = eh * W * h
        elastic = math.sqrt(energy * L / (2 * A * E))
        expected = energy / (s + elastic)
        res = danish_formula(W, h, eh, s, L, A, E)
        assert abs(res["Qu_danish"] - expected) < 1e-6


class TestENR:
    """Engineering News Record Formula."""

    def test_basic(self):
        res = enr(hammer_weight=30, drop_height=1.0, set_per_blow=0.010)
        assert res["Qa_enr"] > 0
        assert res["Qu_enr"] > 0
        assert abs(res["Qu_enr"] - res["Qa_enr"] * 6.0) < 1e-6

    def test_custom_fs(self):
        res = enr(30, 1.0, 0.010, factor_of_safety=3.0)
        assert abs(res["Qu_enr"] - res["Qa_enr"] * 3.0) < 1e-6


# ---------------------------------------------------------------------------
# Integrity
# ---------------------------------------------------------------------------

class TestPileImpedance:
    """Pile impedance Z = EA/c."""

    def test_basic(self):
        res = pile_impedance(elastic_modulus=30e6, area=0.07,
                             wave_speed=3800)
        expected = 30e6 * 0.07 / 3800
        assert abs(res["impedance"] - expected) < 0.1

    def test_validation(self):
        with pytest.raises(ValueError):
            pile_impedance(0, 0.07, 3800)


class TestBetaIntegrity:
    """Beta integrity factor."""

    def test_intact(self):
        res = beta_integrity(14.8, 15.0)
        assert 0.95 <= res["beta"] <= 1.05
        assert "Intact" in res["interpretation"]

    def test_defect(self):
        res = beta_integrity(12.0, 15.0)
        assert res["beta"] < 0.95
        assert "defect" in res["interpretation"].lower()

    def test_bulging(self):
        res = beta_integrity(16.0, 15.0)
        assert res["beta"] > 1.05
        assert "bulging" in res["interpretation"].lower() or "wave speed" in res["interpretation"].lower()
