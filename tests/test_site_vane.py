"""Tests for geoeq.site.vane — Field Vane Shear Test."""

import numpy as np
import pytest
from geoeq.site.vane import vane_su, vane_correction, vane_remolded


class TestVaneSu:
    def test_standard_vane(self):
        # D=0.065m, H=0.130m (H/D=2)
        # K = pi * (D²H/2 + D³/6) = pi * (0.065²*0.13/2 + 0.065³/6)
        D, H = 0.065, 0.130
        K = np.pi * (D**2 * H / 2 + D**3 / 6)
        T = 0.012
        expected = T / K
        result = vane_su(T, D, H)
        assert result == pytest.approx(expected, rel=1e-3)

    def test_hd_ratio_2(self):
        # For H/D=2: Su = 6T / (7*pi*D³)
        D = 0.065
        H = 2 * D
        T = 0.015
        Su_formula = 6 * T / (7 * np.pi * D**3)
        Su_general = vane_su(T, D, H)
        assert Su_general == pytest.approx(Su_formula, rel=1e-3)

    def test_zero_torque(self):
        assert vane_su(0.0, 0.065, 0.130) == 0.0

    def test_array_torque(self):
        T = np.array([0.01, 0.015, 0.020])
        result = vane_su(T, 0.065, 0.130)
        assert len(result) == 3
        assert all(r > 0 for r in result)


class TestVaneCorrection:
    def test_pi_40(self):
        # mu = 1.7 - 0.54*log10(40) ≈ 1.7 - 0.54*1.602 = 0.835
        mu = 1.7 - 0.54 * np.log10(40)
        expected = mu * 50
        assert vane_correction(50, PI=40) == pytest.approx(expected, rel=1e-2)

    def test_pi_100(self):
        # mu = 1.7 - 0.54*log10(100) = 1.7 - 1.08 = 0.62
        mu = 1.7 - 0.54 * 2.0
        expected = mu * 80
        assert vane_correction(80, PI=100) == pytest.approx(expected, rel=1e-2)

    def test_mu_capped_at_1(self):
        # Very low PI → mu would exceed 1, but capped
        result = vane_correction(50, PI=5)
        assert result == pytest.approx(50.0, rel=1e-2)

    def test_mu_minimum_0_5(self):
        # Very high PI → mu could go below 0.5, but capped
        result = vane_correction(100, PI=10000)
        assert result == pytest.approx(50.0, rel=1e-2)


class TestVaneRemolded:
    def test_sensitivity(self):
        res = vane_remolded(T_peak=0.015, T_remolded=0.003, D=0.065, H=0.130)
        assert res["sensitivity"] == pytest.approx(5.0, rel=1e-2)
        assert res["classification"] == "Medium sensitivity"

    def test_insensitive(self):
        res = vane_remolded(T_peak=0.010, T_remolded=0.008, D=0.065, H=0.130)
        assert res["sensitivity"] < 2
        assert res["classification"] == "Insensitive"

    def test_quick_clay(self):
        res = vane_remolded(T_peak=0.020, T_remolded=0.001, D=0.065, H=0.130)
        assert res["sensitivity"] >= 16
        assert res["classification"] == "Quick clay"
