"""Tests for geoeq.site.spt — SPT corrections and correlations."""

import numpy as np
import pytest
from geoeq.site.spt import (
    spt_n60,
    spt_n160,
    spt_n160cs,
    spt_friction_angle,
    spt_su,
    spt_dr,
    spt_modulus,
)


class TestSptN60:
    def test_default_energy_60(self):
        assert spt_n60(30, ER=60) == 30.0

    def test_auto_hammer_72(self):
        assert spt_n60(30, ER=72) == pytest.approx(36.0, rel=1e-2)

    def test_all_corrections(self):
        result = spt_n60(30, ER=72, Cb=1.05, Cs=1.0, Cr=0.95)
        expected = 30 * (72 / 60) * 1.05 * 1.0 * 0.95
        assert result == pytest.approx(expected, rel=1e-6)

    def test_array_input(self):
        N = np.array([10, 20, 30])
        result = spt_n60(N, ER=80)
        expected = N * (80 / 60)
        np.testing.assert_allclose(result, expected)

    def test_zero_N(self):
        assert spt_n60(0) == 0.0

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            spt_n60(-5)


class TestSptN160:
    def test_liao_whitman_at_pa(self):
        # At sigma_v = pa, CN = 1.0
        assert spt_n160(20, sigma_v=100) == pytest.approx(20.0, rel=1e-2)

    def test_liao_whitman_half_pa(self):
        # At sigma_v = 50, CN = sqrt(100/50) = 1.414
        assert spt_n160(20, sigma_v=50) == pytest.approx(28.28, rel=1e-2)

    def test_skempton_method(self):
        # CN = 2 / (1 + 100/100) = 1.0
        assert spt_n160(20, sigma_v=100, method="skempton") == pytest.approx(20.0, rel=1e-2)

    def test_peck_method(self):
        CN = 0.77 * np.log10(2000 / 100)
        expected = 20 * CN
        assert spt_n160(20, sigma_v=100, method="peck") == pytest.approx(expected, rel=1e-2)

    def test_cn_capped_at_2(self):
        # Very low stress → CN would be > 2
        result = spt_n160(10, sigma_v=10)
        assert result == pytest.approx(20.0, rel=1e-2)  # 10 * 2.0 (capped)

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            spt_n160(20, sigma_v=100, method="bogus")

    def test_array_input(self):
        N60 = np.array([10, 20, 30])
        sv = np.array([50, 100, 200])
        result = spt_n160(N60, sv)
        assert len(result) == 3


class TestSptN160cs:
    def test_low_fines(self):
        # FC ≤ 5%: alpha=0, beta=1 → no correction
        assert spt_n160cs(15, FC=3) == pytest.approx(15.0, rel=1e-2)

    def test_high_fines(self):
        # FC ≥ 35%: alpha=5.0, beta=1.2
        assert spt_n160cs(15, FC=40) == pytest.approx(5.0 + 1.2 * 15, rel=1e-2)

    def test_intermediate_fines(self):
        # FC = 20% → alpha = exp(1.76 - 190/400), beta = 0.99 + 20^1.5/1000
        alpha = np.exp(1.76 - 190.0 / 400)
        beta = 0.99 + 20**1.5 / 1000
        expected = alpha + beta * 15
        assert spt_n160cs(15, FC=20) == pytest.approx(expected, rel=1e-2)


class TestSptFrictionAngle:
    def test_hatanaka_n20(self):
        # phi = sqrt(20*20) + 20 = 20 + 20 = 40
        assert spt_friction_angle(20, method="hatanaka") == pytest.approx(40.0, rel=1e-2)

    def test_hatanaka_n5(self):
        # phi = sqrt(100) + 20 = 10 + 20 = 30
        assert spt_friction_angle(5, method="hatanaka") == pytest.approx(30.0, rel=1e-2)

    def test_peck_method(self):
        result = spt_friction_angle(15, method="peck")
        assert 26 < result < 45

    def test_kulhawy_method(self):
        result = spt_friction_angle(20, sigma_v=100, method="kulhawy")
        assert 25 < result < 50

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            spt_friction_angle(20, method="bogus")


class TestSptSu:
    def test_hara(self):
        # Su = 29 * 10^0.72
        expected = 29.0 * 10**0.72
        assert spt_su(10, method="hara") == pytest.approx(expected, rel=1e-2)

    def test_stroud_default(self):
        # default PI=None → f1=5.0
        assert spt_su(10, method="stroud") == pytest.approx(50.0, rel=1e-2)

    def test_stroud_high_PI(self):
        # PI=50 → f1=6.0
        assert spt_su(10, method="stroud", PI=50) == pytest.approx(60.0, rel=1e-2)

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            spt_su(10, method="bogus")


class TestSptDr:
    def test_meyerhof(self):
        # Dr = sqrt(20/41) * 100
        expected = np.sqrt(20 / 41.0) * 100
        assert spt_dr(20, method="meyerhof") == pytest.approx(expected, rel=1e-2)

    def test_skempton(self):
        expected = np.sqrt(20 / 55.0) * 100
        assert spt_dr(20, method="skempton") == pytest.approx(expected, rel=1e-2)

    def test_kulhawy(self):
        expected = np.sqrt(20 / 60.0) * 100
        assert spt_dr(20, method="kulhawy") == pytest.approx(expected, rel=1e-2)

    def test_capped_at_100(self):
        assert spt_dr(100, method="meyerhof") == pytest.approx(100.0, rel=1e-2)

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            spt_dr(20, method="bogus")


class TestSptModulus:
    def test_sand(self):
        # Es = 500 * (20 + 15) = 17500
        assert spt_modulus(20, "sand") == pytest.approx(17500.0, rel=1e-6)

    def test_gravel(self):
        assert spt_modulus(20, "gravel") == pytest.approx(1200 * 26, rel=1e-6)

    def test_clay_soft(self):
        assert spt_modulus(10, "clay_soft") == pytest.approx(300 * 16, rel=1e-6)

    def test_invalid_soil(self):
        with pytest.raises(ValueError):
            spt_modulus(20, "marble")

    def test_array(self):
        result = spt_modulus(np.array([10, 20, 30]), "sand")
        expected = 500 * (np.array([10, 20, 30]) + 15)
        np.testing.assert_allclose(result, expected)
