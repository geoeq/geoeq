"""Tests for geoeq.site.pmt — Pressuremeter Test."""

import numpy as np
import pytest
from geoeq.site.pmt import (
    pmt_parameters,
    pmt_modulus,
    pmt_su,
    pmt_bearing,
    pmt_settlement,
    pmt_ko,
)


class TestPmtParameters:
    def test_basic_extraction(self):
        # Linear phase clearly defined
        p = [0, 50, 100, 150, 200, 250, 300, 400, 500, 700, 1000]
        v = [0, 10, 20, 30, 40, 50, 60, 90, 140, 250, 535]
        res = pmt_parameters(p, v)
        assert res["Em"] > 0
        assert res["pL"] > 0
        assert res["p0"] >= 0
        assert res["pf"] > res["p0"]

    def test_short_data_raises(self):
        with pytest.raises(ValueError):
            pmt_parameters([0, 50, 100], [0, 10, 20])


class TestPmtModulus:
    def test_basic(self):
        assert pmt_modulus(15000, alpha=0.5) == 30000.0

    def test_different_alpha(self):
        assert pmt_modulus(15000, alpha=0.33) == pytest.approx(45454.5, rel=1e-2)

    def test_zero_Em_raises(self):
        with pytest.raises(ValueError):
            pmt_modulus(0)


class TestPmtSu:
    def test_basic(self):
        # Su = (600 - 100) / 5.5 = 90.9
        assert pmt_su(pL=600, p0=100) == pytest.approx(90.9, rel=1e-2)

    def test_with_sigma_h0(self):
        # Su = (600 - 80) / 5.5
        assert pmt_su(pL=600, p0=100, sigma_h0=80) == pytest.approx(94.5, rel=1e-2)

    def test_non_negative(self):
        # If pL < p0, Su = 0
        assert pmt_su(pL=50, p0=100) == 0.0


class TestPmtBearing:
    def test_strip(self):
        res = pmt_bearing(pL=800, p0=100, sigma_v=50, shape="strip")
        # qu = 50 + 0.8*(800-100) = 50 + 560 = 610
        assert res["qu"] == pytest.approx(610.0, rel=1e-6)
        assert res["k"] == 0.8

    def test_square(self):
        res = pmt_bearing(pL=800, p0=100, sigma_v=50, shape="square")
        # qu = 50 + 1.2*(800-100) = 50 + 840 = 890
        assert res["qu"] == pytest.approx(890.0, rel=1e-6)
        assert res["k"] == 1.2

    def test_invalid_shape(self):
        with pytest.raises(ValueError):
            pmt_bearing(800, 100, 50, shape="hexagon")


class TestPmtSettlement:
    def test_positive(self):
        s = pmt_settlement(q=150, B=2.0, Em=15000)
        assert s > 0

    def test_increases_with_pressure(self):
        s1 = pmt_settlement(q=100, B=2.0, Em=15000)
        s2 = pmt_settlement(q=200, B=2.0, Em=15000)
        assert s2 > s1

    def test_increases_with_width(self):
        s1 = pmt_settlement(q=150, B=1.0, Em=15000)
        s2 = pmt_settlement(q=150, B=3.0, Em=15000)
        assert s2 > s1


class TestPmtKo:
    def test_basic(self):
        # K0 = (80-30) / (150-30) = 50/120
        assert pmt_ko(80, 150, u0=30) == pytest.approx(50 / 120, rel=1e-6)

    def test_no_pore_pressure(self):
        assert pmt_ko(60, 100) == pytest.approx(0.6, rel=1e-6)

    def test_negative_effective_stress_raises(self):
        with pytest.raises(ValueError):
            pmt_ko(50, 30, u0=40)
