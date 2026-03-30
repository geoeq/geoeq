"""Tests for GeoEq (`geoeq.soil.properties`)."""

import pytest
from geoeq.soil import properties as sp


# -------------------------------------------------------------------
# Void ratio / porosity
# -------------------------------------------------------------------


class TestVoidRatioPorosity:
    def test_void_ratio_from_porosity(self):
        assert sp.void_ratio_from_porosity(0.5) == pytest.approx(1.0)
        assert sp.void_ratio_from_porosity(0.4) == pytest.approx(0.4 / 0.6)

    def test_porosity_from_void_ratio(self):
        assert sp.porosity_from_void_ratio(1.0) == pytest.approx(0.5)
        assert sp.porosity_from_void_ratio(0.7) == pytest.approx(0.7 / 1.7)

    def test_roundtrip(self):
        e = 0.85
        n = sp.porosity_from_void_ratio(e)
        e_back = sp.void_ratio_from_porosity(n)
        assert e_back == pytest.approx(e)

    def test_porosity_zero_raises(self):
        with pytest.raises(ValueError):
            sp.void_ratio_from_porosity(0.0)

    def test_porosity_one_raises(self):
        with pytest.raises(ValueError):
            sp.void_ratio_from_porosity(1.0)

    def test_negative_void_ratio_raises(self):
        with pytest.raises(ValueError):
            sp.porosity_from_void_ratio(-0.1)


# -------------------------------------------------------------------
# Degree of saturation and water content
# -------------------------------------------------------------------


class TestSaturation:
    def test_basic(self):
        S = sp.degree_of_saturation(w=0.25, Gs=2.65, e=0.80)
        assert S == pytest.approx(2.65 * 0.25 / 0.80)

    def test_full_saturation(self):
        # e = w * Gs when S = 1
        e = 0.25 * 2.65
        S = sp.degree_of_saturation(w=0.25, Gs=2.65, e=e)
        assert S == pytest.approx(1.0)

    def test_oversaturated_raises(self):
        with pytest.raises(ValueError):
            sp.degree_of_saturation(w=0.9, Gs=2.65, e=0.5)

    def test_zero_water_content(self):
        S = sp.degree_of_saturation(w=0.0, Gs=2.65, e=0.80)
        assert S == pytest.approx(0.0)


class TestWaterContent:
    def test_basic(self):
        w = sp.water_content(S=0.80, Gs=2.65, e=0.72)
        assert w == pytest.approx(0.80 * 0.72 / 2.65)

    def test_roundtrip(self):
        w_in = 0.20
        Gs = 2.65
        e = 0.75
        S = sp.degree_of_saturation(w_in, Gs, e)
        w_out = sp.water_content(S, Gs, e)
        assert w_out == pytest.approx(w_in)


# -------------------------------------------------------------------
# Unit weights
# -------------------------------------------------------------------


class TestUnitWeights:
    def test_dry_unit_weight(self):
        gd = sp.dry_unit_weight(Gs=2.65, e=0.70)
        assert gd == pytest.approx(2.65 * 9.81 / 1.70)

    def test_saturated_unit_weight(self):
        gs = sp.saturated_unit_weight(Gs=2.65, e=0.70)
        assert gs == pytest.approx((2.65 + 0.70) * 9.81 / 1.70)

    def test_bulk_equals_saturated_when_S1(self):
        gb = sp.bulk_unit_weight(Gs=2.65, e=0.70, S=1.0)
        gs = sp.saturated_unit_weight(Gs=2.65, e=0.70)
        assert gb == pytest.approx(gs)

    def test_bulk_equals_dry_when_S0(self):
        gb = sp.bulk_unit_weight(Gs=2.65, e=0.70, S=0.0)
        gd = sp.dry_unit_weight(Gs=2.65, e=0.70)
        assert gb == pytest.approx(gd)

    def test_submerged_unit_weight(self):
        gp = sp.submerged_unit_weight(Gs=2.65, e=0.70)
        assert gp == pytest.approx((2.65 - 1.0) * 9.81 / 1.70)

    def test_submerged_equals_sat_minus_gamma_w(self):
        gs = sp.saturated_unit_weight(Gs=2.65, e=0.70)
        gp = sp.submerged_unit_weight(Gs=2.65, e=0.70)
        assert gp == pytest.approx(gs - 9.81)


# -------------------------------------------------------------------
# Density
# -------------------------------------------------------------------


class TestDensity:
    def test_dry_density(self):
        rd = sp.dry_density(Gs=2.65, e=0.70)
        assert rd == pytest.approx(2.65 * 1000 / 1.70)

    def test_bulk_density(self):
        rb = sp.bulk_density(Gs=2.65, e=0.70, S=0.80)
        assert rb == pytest.approx((2.65 + 0.80 * 0.70) * 1000 / 1.70)


# -------------------------------------------------------------------
# Void ratio from other quantities
# -------------------------------------------------------------------


class TestVoidRatioFrom:
    def test_from_water_content(self):
        e = sp.void_ratio_from_water_content(w=0.25, Gs=2.65, S=1.0)
        assert e == pytest.approx(0.25 * 2.65)

    def test_from_dry_unit_weight(self):
        # Set up a known gamma_d from Gs=2.65, e=0.72
        gd = sp.dry_unit_weight(Gs=2.65, e=0.72)
        e = sp.void_ratio_from_dry_unit_weight(gd, Gs=2.65)
        assert e == pytest.approx(0.72)


# -------------------------------------------------------------------
# Relative density
# -------------------------------------------------------------------


class TestRelativeDensity:
    def test_basic(self):
        Dr = sp.relative_density(e=0.60, e_max=0.90, e_min=0.50)
        assert Dr == pytest.approx(0.75)

    def test_densest(self):
        Dr = sp.relative_density(e=0.50, e_max=0.90, e_min=0.50)
        assert Dr == pytest.approx(1.0)

    def test_loosest(self):
        Dr = sp.relative_density(e=0.90, e_max=0.90, e_min=0.50)
        assert Dr == pytest.approx(0.0)

    def test_invalid_emax_emin(self):
        with pytest.raises(ValueError):
            sp.relative_density(e=0.60, e_max=0.50, e_min=0.90)


# -------------------------------------------------------------------
# Atterberg limits
# -------------------------------------------------------------------


class TestAtterberg:
    def test_plasticity_index(self):
        assert sp.plasticity_index(LL=45, PL=22) == pytest.approx(23)

    def test_plasticity_index_invalid(self):
        with pytest.raises(ValueError):
            sp.plasticity_index(LL=20, PL=30)

    def test_liquidity_index(self):
        LI = sp.liquidity_index(w=35, PL=22, PI=26)
        assert LI == pytest.approx((35 - 22) / 26)

    def test_consistency_index(self):
        CI = sp.consistency_index(w=35, LL=48, PI=26)
        assert CI == pytest.approx((48 - 35) / 26)

    def test_liquidity_plus_consistency_equals_one(self):
        w, LL, PL = 35.0, 48.0, 22.0
        PI = sp.plasticity_index(LL, PL)
        LI = sp.liquidity_index(w, PL, PI)
        CI = sp.consistency_index(w, LL, PI)
        assert LI + CI == pytest.approx(1.0)
