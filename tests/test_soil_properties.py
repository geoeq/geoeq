"""Tests for GeoEq (`geoeq.soil.properties`)."""

import pytest
from geoeq.soil import properties as sp


class TestVoidRatioPorosity:
    def test_void_ratio(self):
        assert sp.void_ratio(n=0.5) == pytest.approx(1.0)
        assert sp.void_ratio(n=0.4) == pytest.approx(0.4 / 0.6)
        assert sp.void_ratio(Vv=0.4, Vs=0.6) == pytest.approx(0.4 / 0.6)

    def test_porosity(self):
        assert sp.porosity(e=1.0) == pytest.approx(0.5)
        assert sp.porosity(e=0.7) == pytest.approx(0.7 / 1.7)
        assert sp.porosity(Vv=0.4, V=1.0) == pytest.approx(0.4)

    def test_roundtrip(self):
        e = 0.85
        n = sp.porosity(e=e)
        e_back = sp.void_ratio(n=n)
        assert e_back == pytest.approx(e)


class TestSaturation:
    def test_basic(self):
        S = sp.saturation(w=0.25, Gs=2.65, e=0.80)
        assert S == pytest.approx(2.65 * 0.25 / 0.80)

    def test_from_volumes(self):
        assert sp.saturation(Vw=0.25, Vv=0.50) == pytest.approx(0.5)

    def test_full_saturation(self):
        e = 0.25 * 2.65
        S = sp.saturation(w=0.25, Gs=2.65, e=e)
        assert S == pytest.approx(1.0)

    def test_oversaturated_raises(self):
        with pytest.raises(ValueError):
            sp.saturation(w=0.9, Gs=2.65, e=0.5)

    def test_zero_water_content(self):
        S = sp.saturation(w=0.0, Gs=2.65, e=0.80)
        assert S == pytest.approx(0.0)


class TestWaterContent:
    def test_basic(self):
        w = sp.water_content(S=0.80, Gs=2.65, e=0.72)
        assert w == pytest.approx(0.80 * 0.72 / 2.65)

    def test_from_masses(self):
        assert sp.water_content(Mw=20.0, Ms=100.0) == pytest.approx(0.2)

    def test_roundtrip(self):
        w_in = 0.20
        Gs = 2.65
        e = 0.75
        S = sp.saturation(w=w_in, Gs=Gs, e=e)
        w_out = sp.water_content(S=S, Gs=Gs, e=e)
        assert w_out == pytest.approx(w_in)


class TestDensityAndUnitWeights:
    def test_dry_unit_weight(self):
        gd = sp.density(Gs=2.65, e=0.70, kind="dry", unit="kN/m3")
        assert gd == pytest.approx(2.65 * 9.81 / 1.70)

    def test_saturated_unit_weight(self):
        gs = sp.density(Gs=2.65, e=0.70, kind="saturated", unit="kN/m3")
        assert gs == pytest.approx((2.65 + 0.70) * 9.81 / 1.70)

    def test_bulk_unit_weight(self):
        gb = sp.density(Gs=2.65, e=0.70, S=1.0, kind="bulk", unit="kN/m3")
        gs = sp.density(Gs=2.65, e=0.70, kind="saturated", unit="kN/m3")
        assert gb == pytest.approx(gs)

    def test_submerged_unit_weight(self):
        gp = sp.density(Gs=2.65, e=0.70, kind="submerged", unit="kN/m3")
        assert gp == pytest.approx((2.65 - 1.0) * 9.81 / 1.70)

    def test_dry_density(self):
        rd = sp.density(Gs=2.65, e=0.70, kind="dry", unit="kg/m3")
        assert rd == pytest.approx(2.65 * 1000 / 1.70)

    def test_bulk_density(self):
        rb = sp.density(Gs=2.65, e=0.70, S=0.80, kind="bulk", unit="kg/m3")
        assert rb == pytest.approx((2.65 + 0.80 * 0.70) * 1000 / 1.70)
        
    def test_mass_volume_density(self):
        # 100 kg / 0.1 m3 = 1000 kg/m3
        rd = sp.density(mass=100.0, volume=0.1, unit="kg/m3")
        assert rd == pytest.approx(1000.0)
        
        # In kN/m3: 1000 kg/m3 * (9.81 / 1000) = 9.81
        rw = sp.density(mass=100.0, volume=0.1, unit="kN/m3")
        assert rw == pytest.approx(9.81)


class TestVoidRatioFrom:
    def test_from_water_content(self):
        e = sp.void_ratio(w=0.25, Gs=2.65, S=1.0)
        assert e == pytest.approx(0.25 * 2.65)


class TestRelativeDensity:
    def test_basic_void(self):
        Dr = sp.relative_density(e=0.60, e_max=0.90, e_min=0.50, kind="void")
        assert Dr == pytest.approx(0.75)

    def test_densest_void(self):
        Dr = sp.relative_density(e=0.50, e_max=0.90, e_min=0.50, kind="void")
        assert Dr == pytest.approx(1.0)

    def test_invalid_emax_emin(self):
        with pytest.raises(ValueError):
            sp.relative_density(e=0.60, e_max=0.50, e_min=0.90, kind="void")
            
    def test_density_kind(self):
        # Let's say rho_min = 1400, rho_max = 1800, rho = 1600
        # Dr = ((1600 - 1400) / (1800 - 1400)) * (1800 / 1600)
        # Dr = (200 / 400) * (1800 / 1600) = 0.5 * 1.125 = 0.5625
        Dr = sp.relative_density(rho=1600, rho_max=1800, rho_min=1400, kind="density")
        assert Dr == pytest.approx(0.5625)


class TestAtterberg:
    def test_plasticity_index(self):
        assert sp.atterberg(LL=45, PL=22, kind="PI") == pytest.approx(23)

    def test_plasticity_index_invalid(self):
        with pytest.raises(ValueError):
            sp.atterberg(LL=20, PL=30, kind="PI")

    def test_liquidity_index(self):
        LI = sp.atterberg(LL=48, PL=22, w=35, kind="LI")
        assert LI == pytest.approx((35 - 22) / 26)

    def test_consistency_index(self):
        CI = sp.atterberg(LL=48, PL=22, w=35, kind="CI")
        assert CI == pytest.approx((48 - 35) / 26)

    def test_liquidity_plus_consistency_equals_one(self):
        w, LL, PL = 35.0, 48.0, 22.0
        LI = sp.atterberg(LL=LL, PL=PL, w=w, kind="LI")
        CI = sp.atterberg(LL=LL, PL=PL, w=w, kind="CI")
        assert LI + CI == pytest.approx(1.0)
        
    def test_all_indices(self):
        res = sp.atterberg(LL=48, PL=22, w=35, kind="all")
        assert res["PI"] == pytest.approx(26)
        assert res["LI"] == pytest.approx((35 - 22) / 26)
        assert res["CI"] == pytest.approx((48 - 35) / 26)
