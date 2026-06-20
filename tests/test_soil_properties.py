"""Tests for GeoEq soil properties (``geoeq.soil.properties``).

Every assertion is verified against a hand-calculated textbook example.
Primary reference: Das, B. M. (2021). *Principles of Geotechnical
Engineering*, 10th ed.
"""

import pytest
import numpy as np
from geoeq.soil import properties as sp


# ===================================================================
# Phase relations
# ===================================================================

class TestVoidRatio:
    """Tests for sp.void_ratio() against Das (2021) Ch. 3."""

    def test_from_porosity_basic(self):
        # n = 0.5 → e = 0.5 / (1 - 0.5) = 1.0
        assert sp.void_ratio(n=0.5) == pytest.approx(1.0)

    def test_from_porosity_das(self):
        # Das Ch. 3: n = 0.4 → e = 0.4 / 0.6 = 0.6667
        assert sp.void_ratio(n=0.4) == pytest.approx(0.6667, rel=1e-2)

    def test_from_volumes(self):
        # e = Vv / Vs = 0.4 / 0.6 = 0.6667  [Das Eq. 3.1]
        assert sp.void_ratio(Vv=0.4, Vs=0.6) == pytest.approx(0.6667, rel=1e-2)

    def test_from_w_Gs_S(self):
        # Das Eq. 3.10: Se = wGs → e = wGs/S
        # w=0.25, Gs=2.65, S=1.0 → e = 0.25*2.65/1.0 = 0.6625
        assert sp.void_ratio(w=0.25, Gs=2.65, S=1.0) == pytest.approx(0.6625)

    def test_roundtrip_with_porosity(self):
        e = 0.85
        n = sp.porosity(e=e)
        e_back = sp.void_ratio(n=n)
        assert e_back == pytest.approx(e)

    def test_invalid_porosity_raises(self):
        with pytest.raises(ValueError, match="porosity"):
            sp.void_ratio(n=1.0)
        with pytest.raises(ValueError, match="porosity"):
            sp.void_ratio(n=0.0)

    def test_insufficient_inputs_raises(self):
        with pytest.raises(ValueError, match="Insufficient"):
            sp.void_ratio()
        with pytest.raises(ValueError, match="Insufficient"):
            sp.void_ratio(w=0.2)

    def test_array_input(self):
        n_arr = np.array([0.3, 0.4, 0.5])
        result = sp.void_ratio(n=n_arr)
        expected = n_arr / (1.0 - n_arr)
        np.testing.assert_allclose(result, expected)
        assert result.shape == (3,)

    def test_returns_float_for_scalar(self):
        result = sp.void_ratio(n=0.4)
        assert isinstance(result, float)


class TestPorosity:
    """Tests for sp.porosity() against Das (2021) Ch. 3."""

    def test_from_void_ratio(self):
        # e = 1.0 → n = 1 / 2 = 0.5  [Das Eq. 3.3]
        assert sp.porosity(e=1.0) == pytest.approx(0.5)

    def test_from_void_ratio_typical(self):
        # e = 0.7 → n = 0.7 / 1.7 = 0.4118
        assert sp.porosity(e=0.7) == pytest.approx(0.4118, rel=1e-2)

    def test_from_volumes(self):
        assert sp.porosity(Vv=0.4, V=1.0) == pytest.approx(0.4)

    def test_edge_case_zero_void_ratio(self):
        # e = 0 → n = 0 (solid rock)
        assert sp.porosity(e=0.0) == pytest.approx(0.0)

    def test_insufficient_raises(self):
        with pytest.raises(ValueError, match="Insufficient"):
            sp.porosity()

    def test_array_input(self):
        e_arr = np.array([0.5, 0.7, 1.0])
        result = sp.porosity(e=e_arr)
        expected = e_arr / (1.0 + e_arr)
        np.testing.assert_allclose(result, expected)
        assert result.shape == (3,)


class TestSaturation:
    """Tests for sp.saturation() against Das (2021) Ch. 3."""

    def test_from_w_Gs_e(self):
        # S = wGs/e = 0.25 * 2.65 / 0.80 = 0.828125  [Das Eq. 3.10]
        S = sp.saturation(w=0.25, Gs=2.65, e=0.80)
        assert S == pytest.approx(0.828125)

    def test_from_volumes(self):
        assert sp.saturation(Vw=0.25, Vv=0.50) == pytest.approx(0.5)

    def test_full_saturation(self):
        # w=0.25, Gs=2.65 → e for S=1: e = wGs = 0.6625
        e = 0.25 * 2.65
        S = sp.saturation(w=0.25, Gs=2.65, e=e)
        assert S == pytest.approx(1.0)

    def test_zero_water_content(self):
        S = sp.saturation(w=0.0, Gs=2.65, e=0.80)
        assert S == pytest.approx(0.0)

    def test_oversaturated_raises(self):
        with pytest.raises(ValueError, match="saturation"):
            sp.saturation(w=0.9, Gs=2.65, e=0.5)

    def test_insufficient_raises(self):
        with pytest.raises(ValueError, match="Insufficient"):
            sp.saturation(w=0.2)

    def test_array_input(self):
        w_arr = np.array([0.10, 0.20, 0.30])
        S = sp.saturation(w=w_arr, Gs=2.65, e=0.80)
        expected = w_arr * 2.65 / 0.80
        np.testing.assert_allclose(S, expected)


class TestWaterContent:
    """Tests for sp.water_content() against Das (2021) Ch. 3."""

    def test_from_masses(self):
        # w = Mw/Ms = 20/100 = 0.20  [Das Eq. 3.7]
        assert sp.water_content(Mw=20.0, Ms=100.0) == pytest.approx(0.2)

    def test_from_weights(self):
        assert sp.water_content(Ww=15.0, Ws=75.0) == pytest.approx(0.2)

    def test_from_S_Gs_e(self):
        # w = Se/Gs = 0.80 * 0.72 / 2.65  [Das Eq. 3.10]
        w = sp.water_content(S=0.80, Gs=2.65, e=0.72)
        assert w == pytest.approx(0.80 * 0.72 / 2.65)

    def test_roundtrip_with_saturation(self):
        w_in = 0.20
        Gs, e = 2.65, 0.75
        S = sp.saturation(w=w_in, Gs=Gs, e=e)
        w_out = sp.water_content(S=S, Gs=Gs, e=e)
        assert w_out == pytest.approx(w_in)

    def test_insufficient_raises(self):
        with pytest.raises(ValueError, match="Insufficient"):
            sp.water_content()

    def test_array_input(self):
        Mw = np.array([10, 20, 30], dtype=float)
        Ms = np.array([100, 100, 100], dtype=float)
        w = sp.water_content(Mw=Mw, Ms=Ms)
        np.testing.assert_allclose(w, [0.1, 0.2, 0.3])
        assert w.shape == (3,)


class TestSpecificGravity:
    """Tests for sp.specific_gravity()."""

    def test_basic(self):
        # Gs = Ms / (Vs * rho_w) = 265 / (100 * 1.0) = 2.65
        assert sp.specific_gravity(Ms=265.0, Vs=100.0) == pytest.approx(2.65)

    def test_quartz(self):
        assert sp.specific_gravity(Ms=266.0, Vs=100.0) == pytest.approx(2.66)

    def test_missing_input_raises(self):
        with pytest.raises(ValueError, match="Must provide"):
            sp.specific_gravity(Ms=265.0)

    def test_invalid_volume_raises(self):
        with pytest.raises(ValueError, match="Vs"):
            sp.specific_gravity(Ms=265.0, Vs=-1.0)

    def test_array_input(self):
        Ms = np.array([265.0, 270.0, 275.0])
        Vs = np.array([100.0, 100.0, 100.0])
        Gs = sp.specific_gravity(Ms=Ms, Vs=Vs)
        np.testing.assert_allclose(Gs, [2.65, 2.70, 2.75])

    def test_returns_float_for_scalar(self):
        result = sp.specific_gravity(Ms=265.0, Vs=100.0)
        assert isinstance(result, float)


# ===================================================================
# Density and unit weight
# ===================================================================

class TestDensity:
    """Tests for sp.density() against Das (2021) Ch. 3."""

    def test_dry_unit_weight_das_example(self):
        # Das Example 3.2: Gs=2.67, e=0.72
        # γ_d = Gs * γ_w / (1+e) = 2.67 * 9.81 / 1.72 = 15.24 kN/m³
        gd = sp.density(Gs=2.67, e=0.72, kind="dry")
        assert gd == pytest.approx(15.24, rel=1e-2)

    def test_saturated_unit_weight_das_example(self):
        # Das Example 3.2: Gs=2.67, e=0.72
        # γ_sat = (Gs+e)*γ_w / (1+e) = (2.67+0.72)*9.81 / 1.72 = 19.33 kN/m³
        gs = sp.density(Gs=2.67, e=0.72, kind="saturated")
        assert gs == pytest.approx(19.33, rel=1e-2)

    def test_submerged_unit_weight(self):
        # γ' = (Gs-1)*γ_w / (1+e) = (2.65-1)*9.81 / 1.70 = 9.52 kN/m³
        gp = sp.density(Gs=2.65, e=0.70, kind="submerged")
        assert gp == pytest.approx(9.52, rel=1e-2)

    def test_bulk_unit_weight(self):
        # γ = (Gs+Se)*γ_w / (1+e) for S=1.0 should equal γ_sat
        gb = sp.density(Gs=2.65, e=0.70, S=1.0, kind="bulk")
        gs = sp.density(Gs=2.65, e=0.70, kind="saturated")
        assert gb == pytest.approx(gs)

    def test_bulk_partial_saturation(self):
        # γ = (Gs + S*e)*γ_w / (1+e) = (2.65 + 0.8*0.70)*9.81/1.70
        #   = (2.65 + 0.56)*9.81/1.70 = 3.21*9.81/1.70 = 18.52
        gb = sp.density(Gs=2.65, e=0.70, S=0.80, kind="bulk")
        assert gb == pytest.approx(18.52, rel=1e-2)

    def test_bulk_requires_S(self):
        with pytest.raises(ValueError, match="saturation"):
            sp.density(Gs=2.65, e=0.70, kind="bulk")

    def test_kind_all(self):
        res = sp.density(Gs=2.65, e=0.70, S=0.8, kind="all")
        assert "dry" in res
        assert "saturated" in res
        assert "submerged" in res
        assert "bulk" in res

    def test_dry_density_kg_m3(self):
        # γ_d = Gs * ρ_w / (1+e) = 2.65 * 1000 / 1.70 = 1558.8 kg/m³
        rd = sp.density(Gs=2.65, e=0.70, kind="dry", unit="kg/m3")
        assert rd == pytest.approx(1558.8, rel=1e-2)

    def test_mass_volume(self):
        # 100 kg / 0.1 m³ = 1000 kg/m³
        rd = sp.density(mass=100.0, volume=0.1, unit="kg/m3")
        assert rd == pytest.approx(1000.0)

    def test_mass_volume_kn(self):
        # 1000 kg/m³ * (9.81/1000) = 9.81 kN/m³
        rw = sp.density(mass=100.0, volume=0.1, unit="kN/m3")
        assert rw == pytest.approx(9.81)

    def test_unknown_unit_raises(self):
        with pytest.raises(ValueError, match="Unsupported unit"):
            sp.density(Gs=2.65, e=0.7, kind="dry", unit="lb/ft2")

    def test_unknown_kind_raises(self):
        with pytest.raises(ValueError, match="Unknown kind"):
            sp.density(Gs=2.65, e=0.7, kind="wet")

    def test_insufficient_raises(self):
        with pytest.raises(ValueError, match="Must provide"):
            sp.density(Gs=2.65)

    def test_array_input(self):
        e_arr = np.array([0.5, 0.7, 0.9])
        result = sp.density(Gs=2.65, e=e_arr, kind="dry")
        assert result.shape == (3,)
        assert np.all(result > 0)

    def test_edge_case_zero_void_ratio(self):
        # e = 0 → γ_d = Gs * γ_w = 2.65 * 9.81 = 25.997
        gd = sp.density(Gs=2.65, e=0.0, kind="dry")
        assert gd == pytest.approx(2.65 * 9.81, rel=1e-2)

    def test_returns_float_for_scalar(self):
        result = sp.density(Gs=2.65, e=0.7, kind="dry")
        assert isinstance(result, float)


# ===================================================================
# Relative density
# ===================================================================

class TestRelativeDensity:
    """Tests for sp.relative_density() against Das (2021) Eq. 3.22."""

    def test_basic_void(self):
        # Dr = (e_max - e) / (e_max - e_min) = (0.90-0.60)/(0.90-0.50) = 0.75
        Dr = sp.relative_density(e=0.60, e_max=0.90, e_min=0.50)
        assert Dr == pytest.approx(0.75)

    def test_densest_state(self):
        # e = e_min → Dr = 1.0
        Dr = sp.relative_density(e=0.50, e_max=0.90, e_min=0.50)
        assert Dr == pytest.approx(1.0)

    def test_loosest_state(self):
        # e = e_max → Dr = 0.0
        Dr = sp.relative_density(e=0.90, e_max=0.90, e_min=0.50)
        assert Dr == pytest.approx(0.0)

    def test_invalid_emax_emin_raises(self):
        with pytest.raises(ValueError, match="e_max"):
            sp.relative_density(e=0.60, e_max=0.50, e_min=0.90)

    def test_density_kind(self):
        # Dr = ((ρ-ρ_min)/(ρ_max-ρ_min)) * (ρ_max/ρ)
        # = ((1600-1400)/(1800-1400)) * (1800/1600) = 0.5 * 1.125 = 0.5625
        Dr = sp.relative_density(
            rho=1600, rho_max=1800, rho_min=1400, kind="density"
        )
        assert Dr == pytest.approx(0.5625)

    def test_missing_inputs_raises(self):
        with pytest.raises(ValueError, match="Must provide"):
            sp.relative_density(e=0.6, e_max=0.9, kind="void")

    def test_invalid_kind_raises(self):
        with pytest.raises(ValueError, match="kind"):
            sp.relative_density(e=0.6, e_max=0.9, e_min=0.5, kind="other")

    def test_array_input(self):
        e = np.array([0.50, 0.60, 0.70, 0.90])
        Dr = sp.relative_density(e=e, e_max=0.90, e_min=0.50)
        expected = (0.90 - e) / (0.90 - 0.50)
        np.testing.assert_allclose(Dr, expected)


# ===================================================================
# Atterberg limits & index properties
# ===================================================================

class TestAtterberg:
    """Tests for sp.atterberg() against Das (2021) Ch. 4."""

    def test_plasticity_index(self):
        # PI = LL - PL = 45 - 22 = 23
        assert sp.atterberg(LL=45, PL=22, kind="PI") == pytest.approx(23.0)

    def test_plasticity_index_invalid(self):
        with pytest.raises(ValueError, match="Plastic limit"):
            sp.atterberg(LL=20, PL=30, kind="PI")

    def test_liquidity_index(self):
        # LI = (w - PL) / PI = (35 - 22) / 26 = 0.5
        LI = sp.atterberg(LL=48, PL=22, w=35, kind="LI")
        assert LI == pytest.approx(0.5)

    def test_consistency_index(self):
        # CI = (LL - w) / PI = (48 - 35) / 26 = 0.5
        CI = sp.atterberg(LL=48, PL=22, w=35, kind="CI")
        assert CI == pytest.approx(0.5)

    def test_li_plus_ci_equals_one(self):
        LI = sp.atterberg(LL=48, PL=22, w=35, kind="LI")
        CI = sp.atterberg(LL=48, PL=22, w=35, kind="CI")
        assert LI + CI == pytest.approx(1.0)

    def test_all_indices(self):
        res = sp.atterberg(LL=48, PL=22, w=35, kind="all")
        assert res["PI"] == pytest.approx(26.0)
        assert res["LI"] == pytest.approx(0.5)
        assert res["CI"] == pytest.approx(0.5)

    def test_insufficient_raises(self):
        with pytest.raises(ValueError, match="Insufficient"):
            sp.atterberg(kind="all")

    def test_kind_needs_w(self):
        with pytest.raises(ValueError, match="Cannot calculate"):
            sp.atterberg(LL=45, PL=22, kind="LI")

    def test_array_input(self):
        LL = np.array([40, 50, 60])
        PL = np.array([20, 25, 30])
        PI = sp.atterberg(LL=LL, PL=PL, kind="PI")
        np.testing.assert_allclose(PI, [20, 25, 30])


class TestActivity:
    """Tests for sp.activity() — Skempton (1953)."""

    def test_inactive_clay(self):
        # A = PI / CF = 30 / 40 = 0.75 → boundary of Inactive/Normal
        A = sp.activity(PI=30, clay_fraction=40)
        assert A == pytest.approx(0.75)

    def test_active_clay(self):
        # A = 50 / 20 = 2.5 → Active
        A = sp.activity(PI=50, clay_fraction=20)
        assert A == pytest.approx(2.5)

    def test_normal_clay(self):
        # A = 25 / 25 = 1.0 → Normal
        A = sp.activity(PI=25, clay_fraction=25)
        assert A == pytest.approx(1.0)

    def test_invalid_cf_raises(self):
        with pytest.raises(ValueError, match="clay_fraction"):
            sp.activity(PI=20, clay_fraction=0)

    def test_negative_pi_raises(self):
        with pytest.raises(ValueError, match="PI"):
            sp.activity(PI=-5, clay_fraction=20)

    def test_array_input(self):
        PI = np.array([15, 30, 50])
        CF = np.array([20, 30, 40])
        A = sp.activity(PI=PI, clay_fraction=CF)
        np.testing.assert_allclose(A, [0.75, 1.0, 1.25])
        assert A.shape == (3,)


class TestSensitivity:
    """Tests for sp.sensitivity() — Das Ch. 4."""

    def test_medium_sensitive(self):
        # St = 100 / 25 = 4.0
        St = sp.sensitivity(Su_undisturbed=100, Su_remolded=25)
        assert St == pytest.approx(4.0)

    def test_insensitive(self):
        St = sp.sensitivity(Su_undisturbed=50, Su_remolded=50)
        assert St == pytest.approx(1.0)

    def test_quick_clay(self):
        # St = 80 / 5 = 16 → Quick clay
        St = sp.sensitivity(Su_undisturbed=80, Su_remolded=5)
        assert St == pytest.approx(16.0)

    def test_invalid_su_raises(self):
        with pytest.raises(ValueError, match="Su_undisturbed"):
            sp.sensitivity(Su_undisturbed=-10, Su_remolded=5)

    def test_invalid_su_rem_raises(self):
        with pytest.raises(ValueError, match="Su_remolded"):
            sp.sensitivity(Su_undisturbed=50, Su_remolded=0)

    def test_array_input(self):
        Su_und = np.array([100, 80, 50])
        Su_rem = np.array([25, 10, 50])
        St = sp.sensitivity(Su_undisturbed=Su_und, Su_remolded=Su_rem)
        np.testing.assert_allclose(St, [4.0, 8.0, 1.0])


class TestLiquidityIndex:
    """Tests for sp.liquidity_index() — Das Ch. 4."""

    def test_plastic_range(self):
        # LI = (35 - 22) / 26 = 0.5
        LI = sp.liquidity_index(w=35, PL=22, PI=26)
        assert LI == pytest.approx(0.5)

    def test_at_plastic_limit(self):
        # w = PL → LI = 0
        LI = sp.liquidity_index(w=22, PL=22, PI=26)
        assert LI == pytest.approx(0.0)

    def test_at_liquid_limit(self):
        # w = LL = PL + PI → LI = PI / PI = 1.0
        LI = sp.liquidity_index(w=48, PL=22, PI=26)
        assert LI == pytest.approx(1.0)

    def test_above_liquid_limit(self):
        # w > LL → LI > 1 (liquid state)
        LI = sp.liquidity_index(w=60, PL=22, PI=26)
        assert LI > 1.0

    def test_below_plastic_limit(self):
        # w < PL → LI < 0 (semi-solid state)
        LI = sp.liquidity_index(w=15, PL=22, PI=26)
        assert LI < 0.0

    def test_zero_pi_raises(self):
        with pytest.raises(ValueError, match="PI"):
            sp.liquidity_index(w=30, PL=22, PI=0)

    def test_array_input(self):
        w = np.array([22, 35, 48])
        LI = sp.liquidity_index(w=w, PL=22, PI=26)
        np.testing.assert_allclose(LI, [0.0, 0.5, 1.0])
