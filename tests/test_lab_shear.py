"""Tests for GeoEq shear strength lab tests (``geoeq.lab.shear``).

References:
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 8.
    ASTM D3080, D2850, D4767, D2166.
"""

import pytest
import numpy as np

import matplotlib
matplotlib.use("Agg")

from geoeq.lab.shear import (
    direct_shear,
    triaxial,
    unconfined,
    mohr_circle,
    direct_shear_plot,
)


# ===================================================================
# Direct Shear
# ===================================================================

class TestDirectShear:
    """Tests for direct_shear() — Mohr–Coulomb fitting."""

    def test_basic_three_points(self):
        # Hand calc: σ = [50, 100, 150], τ = [38, 62, 86]
        # Slope = (86-38)/(150-50) = 0.48 → φ = arctan(0.48) = 25.6°
        # Intercept ≈ 38 - 0.48*50 = 14
        res = direct_shear([50, 100, 150], [38, 62, 86])
        assert res["phi"] == pytest.approx(25.6, abs=1.0)
        assert res["c"] == pytest.approx(14.0, abs=2.0)
        assert res["r_squared"] > 0.99

    def test_cohesionless_soil(self):
        # Pure sand: τ = σ * tan(30°) = σ * 0.577
        sigma = [50, 100, 200]
        tau = [28.87, 57.74, 115.47]
        res = direct_shear(sigma, tau)
        assert res["phi"] == pytest.approx(30.0, abs=0.5)
        assert res["c"] == pytest.approx(0.0, abs=1.0)

    def test_purely_cohesive(self):
        # φ ≈ 0: τ is roughly constant
        res = direct_shear([50, 100, 200], [25, 26, 27])
        assert res["phi"] < 2.0
        assert res["c"] > 20

    def test_minimum_points_error(self):
        with pytest.raises(ValueError, match="at least 3"):
            direct_shear([50, 100], [30, 50])

    def test_length_mismatch(self):
        with pytest.raises(ValueError):
            direct_shear([50, 100, 150], [30, 50])

    def test_four_points(self):
        res = direct_shear([50, 100, 150, 200], [30, 50, 70, 90])
        assert 0 < res["phi"] < 45
        assert res["c"] >= 0
        assert res["r_squared"] > 0.95


class TestDirectShearPlot:
    def test_returns_dict_with_ax(self):
        res = direct_shear_plot([50, 100, 150], [38, 62, 86])
        assert "ax" in res
        assert "phi" in res
        import matplotlib.pyplot as plt
        plt.close("all")


# ===================================================================
# Triaxial
# ===================================================================

class TestTriaxial:
    """Tests for triaxial() — UU, CU, CD processing."""

    def test_cd_basic(self):
        # σ₃ = [100, 200, 300], Δσ = [220, 380, 540]
        # σ₁ = [320, 580, 840]
        # p = [210, 390, 570], q = [110, 190, 270]
        # sin(φ) = slope of q vs p ≈ (270-110)/(570-210) = 0.444 → φ ≈ 26.4°
        res = triaxial([100, 200, 300], [220, 380, 540], kind="CD")
        assert res["phi"] == pytest.approx(26.4, abs=1.5)
        assert res["c"] >= 0

    def test_uu_returns_su(self):
        # UU: φ = 0, Su = Δσ / 2
        # Δσ = [80, 82, 78] → avg Su = 40
        res = triaxial([100, 200, 300], [80, 82, 78], kind="UU")
        assert res["phi"] == 0.0
        assert res["Su"] == pytest.approx(40.0, abs=2.0)
        assert res["c"] == pytest.approx(40.0, abs=2.0)

    def test_cu_test(self):
        res = triaxial([100, 200, 300], [250, 420, 590], kind="CU")
        assert 0 < res["phi"] < 45
        assert res["c"] >= 0

    def test_invalid_kind(self):
        with pytest.raises(ValueError, match="kind"):
            triaxial([100, 200], [200, 350], kind="XX")

    def test_uu_needs_2_points(self):
        with pytest.raises(ValueError, match="at least 2"):
            triaxial([100], [80], kind="UU")

    def test_cd_needs_3_points(self):
        with pytest.raises(ValueError, match="at least 3"):
            triaxial([100, 200], [200, 350], kind="CD")

    def test_sigma1_computed(self):
        res = triaxial([100, 200], [80, 85], kind="UU")
        np.testing.assert_array_almost_equal(res["sigma1"], [180, 285])
        np.testing.assert_array_almost_equal(res["sigma3"], [100, 200])


# ===================================================================
# Unconfined Compression
# ===================================================================

class TestUnconfined:
    """Tests for unconfined() — Su = qu/2."""

    def test_basic(self):
        # qu = 120 kPa → Su = 60 kPa → "Stiff"
        res = unconfined(120)
        assert res["Su"] == 60.0
        assert res["consistency"] == "Stiff"

    def test_very_soft(self):
        res = unconfined(20)
        assert res["Su"] == 10.0
        assert res["consistency"] == "Very soft"

    def test_hard(self):
        res = unconfined(500)
        assert res["Su"] == 250.0
        assert res["consistency"] == "Hard"

    def test_array_input(self):
        res = unconfined([20, 120, 500])
        assert len(res["Su"]) == 3
        assert len(res["consistency"]) == 3

    def test_invalid_qu(self):
        with pytest.raises(ValueError, match="positive"):
            unconfined(-10)


# ===================================================================
# Mohr Circle
# ===================================================================

class TestMohrCircle:
    """Smoke tests for mohr_circle() plotting."""

    def test_returns_c_phi(self):
        res = mohr_circle([320, 580, 840], [100, 200, 300])
        assert res["c"] is not None
        assert res["phi"] is not None
        assert res["phi"] == pytest.approx(26.4, abs=2.0)
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_two_circles_no_envelope(self):
        res = mohr_circle([200, 400], [50, 100])
        assert res["c"] is None
        assert res["phi"] is None
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_single_circle(self):
        res = mohr_circle(320, 100)
        assert res["ax"] is not None
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_custom_axes(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        res = mohr_circle([320, 580, 840], [100, 200, 300], ax=ax)
        assert res["ax"] is ax
        plt.close("all")
