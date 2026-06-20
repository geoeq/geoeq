"""Tests for GeoEq permeability lab tests (``geoeq.lab.permeability``).

References:
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 5.
    ASTM D2434, D5084.
"""

import pytest
import numpy as np

import matplotlib
matplotlib.use("Agg")

from geoeq.lab.permeability import constant_head, falling_head, permeability_plot


# ===================================================================
# Constant-Head
# ===================================================================

class TestConstantHead:
    """Tests for constant_head() — Darcy's law (Das Eq. 5.11)."""

    def test_basic(self):
        # k = Q*L / (A*h*t) = 500*15 / (30*50*120) = 7500/180000 = 0.04167 cm/s
        k = constant_head(Q=500, L=15, A=30, h=50, t=120)
        assert k == pytest.approx(0.04167, rel=1e-3)

    def test_das_example(self):
        # Das (2021) Example 5.2: Q=620cm³, L=25cm, A=25cm², h=70cm, t=300s
        # k = 620*25 / (25*70*300) = 15500/525000 = 0.02952 cm/s
        k = constant_head(Q=620, L=25, A=25, h=70, t=300)
        assert k == pytest.approx(0.02952, rel=1e-2)

    def test_array_q(self):
        Q = np.array([500, 1000])
        k = constant_head(Q=Q, L=15, A=30, h=50, t=120)
        assert k.shape == (2,)
        assert k[1] == pytest.approx(2 * k[0], rel=1e-6)

    def test_invalid_negative_q(self):
        with pytest.raises(ValueError, match="positive"):
            constant_head(Q=-100, L=15, A=30, h=50, t=120)

    def test_invalid_zero_area(self):
        with pytest.raises(ValueError, match="positive"):
            constant_head(Q=500, L=15, A=0, h=50, t=120)


# ===================================================================
# Falling-Head
# ===================================================================

class TestFallingHead:
    """Tests for falling_head() — Das Eq. 5.13."""

    def test_basic(self):
        # k = (a*L)/(A*t) * ln(h1/h2)
        # = (1.0*15)/(30*600) * ln(100/50) = 15/18000 * 0.6931 = 0.000577 cm/s
        k = falling_head(a=1.0, L=15, A=30, h1=100, h2=50, t=600)
        assert k == pytest.approx(5.776e-4, rel=1e-2)

    def test_h2_greater_than_h1_error(self):
        with pytest.raises(ValueError, match="h2 must be less"):
            falling_head(a=1.0, L=15, A=30, h1=50, h2=100, t=600)

    def test_equal_heads_error(self):
        with pytest.raises(ValueError, match="h2 must be less"):
            falling_head(a=1.0, L=15, A=30, h1=100, h2=100, t=600)

    def test_invalid_negative(self):
        with pytest.raises(ValueError, match="positive"):
            falling_head(a=-1, L=15, A=30, h1=100, h2=50, t=600)

    def test_small_head_ratio(self):
        # h1/h2 ≈ 1.01, small k
        k = falling_head(a=1.0, L=15, A=30, h1=100, h2=99, t=600)
        assert k > 0
        assert k < 1e-3


# ===================================================================
# Permeability Plot (smoke tests)
# ===================================================================

class TestPermeabilityPlot:
    def test_constant_head_plot(self):
        res = permeability_plot(
            Q_values=[10, 20, 30],
            head_gradient=[1.0, 2.0, 3.0],
            test_type="constant",
        )
        assert "ax" in res
        assert "k" in res
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_falling_head_plot(self):
        res = permeability_plot(
            time_values=[0, 60, 120, 180, 240],
            head_values=[100, 85, 72, 62, 53],
            test_type="falling",
        )
        assert "ax" in res
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_invalid_type(self):
        with pytest.raises(ValueError, match="test_type"):
            permeability_plot(test_type="invalid")
