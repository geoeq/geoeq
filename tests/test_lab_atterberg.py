"""Tests for GeoEq Atterberg limits lab tests (``geoeq.lab.atterberg_test``).

References:
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 4.
    ASTM D4318.
"""

import pytest
import numpy as np

import matplotlib
matplotlib.use("Agg")

from geoeq.lab.atterberg_test import liquid_limit_test, flow_curve_plot


# ===================================================================
# Liquid Limit Test
# ===================================================================

class TestLiquidLimitTest:
    """Tests for liquid_limit_test() — flow curve processing."""

    def test_basic(self):
        # 4 trials: LL should be at N=25
        N = [15, 20, 28, 34]
        w = [42.1, 40.8, 38.5, 36.9]
        res = liquid_limit_test(N, w)
        assert 37 < res["LL"] < 41
        assert res["r_squared"] > 0.95

    def test_two_point(self):
        res = liquid_limit_test([20, 30], [42, 38])
        assert 38 < res["LL"] < 42

    def test_known_ll(self):
        # If data is perfectly linear on semi-log, w at N=25 is exact
        # w = 50 - 10*log10(N), at N=25: w = 50 - 10*1.3979 = 36.02
        N = [10, 20, 30, 40]
        w = [50 - 10 * np.log10(n) for n in N]
        res = liquid_limit_test(N, w)
        assert res["LL"] == pytest.approx(36.02, rel=0.01)
        assert res["r_squared"] == pytest.approx(1.0, abs=1e-6)

    def test_negative_slope(self):
        # w should decrease with N → slope is negative
        res = liquid_limit_test([15, 20, 28, 34], [42.1, 40.8, 38.5, 36.9])
        assert res["slope"] < 0

    def test_too_few_points(self):
        with pytest.raises(ValueError, match="at least 2"):
            liquid_limit_test([25], [40])

    def test_length_mismatch(self):
        with pytest.raises(ValueError):
            liquid_limit_test([15, 20, 28], [42, 40])


# ===================================================================
# Flow Curve Plot (smoke tests)
# ===================================================================

class TestFlowCurvePlot:
    def test_returns_dict(self):
        res = flow_curve_plot([15, 20, 28, 34], [42.1, 40.8, 38.5, 36.9])
        assert "LL" in res
        assert "ax" in res
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_custom_axes(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        res = flow_curve_plot([15, 20, 28, 34], [42.1, 40.8, 38.5, 36.9], ax=ax)
        assert res["ax"] is ax
        plt.close("all")
