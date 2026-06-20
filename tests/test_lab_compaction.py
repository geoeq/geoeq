"""Tests for GeoEq compaction lab tests (``geoeq.lab.compaction``).

References:
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 5.
    ASTM D698, D1557.
"""

import pytest
import numpy as np

import matplotlib
matplotlib.use("Agg")

from geoeq.lab.compaction import (
    proctor,
    zav_line,
    saturation_line,
    relative_compaction,
    proctor_plot,
)


# ===================================================================
# Proctor
# ===================================================================

class TestProctor:
    """Tests for proctor() — compaction curve fitting."""

    def test_basic_optimum(self):
        # Symmetric parabola: peak near w=14
        w = [8, 10, 12, 14, 16, 18]
        gd = [17.5, 18.2, 18.8, 19.0, 18.7, 18.1]
        res = proctor(w, gd)
        assert 12.5 < res["w_opt"] < 15.5
        assert 18.5 < res["gamma_d_max"] < 19.5

    def test_das_example_compaction(self):
        # Typical Proctor test data (Das 2021)
        w = [10, 12, 14, 16, 18, 20]
        gd = [16.8, 17.6, 18.2, 18.4, 18.0, 17.4]
        res = proctor(w, gd)
        assert 14 < res["w_opt"] < 18
        assert res["gamma_d_max"] > 18.0

    def test_coeffs_returned(self):
        w = [8, 10, 12, 14, 16, 18]
        gd = [17.5, 18.2, 18.8, 19.0, 18.7, 18.1]
        res = proctor(w, gd)
        assert len(res["coeffs"]) == 3
        assert res["coeffs"][0] < 0  # parabola opens downward

    def test_too_few_points(self):
        with pytest.raises(ValueError, match="at least 3"):
            proctor([10, 14], [17, 18])

    def test_length_mismatch(self):
        with pytest.raises(ValueError):
            proctor([10, 12, 14], [17, 18])

    def test_negative_density(self):
        with pytest.raises(ValueError, match="positive"):
            proctor([10, 12, 14], [-1, 18, 19])


# ===================================================================
# ZAV Line
# ===================================================================

class TestZAVLine:
    """Tests for zav_line() — S = 100% theoretical line."""

    def test_zav_decreasing(self):
        # γ_d should decrease as w increases
        res = zav_line(Gs=2.70)
        assert all(np.diff(res["dry_density"]) < 0)

    def test_zav_specific_value(self):
        # At w=10% (0.10): γ_d = 2.70*9.81 / (1 + 0.10*2.70) = 26.487/1.27 = 20.86
        res = zav_line(Gs=2.70, w_range=[10])
        assert res["dry_density"][0] == pytest.approx(20.86, rel=0.01)

    def test_custom_range(self):
        res = zav_line(Gs=2.65, w_range=[5, 10, 15, 20])
        assert len(res["water_content"]) == 4
        assert len(res["dry_density"]) == 4

    def test_invalid_gs(self):
        with pytest.raises(ValueError, match="positive"):
            zav_line(Gs=-1)


# ===================================================================
# Saturation Line
# ===================================================================

class TestSaturationLine:
    def test_80_percent_below_zav(self):
        zav = zav_line(Gs=2.70, w_range=[10, 15, 20])
        sat80 = saturation_line(Gs=2.70, S=0.8, w_range=[10, 15, 20])
        # S=80% line should be below S=100% (ZAV is the upper bound)
        assert all(sat80["dry_density"] < zav["dry_density"])

    def test_s100_equals_zav(self):
        zav = zav_line(Gs=2.70, w_range=[10, 15, 20])
        sat100 = saturation_line(Gs=2.70, S=1.0, w_range=[10, 15, 20])
        np.testing.assert_allclose(sat100["dry_density"], zav["dry_density"], rtol=1e-6)

    def test_invalid_s(self):
        with pytest.raises(ValueError):
            saturation_line(Gs=2.70, S=1.5)


# ===================================================================
# Relative Compaction
# ===================================================================

class TestRelativeCompaction:
    def test_basic(self):
        # RC = 17.5 / 19.0 * 100 = 92.1%
        assert relative_compaction(17.5, 19.0) == pytest.approx(92.1, rel=0.01)

    def test_100_percent(self):
        assert relative_compaction(19.0, 19.0) == pytest.approx(100.0)

    def test_array_input(self):
        gd = np.array([17.0, 18.0, 19.0])
        rc = relative_compaction(gd, 19.0)
        assert rc.shape == (3,)
        assert rc[-1] == pytest.approx(100.0)

    def test_invalid_negative(self):
        with pytest.raises(ValueError, match="positive"):
            relative_compaction(-1, 19.0)


# ===================================================================
# Proctor Plot (smoke tests)
# ===================================================================

class TestProctorPlot:
    def test_returns_dict(self):
        w = [8, 10, 12, 14, 16, 18]
        gd = [17.5, 18.2, 18.8, 19.0, 18.7, 18.1]
        res = proctor_plot(w, gd, Gs=2.70)
        assert "w_opt" in res
        assert "gamma_d_max" in res
        assert "ax" in res
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_no_zav(self):
        w = [8, 10, 12, 14, 16, 18]
        gd = [17.5, 18.2, 18.8, 19.0, 18.7, 18.1]
        res = proctor_plot(w, gd, show_zav=False, show_sat_lines=False)
        assert "ax" in res
        import matplotlib.pyplot as plt
        plt.close("all")
