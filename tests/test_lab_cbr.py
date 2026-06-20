"""Tests for GeoEq CBR lab tests (``geoeq.lab.cbr``).

References:
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 5.
    ASTM D1883.
"""

import pytest
import numpy as np

import matplotlib
matplotlib.use("Agg")

from geoeq.lab.cbr import cbr_test, cbr_plot


# ===================================================================
# CBR Test
# ===================================================================

class TestCBRTest:
    """Tests for cbr_test() — California Bearing Ratio."""

    def test_basic(self):
        pen = [0, 0.64, 1.27, 2.54, 3.81, 5.08, 7.62, 10.16, 12.70]
        load = [0, 0.8, 1.9, 4.2, 6.1, 7.8, 10.5, 12.1, 13.2]
        res = cbr_test(pen, load)
        assert res["CBR"] > 0
        assert res["CBR_2_54"] > 0
        assert res["CBR_5_08"] > 0

    def test_known_cbr_at_2_54(self):
        # load at 2.54mm = 4.2 kN
        # CBR = 4.2 / 13.24 * 100 = 31.72%
        pen = [0, 2.54, 5.08]
        load = [0, 4.2, 7.8]
        res = cbr_test(pen, load)
        assert res["CBR_2_54"] == pytest.approx(31.72, rel=0.01)

    def test_known_cbr_at_5_08(self):
        # load at 5.08mm = 7.8 kN
        # CBR = 7.8 / 19.96 * 100 = 39.08%
        pen = [0, 2.54, 5.08]
        load = [0, 4.2, 7.8]
        res = cbr_test(pen, load)
        assert res["CBR_5_08"] == pytest.approx(39.08, rel=0.01)

    def test_cbr_is_max(self):
        pen = [0, 2.54, 5.08]
        load = [0, 4.2, 7.8]
        res = cbr_test(pen, load)
        assert res["CBR"] == max(res["CBR_2_54"], res["CBR_5_08"])

    def test_too_few_points(self):
        with pytest.raises(ValueError, match="at least 3"):
            cbr_test([0, 2.54], [0, 4.2])

    def test_length_mismatch(self):
        with pytest.raises(ValueError):
            cbr_test([0, 1, 2, 3], [0, 1, 2])


# ===================================================================
# CBR Plot (smoke tests)
# ===================================================================

class TestCBRPlot:
    def test_returns_dict(self):
        pen = [0, 0.64, 1.27, 2.54, 3.81, 5.08, 7.62, 10.16, 12.70]
        load = [0, 0.8, 1.9, 4.2, 6.1, 7.8, 10.5, 12.1, 13.2]
        res = cbr_plot(pen, load)
        assert "CBR" in res
        assert "ax" in res
        import matplotlib.pyplot as plt
        plt.close("all")
