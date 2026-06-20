"""Tests for GeoEq consolidation lab tests (``geoeq.lab.consolidation``).

References:
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 7.
    ASTM D2435.
"""

import pytest
import numpy as np

import matplotlib
matplotlib.use("Agg")

from geoeq.lab.consolidation import (
    oedometer,
    preconsolidation,
    compression_index,
    cv,
    oedometer_plot,
    cv_plot,
)


# ===================================================================
# Oedometer
# ===================================================================

class TestOedometer:
    """Tests for oedometer() — e-log(p) processing."""

    def test_basic_cc(self):
        # Das-style data: stress doubles, e decreases
        # Between 200→400: Δe = 0.12, Δlog = log(400/200) = 0.301
        # Cc = 0.12 / 0.301 = 0.399; between 400→800: 0.14/0.301 = 0.465
        stress = [25, 50, 100, 200, 400, 800]
        e = [0.88, 0.86, 0.82, 0.74, 0.62, 0.48]
        res = oedometer(stress, e)
        assert res["Cc"] == pytest.approx(0.465, rel=0.05)
        assert res["Cr"] > 0

    def test_cr_from_first_two(self):
        stress = [25, 50, 100, 200, 400, 800]
        e = [0.88, 0.86, 0.82, 0.74, 0.62, 0.48]
        res = oedometer(stress, e)
        # Cr = (0.88-0.86) / log(50/25) = 0.02/0.301 = 0.066
        assert res["Cr"] == pytest.approx(0.066, rel=0.1)

    def test_sorted_output(self):
        # Pass unsorted
        stress = [200, 25, 800, 100, 50, 400]
        e = [0.74, 0.88, 0.48, 0.82, 0.86, 0.62]
        res = oedometer(stress, e)
        assert list(res["stress"]) == sorted(res["stress"])

    def test_too_few_points(self):
        with pytest.raises(ValueError, match="at least 3"):
            oedometer([100, 200], [0.8, 0.7])

    def test_negative_stress(self):
        with pytest.raises(ValueError, match="positive"):
            oedometer([-10, 50, 100], [0.9, 0.85, 0.8])

    def test_length_mismatch(self):
        with pytest.raises(ValueError):
            oedometer([25, 50, 100], [0.88, 0.86])


# ===================================================================
# Preconsolidation
# ===================================================================

class TestPreconsolidation:
    """Tests for preconsolidation() — Casagrande method."""

    def test_pc_reasonable_range(self):
        stress = [25, 50, 100, 200, 400, 800]
        e = [0.88, 0.86, 0.82, 0.74, 0.62, 0.48]
        res = preconsolidation(stress, e)
        assert 30 < res["pc"] < 300

    def test_method_name(self):
        stress = [25, 50, 100, 200, 400, 800]
        e = [0.88, 0.86, 0.82, 0.74, 0.62, 0.48]
        res = preconsolidation(stress, e)
        assert res["method"] == "casagrande"

    def test_too_few_points(self):
        with pytest.raises(ValueError, match="at least 4"):
            preconsolidation([50, 100, 200], [0.86, 0.82, 0.74])

    def test_invalid_method(self):
        with pytest.raises(ValueError, match="not supported"):
            preconsolidation([25, 50, 100, 200], [0.9, 0.88, 0.84, 0.78],
                             method="invalid")


# ===================================================================
# Compression Index
# ===================================================================

class TestCompressionIndex:
    """Tests for compression_index() — empirical correlations."""

    def test_terzaghi(self):
        # Cc = 0.009 * (50 - 10) = 0.36
        assert compression_index("terzaghi", LL=50) == pytest.approx(0.36, rel=1e-3)

    def test_skempton(self):
        # Cc = 0.007 * (50 - 7) = 0.301
        assert compression_index("skempton", LL=50) == pytest.approx(0.301, rel=1e-3)

    def test_rendon(self):
        # Cc = 0.01 * 45 = 0.45
        assert compression_index("rendon", wn=45) == pytest.approx(0.45, rel=1e-3)

    def test_nishida(self):
        # Cc = 1.15 * (0.8 - 0.35) = 0.5175
        assert compression_index("nishida", e0=0.8) == pytest.approx(0.5175, rel=1e-3)

    def test_nagaraj(self):
        # Cc = 0.2343 * 0.8 * 2.65 = 0.4967
        assert compression_index("nagaraj", e0=0.8, Gs=2.65) == pytest.approx(0.4967, rel=1e-2)

    def test_hough(self):
        # Cc = 0.29 * (0.8 - 0.27) = 0.1537
        assert compression_index("hough", e0=0.8) == pytest.approx(0.1537, rel=1e-2)

    def test_missing_parameter(self):
        with pytest.raises(ValueError, match="LL"):
            compression_index("terzaghi")

    def test_unknown_method(self):
        with pytest.raises(ValueError, match="Unknown"):
            compression_index("bogus", LL=50)


# ===================================================================
# Coefficient of Consolidation
# ===================================================================

class TestCv:
    """Tests for cv() — log-time and root-time methods."""

    @pytest.fixture()
    def consol_data(self):
        t = [0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30, 60, 120, 240, 480, 1440]
        d = [0.0, 0.05, 0.12, 0.22, 0.35, 0.50, 0.65, 0.78, 0.88, 0.95,
             0.99, 1.02, 1.04, 1.06]
        return t, d

    def test_log_method_positive(self, consol_data):
        t, d = consol_data
        res = cv(t, d, method="log", H_dr=1.0)
        assert res["cv"] > 0
        assert "t50" in res
        assert res["method"] == "log"

    def test_root_method_positive(self, consol_data):
        t, d = consol_data
        res = cv(t, d, method="root", H_dr=1.0)
        assert res["cv"] > 0
        assert "t90" in res
        assert res["method"] == "root"

    def test_invalid_method(self, consol_data):
        t, d = consol_data
        with pytest.raises(ValueError, match="method"):
            cv(t, d, method="invalid")

    def test_too_few_points(self):
        with pytest.raises(ValueError, match="at least 5"):
            cv([1, 2, 3], [0.1, 0.2, 0.3])

    def test_h_dr_effect(self, consol_data):
        t, d = consol_data
        res1 = cv(t, d, method="log", H_dr=1.0)
        res2 = cv(t, d, method="log", H_dr=2.0)
        # cv scales with H_dr²
        assert res2["cv"] == pytest.approx(res1["cv"] * 4, rel=0.01)


# ===================================================================
# Plots (smoke tests)
# ===================================================================

class TestConsolidationPlots:
    def test_oedometer_plot(self):
        res = oedometer_plot([25, 50, 100, 200, 400, 800],
                             [0.88, 0.86, 0.82, 0.74, 0.62, 0.48])
        assert "Cc" in res
        assert "ax" in res
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_cv_plot_log(self):
        t = [0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30, 60, 120, 240, 480, 1440]
        d = [0.0, 0.05, 0.12, 0.22, 0.35, 0.50, 0.65, 0.78, 0.88, 0.95,
             0.99, 1.02, 1.04, 1.06]
        res = cv_plot(t, d, method="log")
        assert "cv" in res
        assert "ax" in res
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_cv_plot_root(self):
        t = [0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30, 60, 120, 240, 480, 1440]
        d = [0.0, 0.05, 0.12, 0.22, 0.35, 0.50, 0.65, 0.78, 0.88, 0.95,
             0.99, 1.02, 1.04, 1.06]
        res = cv_plot(t, d, method="root")
        assert "cv" in res
        import matplotlib.pyplot as plt
        plt.close("all")
