"""Tests for geoeq.site.cpt — CPT processing and correlations."""

import numpy as np
import pytest
from geoeq.site.cpt import (
    cpt_normalize,
    cpt_ic,
    cpt_sbt,
    cpt_friction_angle,
    cpt_su,
    cpt_dr,
    cpt_modulus,
    cpt_sbt_plot,
)


class TestCptNormalize:
    def test_basic(self):
        res = cpt_normalize(qt=5000, fs=50, sigma_v=100, sigma_v_eff=60)
        assert res["Qt"] == pytest.approx(81.67, rel=1e-2)
        assert res["Fr"] == pytest.approx(50 / 4900 * 100, rel=1e-2)
        assert res["Bq"] == pytest.approx(0.0, abs=1e-6)

    def test_with_pore_pressure(self):
        res = cpt_normalize(qt=5000, fs=50, sigma_v=100, sigma_v_eff=60, u2=200, u0=50)
        assert res["Bq"] == pytest.approx(150 / 4900, rel=1e-2)

    def test_array_input(self):
        qt = np.array([3000, 5000, 10000])
        fs = np.array([30, 50, 80])
        sv = np.array([50, 100, 150])
        sve = np.array([30, 60, 100])
        res = cpt_normalize(qt, fs, sv, sve)
        assert len(res["Qt"]) == 3


class TestCptIc:
    def test_clean_sand(self):
        # Qt=80, Fr=1.0 → Ic ≈ 1.99 (zone 6)
        ic = cpt_ic(Qt=80, Fr=1.0)
        assert ic == pytest.approx(1.986, rel=0.01)

    def test_clay(self):
        # High Fr, low Qt → high Ic
        ic = cpt_ic(Qt=5, Fr=5.0)
        assert ic > 2.95

    def test_array(self):
        Qt = np.array([80, 20, 5])
        Fr = np.array([1.0, 3.0, 5.0])
        ic = cpt_ic(Qt, Fr)
        assert len(ic) == 3


class TestCptSbt:
    def test_sand_zone_6(self):
        res = cpt_sbt(Qt=80, Fr=1.0)
        assert res["zone"] == 6

    def test_clay_zone_3(self):
        res = cpt_sbt(Qt=5, Fr=3.0)
        assert res["zone"] in [3, 4]

    def test_gravelly_zone_7(self):
        res = cpt_sbt(Qt=500, Fr=0.3)
        assert res["zone"] == 7

    def test_array(self):
        res = cpt_sbt(Qt=np.array([80, 5]), Fr=np.array([1.0, 5.0]))
        assert len(res["zone"]) == 2


class TestCptFrictionAngle:
    def test_typical_sand(self):
        # Robertson & Campanella: arctan(0.1 + 0.38*log(5000/60))
        phi = cpt_friction_angle(qt=5000, sigma_v_eff=60)
        assert 30 < phi < 45

    def test_increases_with_qt(self):
        phi1 = cpt_friction_angle(3000, 60)
        phi2 = cpt_friction_angle(10000, 60)
        assert phi2 > phi1


class TestCptSu:
    def test_basic(self):
        # Su = (1000 - 100) / 14 = 64.29
        assert cpt_su(1000, 100, Nkt=14) == pytest.approx(64.286, rel=1e-2)

    def test_non_negative(self):
        # If qt < sigma_v, Su should be 0
        assert cpt_su(50, 100) == 0.0

    def test_array(self):
        qt = np.array([500, 1000, 2000])
        sv = np.array([50, 100, 150])
        result = cpt_su(qt, sv)
        assert len(result) == 3


class TestCptDr:
    def test_typical_sand(self):
        Dr = cpt_dr(10000, 100)
        assert 0 < Dr < 1

    def test_increases_with_qc(self):
        Dr1 = cpt_dr(5000, 100)
        Dr2 = cpt_dr(20000, 100)
        assert Dr2 > Dr1

    def test_capped_at_1(self):
        Dr = cpt_dr(1e6, 50)
        assert Dr <= 1.0


class TestCptModulus:
    def test_with_alpha(self):
        # M = 5 * (5000 - 100) = 24500
        assert cpt_modulus(5000, 100, alpha=5) == pytest.approx(24500.0, rel=1e-6)

    def test_with_Ic_sand(self):
        # Ic=1.5 (sand) → alpha_M = 0.0188 * 10^(0.55*1.5 + 1.68)
        M = cpt_modulus(5000, 100, Ic=1.5)
        assert M > 0

    def test_requires_ic_or_alpha(self):
        with pytest.raises(ValueError):
            cpt_modulus(5000, 100)


class TestCptSbtPlot:
    def test_smoke(self):
        import matplotlib
        matplotlib.use("Agg")
        res = cpt_sbt_plot([80, 20, 5], [1.0, 3.0, 5.0])
        assert "ax" in res
