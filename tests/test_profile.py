"""Tests for ge.profile -- Soil and SoilProfile."""

import numpy as np
import pytest

from geoeq.profile import Soil, SoilProfile, mesh


def test_soil_defaults():
    s = Soil("X")
    assert s.gamma == 18.0
    assert s.gamma_sat == pytest.approx(19.5)


def test_soil_effective_unit_weight():
    s = Soil("Sand", gamma=18, gamma_sat=20)
    assert s.gamma_effective() == pytest.approx(20 - 9.81)


def test_profile_contiguity_check():
    with pytest.raises(ValueError, match="contiguous"):
        SoilProfile([(0, 2, Soil("a")), (3, 5, Soil("b"))])


def test_profile_layer_order_check():
    with pytest.raises(ValueError, match="bottom"):
        SoilProfile([(0, 0, Soil("a"))])


def test_profile_empty():
    with pytest.raises(ValueError, match="at least one"):
        SoilProfile([])


def test_total_stress_dry_profile():
    p = SoilProfile([(0, 5, Soil("S", gamma=20))], water_table=None)
    assert p.total_stress(5) == pytest.approx(100.0)
    assert p.pore_pressure(5) == 0.0
    assert p.effective_stress(5) == pytest.approx(100.0)


def test_total_stress_multi_layer_with_WT():
    # Das example: 2 m fill (gamma=17) over 6 m sand (gamma_sat=20),
    # water table at 2 m. At z = 8 m:
    #   sigma = 17*2 + 20*6 = 34 + 120 = 154 kPa
    #   u     = 9.81 * 6 = 58.86 kPa
    #   sigma'= 154 - 58.86 = 95.14 kPa
    p = SoilProfile([
        (0, 2, Soil("Fill", gamma=17)),
        (2, 8, Soil("Sand", gamma=18, gamma_sat=20)),
    ], water_table=2.0)
    assert p.total_stress(8) == pytest.approx(154.0)
    assert p.pore_pressure(8) == pytest.approx(58.86, abs=0.05)
    assert p.effective_stress(8) == pytest.approx(95.14, abs=0.05)


def test_stress_array_input():
    p = SoilProfile([(0, 10, Soil("S", gamma=18, gamma_sat=20))], water_table=5)
    z = np.array([0, 5, 10])
    sigma = p.total_stress(z)
    assert sigma[0] == 0.0
    assert sigma[1] == pytest.approx(90.0)  # 18*5
    assert sigma[2] == pytest.approx(190.0)  # 18*5 + 20*5


def test_stress_at_returns_dict():
    p = SoilProfile([(0, 5, Soil("S", gamma=20))], water_table=None)
    d = p.stress_at(5)
    assert set(d.keys()) == {"sigma", "u", "sigma_eff"}
    assert d["sigma_eff"] == pytest.approx(100.0)


def test_layer_at():
    p = SoilProfile([
        (0, 2, Soil("A")), (2, 5, Soil("B"))])
    assert p.layer_at(1).name == "A"
    assert p.layer_at(3).name == "B"
    with pytest.raises(ValueError):
        p.layer_at(10)


def test_add_layer():
    p = SoilProfile([(0, 2, Soil("A"))])
    p.add_layer(2, 5, Soil("B"))
    assert p.n_layers == 2
    with pytest.raises(ValueError, match="must start"):
        p.add_layer(10, 12, Soil("X"))


def test_mesh():
    p = SoilProfile([(0, 5, Soil("a"))])
    grid = mesh(p, dz=1.0)
    assert grid[0] == 0
    assert grid[-1] == pytest.approx(5.0)
    assert len(grid) == 6


def test_to_dataframe():
    pd = pytest.importorskip("pandas")
    p = SoilProfile([(0, 2, Soil("A")), (2, 5, Soil("B"))])
    df = p.to_dataframe()
    assert len(df) == 2
    assert "thickness" in df.columns


def test_plot_returns_figure():
    p = SoilProfile([
        (0, 2, Soil("Fill", gamma=17)),
        (2, 8, Soil("Clay", gamma=18, gamma_sat=19)),
    ], water_table=2)
    fig = p.plot(dz=0.5)
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close(fig)
