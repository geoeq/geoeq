"""Tests for USDA soil texture classification (geoeq.soil.texture)."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

import geoeq as ge


# -------------------------------------------------------------------
# classify_usda — one representative point per class
# -------------------------------------------------------------------

@pytest.mark.parametrize(
    "sand, silt, clay, expected",
    [
        (92, 5, 3, "Sand"),
        (82, 12, 6, "Loamy Sand"),
        (60, 30, 10, "Sandy Loam"),
        (40, 40, 20, "Loam"),
        (20, 65, 15, "Silt Loam"),
        (5, 88, 7, "Silt"),
        (60, 15, 25, "Sandy Clay Loam"),
        (33, 33, 34, "Clay Loam"),
        (10, 55, 35, "Silty Clay Loam"),
        (55, 10, 35, "Sandy Clay"),
        (5, 50, 45, "Silty Clay"),
        (20, 20, 60, "Clay"),
    ],
)
def test_classify_usda_all_classes(sand, silt, clay, expected):
    result = ge.classify_usda(sand=sand, silt=silt, clay=clay)
    assert result["texture"] == expected


def test_classify_usda_returns_inputs():
    result = ge.classify_usda(sand=40, silt=40, clay=20)
    assert result["sand"] == 40
    assert result["silt"] == 40
    assert result["clay"] == 20


# -------------------------------------------------------------------
# boundary behaviour (exact USDA inequalities)
# -------------------------------------------------------------------

def test_sand_loamy_sand_boundary():
    # silt + 1.5*clay = 15 is exactly NOT sand (strict <)
    assert ge.classify_usda(sand=85, silt=15, clay=0)["texture"] == "Loamy Sand"
    assert ge.classify_usda(sand=86, silt=14, clay=0)["texture"] == "Sand"


def test_clay_minimum_40_percent():
    assert ge.classify_usda(sand=30, silt=30, clay=40)["texture"] == "Clay"
    assert ge.classify_usda(sand=31, silt=30, clay=39)["texture"] == "Clay Loam"


def test_silt_requires_80_percent():
    assert ge.classify_usda(sand=10, silt=80, clay=10)["texture"] == "Silt"
    assert ge.classify_usda(sand=11, silt=79, clay=10)["texture"] == "Silt Loam"


def test_silty_clay_boundary():
    assert ge.classify_usda(sand=0, silt=40, clay=60)["texture"] == "Silty Clay"
    assert ge.classify_usda(sand=0, silt=60, clay=40)["texture"] == "Silty Clay"


# -------------------------------------------------------------------
# validation
# -------------------------------------------------------------------

def test_sum_must_be_100():
    with pytest.raises(ValueError, match="must equal 100"):
        ge.classify_usda(sand=50, silt=30, clay=10)


def test_negative_fraction_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        ge.classify_usda(sand=-5, silt=60, clay=45)


def test_tolerance_accepts_rounding():
    # 33.3 + 33.3 + 33.4 = 100.0 exactly; 33+33+33.5 = 99.5 within atol
    result = ge.classify_usda(sand=33.0, silt=33.0, clay=33.5)
    assert result["texture"] == "Clay Loam"


# -------------------------------------------------------------------
# full-triangle sweep — every valid point must classify
# -------------------------------------------------------------------

def test_no_unclassified_gaps():
    for sand in range(0, 101, 2):
        for clay in range(0, 101 - sand, 2):
            silt = 100 - sand - clay
            result = ge.classify_usda(sand=sand, silt=silt, clay=clay)
            assert result["texture"] != "Unclassified", (
                f"gap at sand={sand}, silt={silt}, clay={clay}"
            )


# -------------------------------------------------------------------
# texture_triangle plot
# -------------------------------------------------------------------

def test_texture_triangle_returns_axes():
    ax = ge.texture_triangle()
    assert ax is not None
    plt.close("all")


def test_texture_triangle_with_points():
    ax = ge.texture_triangle(sand=[40, 92], silt=[40, 5], clay=[20, 3],
                             labels=["A", "B"])
    assert ax is not None
    plt.close("all")


def test_texture_triangle_single_point_scalar():
    ax = ge.texture_triangle(sand=40, silt=40, clay=20, labels="S1")
    assert ax is not None
    plt.close("all")


def test_texture_triangle_save(tmp_path):
    out = tmp_path / "tri.png"
    ge.texture_triangle(sand=33, silt=33, clay=34, save=str(out))
    assert out.exists()
    plt.close("all")


def test_texture_triangle_invalid_point_raises():
    with pytest.raises(ValueError):
        ge.texture_triangle(sand=50, silt=30, clay=10)
    plt.close("all")
