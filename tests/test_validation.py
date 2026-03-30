"""Tests for GeoEq (`geoeq.core.validation`)."""

import pytest
from geoeq.core import validation as v


class TestCheckPositive:
    def test_passes(self):
        v.check_positive(1.0, "x")
        v.check_positive(0.001, "x")

    def test_zero_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            v.check_positive(0.0, "x")

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            v.check_positive(-1.0, "x")


class TestCheckNonNegative:
    def test_passes(self):
        v.check_non_negative(0.0, "x")
        v.check_non_negative(5.0, "x")

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="must be non-negative"):
            v.check_non_negative(-0.1, "x")


class TestCheckRange:
    def test_inclusive_passes(self):
        v.check_range(0.0, "x", 0.0, 1.0)
        v.check_range(1.0, "x", 0.0, 1.0)

    def test_inclusive_fails(self):
        with pytest.raises(ValueError):
            v.check_range(1.1, "x", 0.0, 1.0)

    def test_exclusive_passes(self):
        v.check_range(0.5, "x", 0.0, 1.0, inclusive=False)

    def test_exclusive_boundary_fails(self):
        with pytest.raises(ValueError):
            v.check_range(0.0, "x", 0.0, 1.0, inclusive=False)


class TestCheckFraction:
    def test_passes(self):
        v.check_fraction(0.0, "x")
        v.check_fraction(0.5, "x")
        v.check_fraction(1.0, "x")

    def test_above_one_raises(self):
        with pytest.raises(ValueError):
            v.check_fraction(1.1, "x")

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            v.check_fraction(-0.1, "x")
