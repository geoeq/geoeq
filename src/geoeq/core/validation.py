"""
Input validation utilities for GeoEq.

Every public function should validate its inputs before computing.
These helpers raise ``ValueError`` with clear messages when inputs
fall outside physically meaningful ranges.
"""


def check_positive(value: float, name: str) -> None:
    """Raise ``ValueError`` if *value* is not strictly positive."""
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value}")


def check_non_negative(value: float, name: str) -> None:
    """Raise ``ValueError`` if *value* is negative."""
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value}")


def check_range(
    value: float,
    name: str,
    low: float,
    high: float,
    inclusive: bool = True,
) -> None:
    """Raise ``ValueError`` if *value* is outside [low, high] (or (low, high))."""
    if inclusive:
        if value < low or value > high:
            raise ValueError(f"{name} must be in [{low}, {high}], got {value}")
    else:
        if value <= low or value >= high:
            raise ValueError(f"{name} must be in ({low}, {high}), got {value}")


def check_fraction(value: float, name: str) -> None:
    """Raise ``ValueError`` if *value* is not in [0, 1]."""
    check_range(value, name, 0.0, 1.0)


def check_angle_degrees(
    value: float,
    name: str,
    low: float = 0.0,
    high: float = 90.0,
) -> None:
    """Raise ``ValueError`` if an angle in degrees is outside [low, high]."""
    check_range(value, name, low, high)
