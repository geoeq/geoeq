"""
Input validation utilities for GeoEq.

Every public function should validate its inputs before computing.
These helpers raise ``ValueError`` with clear messages when inputs
fall outside physically meaningful ranges.

All validators accept both scalars and numpy arrays.
"""

import numpy as np


def check_positive(value, name: str) -> None:
    """Raise ``ValueError`` if any element of *value* is not strictly positive."""
    arr = np.asarray(value, dtype=float)
    if np.any(arr <= 0):
        raise ValueError(f"{name} must be positive, got {value}")


def check_non_negative(value, name: str) -> None:
    """Raise ``ValueError`` if any element of *value* is negative."""
    arr = np.asarray(value, dtype=float)
    if np.any(arr < 0):
        raise ValueError(f"{name} must be non-negative, got {value}")


def check_range(
    value,
    name: str,
    low: float,
    high: float,
    inclusive: bool = True,
) -> None:
    """Raise ``ValueError`` if any element of *value* is outside [low, high] (or (low, high))."""
    arr = np.asarray(value, dtype=float)
    if inclusive:
        if np.any(arr < low) or np.any(arr > high):
            raise ValueError(f"{name} must be in [{low}, {high}], got {value}")
    else:
        if np.any(arr <= low) or np.any(arr >= high):
            raise ValueError(f"{name} must be in ({low}, {high}), got {value}")


def check_fraction(value, name: str) -> None:
    """Raise ``ValueError`` if any element of *value* is not in [0, 1]."""
    check_range(value, name, 0.0, 1.0)


def check_angle_degrees(
    value,
    name: str,
    low: float = 0.0,
    high: float = 90.0,
) -> None:
    """Raise ``ValueError`` if an angle in degrees is outside [low, high]."""
    check_range(value, name, low, high)
