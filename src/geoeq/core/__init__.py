"""
GeoEq core (``geoeq.core``) -- constants, validation, and types.
"""

from .constants import GRAVITY, GAMMA_WATER, DENSITY_WATER, ATMOSPHERIC_PRESSURE
from .validation import (
    check_positive,
    check_non_negative,
    check_range,
    check_fraction,
    check_angle_degrees,
)
from .types import Soil
