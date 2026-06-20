"""
GeoEq soil property formulas (``geoeq.soil``).

Submodules
----------
properties
    Phase relations, density, Atterberg limits, index properties.
classification
    USCS and AASHTO classification systems, plasticity chart.
"""

from . import properties
from . import classification
from geoeq.core.types import Soil
