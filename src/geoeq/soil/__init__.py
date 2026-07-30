"""
GeoEq soil property formulas (``geoeq.soil``).

Submodules
----------
properties
    Phase relations, density, Atterberg limits, index properties.
classification
    USCS and AASHTO classification systems, plasticity chart.
texture
    USDA soil texture classification and texture triangle plot.
"""

from . import properties
from . import classification
from . import texture
from geoeq.core.types import Soil
