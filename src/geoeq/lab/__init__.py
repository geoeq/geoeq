"""
GeoEq laboratory testing module (``geoeq.lab``).

Submodules
----------
shear
    Direct shear, triaxial, unconfined compression, Mohr circle.
consolidation
    Oedometer, preconsolidation pressure, compression index, cv.
compaction
    Proctor test, zero air voids line, relative compaction.
permeability
    Constant-head and falling-head permeability tests.
atterberg
    Liquid limit (flow curve), plastic limit.
cbr
    California Bearing Ratio test processing.
"""

from . import shear
from . import consolidation
from . import compaction
from . import permeability
from . import atterberg_test
from . import cbr
