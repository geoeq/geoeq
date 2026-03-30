"""
GeoEq -- A Python library for geotechnical engineering calculations.

Install the package as ``geoeq`` (``pip install geoeq``); import ``geoeq`` in code.

Version 0.0.3 provides soil property and phase-relationship formulas.
Future versions will add effective stress, bearing capacity, settlement,
consolidation, earth pressure, slope stability, pile capacity,
and liquefaction assessment.

Quick start
-----------
>>> from geoeq.soil import properties as sp
>>> sp.dry_unit_weight(Gs=2.65, e=0.72)
15.11...
"""

__version__ = "0.0.3"

from . import core
from . import soil
