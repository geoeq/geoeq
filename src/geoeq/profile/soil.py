"""
Soil layer object used inside a SoilProfile.

A ``Soil`` is a plain dataclass holding the engineering properties of
one homogeneous deposit. The convention follows Das (2010):

* ``gamma``     -- moist/bulk unit weight above the water table (kN/m^3)
* ``gamma_sat`` -- saturated unit weight below the water table (kN/m^3)
* ``phi``       -- effective friction angle (degrees)
* ``c``         -- effective cohesion (kPa)
* ``e``         -- void ratio (-)
* ``Gs``        -- specific gravity of solids (-)
* ``k``         -- hydraulic conductivity (m/s)
* ``Es``        -- elastic modulus (kPa)
* ``mu``        -- Poisson's ratio (-)
* ``Cc``        -- compression index (-)
* ``Cr``        -- recompression index (-)
* ``OCR``       -- over-consolidation ratio (-)
* ``Su``        -- undrained shear strength (kPa)

Only ``name`` is required. ``gamma`` defaults to 18 kN/m^3 (a generic
soil) which is enough to compute stresses if nothing else is supplied.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Soil:
    """A single soil layer's engineering properties.

    Examples
    --------
    >>> from geoeq.profile import Soil
    >>> clay = Soil("Soft Clay", gamma=17, gamma_sat=18.5, phi=0, c=25, e=0.9)
    >>> clay.name
    'Soft Clay'
    """

    name: str = "Soil"
    gamma: float = 18.0
    gamma_sat: Optional[float] = None
    phi: float = 30.0
    c: float = 0.0
    e: Optional[float] = None
    Gs: float = 2.70
    k: Optional[float] = None
    Es: Optional[float] = None
    mu: float = 0.30
    Cc: Optional[float] = None
    Cr: Optional[float] = None
    OCR: float = 1.0
    Su: Optional[float] = None
    description: str = ""

    def __post_init__(self) -> None:
        # Default gamma_sat: estimate from gamma if not given (assumes saturation).
        if self.gamma_sat is None:
            self.gamma_sat = self.gamma + 1.5  # crude default; user should override

    def gamma_effective(self) -> float:
        """Submerged (buoyant) unit weight gamma' = gamma_sat - gamma_w."""
        from geoeq.core.constants import GAMMA_WATER
        return float(self.gamma_sat) - GAMMA_WATER

    def __repr__(self) -> str:  # pragma: no cover -- cosmetic
        return (f"Soil({self.name!r}, gamma={self.gamma}, "
                f"phi={self.phi}, c={self.c})")
