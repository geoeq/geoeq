"""
Common data structures for GeoEq.

Using dataclasses keeps data structured and inspectable rather than
passing around bare dictionaries.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class Soil:
    """
    Represents a soil with basic engineering properties.

    Parameters
    ----------
    gamma : float
        Total unit weight (kN/m^3).  Default 18.0.
    phi : float
        Effective friction angle (degrees).  Default 30.0.
    c : float
        Effective cohesion (kPa).  Default 0.0.
    e : float or None
        Void ratio.  Default ``None``.
    Gs : float
        Specific gravity of solids.  Default 2.65.
    w : float or None
        Water content as a decimal (e.g. 0.25 for 25 %).  Default ``None``.
    gamma_d : float or None
        Dry unit weight (kN/m^3).  Default ``None``.
    Su : float or None
        Undrained shear strength (kPa).  Default ``None``.
    description : str
        Free-text label.  Default ``""``.
    """

    gamma: float = 18.0
    phi: float = 30.0
    c: float = 0.0
    e: Optional[float] = None
    Gs: float = 2.65
    w: Optional[float] = None
    gamma_d: Optional[float] = None
    Su: Optional[float] = None
    description: str = ""
