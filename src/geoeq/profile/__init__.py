"""
ge.profile -- Soil profiles and layered data structures.

The connective tissue between lab tests, site investigation, and design.
A ``SoilProfile`` is a vertical stack of ``Soil`` layers with an optional
water table. It computes total, pore, and effective stress at any depth
and exposes plotting and DataFrame export.

References
----------
Das, B. M. (2010). *Principles of Geotechnical Engineering* (7th ed.), Ch. 5.
Budhu, M. (2011). *Soil Mechanics and Foundations* (3rd ed.).
"""

from geoeq.profile.soil import Soil
from geoeq.profile.profile import SoilProfile, mesh, log_plot

__all__ = ["Soil", "SoilProfile", "mesh", "log_plot"]
