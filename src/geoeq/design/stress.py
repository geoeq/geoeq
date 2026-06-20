"""
Effective stress and pore-water pressure.

References
----------
* Terzaghi, K. (1943). *Theoretical Soil Mechanics*. Wiley.
* Das, B. M. (2010). *Principles of Geotechnical Engineering* (7th ed.), Ch. 5.
* Budhu, M. (2011). *Soil Mechanics and Foundations* (3rd ed.), Ch. 7.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Union

import numpy as np

from geoeq.core.constants import GAMMA_WATER
from geoeq.core.validation import check_positive, check_non_negative


def total_stress(
    gamma: Union[float, Sequence[float]],
    depth: Union[float, Sequence[float]],
) -> float:
    """Total vertical stress sigma_v = sum(gamma_i * H_i).

    Parameters
    ----------
    gamma : float or sequence
        Unit weight(s) of each layer (kN/m^3).
    depth : float or sequence
        Either a single depth (with scalar gamma) -- returns gamma*depth --
        or layer thicknesses matching ``gamma`` -- returns sum(gamma_i * H_i).

    Returns
    -------
    sigma : float
        Total vertical stress (kPa).

    Reference
    ---------
    Das (2010) Eq. 5.1; Terzaghi (1943).
    """
    g = np.atleast_1d(np.asarray(gamma, dtype=float))
    d = np.atleast_1d(np.asarray(depth, dtype=float))
    check_positive(g, "gamma")
    check_non_negative(d, "depth")
    if g.size == 1 and d.size == 1:
        return float(g[0] * d[0])
    if g.size != d.size:
        raise ValueError(
            "gamma and depth must have equal length (one entry per layer).")
    return float(np.sum(g * d))


def pore_pressure(
    z: Union[float, Iterable[float]],
    z_w: float = 0.0,
    gamma_w: float = GAMMA_WATER,
) -> Union[float, np.ndarray]:
    """Hydrostatic pore pressure u = gamma_w * (z - z_w), clipped at 0.

    Parameters
    ----------
    z : float or array
        Depth(s) of interest (m, positive downward).
    z_w : float
        Depth of the water table (m). Default 0 (water table at surface).
    gamma_w : float
        Unit weight of water (kN/m^3). Default 9.81.

    Returns
    -------
    u : float or array
        Pore pressure (kPa). Zero above the water table.

    Reference
    ---------
    Das (2010) Eq. 5.2.
    """
    check_positive(gamma_w, "gamma_w")
    z_arr = np.asarray(z, dtype=float)
    u = gamma_w * np.maximum(0.0, z_arr - z_w)
    return float(u) if u.ndim == 0 else u


def effective_stress(
    sigma: Union[float, np.ndarray],
    u: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    """Terzaghi effective stress sigma' = sigma - u (kPa).

    Reference
    ---------
    Terzaghi (1943); Das (2010) Eq. 5.5.
    """
    sigma = np.asarray(sigma, dtype=float)
    u = np.asarray(u, dtype=float)
    res = sigma - u
    return float(res) if res.ndim == 0 else res


def capillary_rise(D10: float, e: float, C: float = 0.2) -> float:
    """Capillary rise height in soil (Hazen-type estimate).

    h_c [cm] = C / (e * D10 [cm])

    Parameters
    ----------
    D10 : float
        Effective grain size, the 10% passing diameter (mm).
    e : float
        Void ratio (-).
    C : float
        Empirical constant; 0.1 < C < 0.5 typical. Default 0.2.

    Returns
    -------
    h_c : float
        Capillary rise (m).

    Reference
    ---------
    Hazen (1930); Das (2010) Eq. 5.8.
    """
    check_positive(D10, "D10")
    check_positive(e, "e")
    D10_cm = D10 / 10.0  # mm to cm
    h_c_cm = C / (e * D10_cm)
    return h_c_cm / 100.0  # to m


def stress_plot(profile, dz: float = 0.1, **kwargs):
    """Plot total, pore, and effective stress vs depth.

    Delegates to ``profile.plot()``. See ``SoilProfile.plot``.
    """
    return profile.plot(dz=dz, **kwargs)
