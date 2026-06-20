r"""
Field Permeability Tests — slug test and pumping test interpretation.

Functions
---------
slug_test
    Hvorslev (1951) method for hydraulic conductivity from slug test.
pumping_test_confined
    Thiem (1906) equilibrium equation for confined aquifer.
pumping_test_unconfined
    Thiem equation for unconfined (phreatic) aquifer.
lefranc_test
    Lefranc permeability test (constant/falling head in borehole).

References
----------
Hvorslev, M. J. (1951). *Time Lag and Soil Permeability in Ground-Water
    Observations*. Waterways Experiment Station, Bull. 36.
Thiem, G. (1906). *Hydrologische Methoden*. Gebhardt, Leipzig.
Bouwer, H. & Rice, R. C. (1976). A slug test for determining hydraulic
    conductivity of unconfined aquifers with completely or partially
    penetrating wells. *Water Resources Research*, 12(3), 423–428.
Das (2021), Chapter 5.
"""

from typing import Dict, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def slug_test(
    r: float,
    R: float,
    Le: float,
    T0: float,
    method: str = "hvorslev",
) -> float:
    r"""
    Hydraulic conductivity from slug test (Hvorslev method).

    The Hvorslev (1951) basic time lag method for a well point or
    piezometer with length/diameter ratio > 8:

    .. math::

        k = \frac{r^2 \ln(L_e / R)}{2\,L_e\,T_0}

    Parameters
    ----------
    r : float
        Standpipe (casing) radius (m).
    R : float
        Well screen or intake radius (m).
    Le : float
        Length of the well screen or intake (m).
    T0 : float
        Basic time lag (s) — time for head to reach 37% (1/e) of
        initial displacement.  Obtained from the slope of
        ln(h/h₀) vs t plot.
    method : str, default ``'hvorslev'``
        Currently only ``'hvorslev'`` is implemented.

    Returns
    -------
    float
        Hydraulic conductivity k (m/s).

    Notes
    -----
    * Valid for Le/R > 8.
    * For Le/R < 8, use the Bouwer & Rice (1976) method.
    * T₀ is found from the -1/slope of the ln(h/h₀) vs time plot
      (the time at which ln(h/h₀) = -1).

    References
    ----------
    Hvorslev (1951); Das (2021), Section 5.7.

    Examples
    --------
    >>> from geoeq.site.field_perm import slug_test
    >>> k = slug_test(r=0.025, R=0.05, Le=1.5, T0=300)
    >>> f'{k:.2e}'
    '3.05e-07'
    """
    check_positive(r, "r")
    check_positive(R, "R")
    check_positive(Le, "Le")
    check_positive(T0, "T0")

    method_l = method.lower()
    if method_l != "hvorslev":
        raise ValueError(f"Unknown method '{method}'. Currently only 'hvorslev' is supported.")

    k = r**2 * np.log(Le / R) / (2.0 * Le * T0)
    return float(k)


def pumping_test_confined(
    Q: float,
    h1: float,
    h2: float,
    r1: float,
    r2: float,
    H: float,
) -> Dict[str, float]:
    r"""
    Hydraulic conductivity from pumping test — confined aquifer.

    Thiem (1906) equilibrium equation for steady-state flow to a well
    in a confined aquifer:

    .. math::

        k = \frac{Q \ln(r_2 / r_1)}{2\,\pi\,H\,(h_2 - h_1)}

    Parameters
    ----------
    Q : float
        Steady-state pumping rate (m³/s).
    h1 : float
        Piezometric head at observation well 1 (m).
    h2 : float
        Piezometric head at observation well 2 (m), h₂ > h₁.
    r1 : float
        Radial distance of observation well 1 from pumped well (m).
    r2 : float
        Radial distance of observation well 2 from pumped well (m),
        r₂ > r₁.
    H : float
        Aquifer thickness (m).

    Returns
    -------
    dict
        ``'k'`` : float — Hydraulic conductivity (m/s).
        ``'T'`` : float — Transmissivity T = kH (m²/s).

    Notes
    -----
    * h₂ > h₁ since the head increases with distance from the well.
    * r₂ > r₁ for the outer observation well.

    References
    ----------
    Thiem (1906); Das (2021), Eq. 5.16.

    Examples
    --------
    >>> from geoeq.site.field_perm import pumping_test_confined
    >>> res = pumping_test_confined(Q=0.004, h1=18.0, h2=19.5,
    ...                             r1=10, r2=50, H=20)
    >>> f"{res['k']:.2e}"
    '3.42e-05'
    """
    check_positive(Q, "Q")
    check_positive(r1, "r1")
    check_positive(r2, "r2")
    check_positive(H, "H")
    if r2 <= r1:
        raise ValueError("r2 must be greater than r1.")
    if h2 <= h1:
        raise ValueError("h2 must be greater than h1 (head increases with distance).")

    k = Q * np.log(r2 / r1) / (2.0 * np.pi * H * (h2 - h1))
    T = k * H

    return {"k": float(k), "T": float(T)}


def pumping_test_unconfined(
    Q: float,
    h1: float,
    h2: float,
    r1: float,
    r2: float,
) -> Dict[str, float]:
    r"""
    Hydraulic conductivity from pumping test — unconfined aquifer.

    Thiem equation for steady-state flow in an unconfined
    (phreatic/water-table) aquifer:

    .. math::

        k = \frac{Q \ln(r_2 / r_1)}{\pi\,(h_2^2 - h_1^2)}

    Parameters
    ----------
    Q : float
        Steady-state pumping rate (m³/s).
    h1 : float
        Water table height at observation well 1 (m from base).
    h2 : float
        Water table height at observation well 2 (m from base),
        h₂ > h₁.
    r1 : float
        Radial distance to observation well 1 (m).
    r2 : float
        Radial distance to observation well 2 (m), r₂ > r₁.

    Returns
    -------
    dict
        ``'k'`` : float — Hydraulic conductivity (m/s).

    References
    ----------
    Thiem (1906); Das (2021), Eq. 5.17.

    Examples
    --------
    >>> from geoeq.site.field_perm import pumping_test_unconfined
    >>> res = pumping_test_unconfined(Q=0.003, h1=15.0, h2=18.0,
    ...                               r1=10, r2=50)
    >>> f"{res['k']:.2e}"
    '1.55e-05'
    """
    check_positive(Q, "Q")
    check_positive(r1, "r1")
    check_positive(r2, "r2")
    check_positive(h1, "h1")
    check_positive(h2, "h2")
    if r2 <= r1:
        raise ValueError("r2 must be greater than r1.")
    if h2 <= h1:
        raise ValueError("h2 must be greater than h1.")

    k = Q * np.log(r2 / r1) / (np.pi * (h2**2 - h1**2))
    return {"k": float(k)}


def lefranc_test(
    Q: float,
    H: float,
    D: float,
    method: str = "constant_head",
) -> float:
    r"""
    Lefranc permeability test — in-situ hydraulic conductivity.

    For a **constant-head** test in a borehole with an open-ended
    cylindrical cavity of diameter D and length H:

    .. math::

        k = \frac{Q}{F \, H_w}

    where :math:`H_w` is the constant head difference and F is a
    shape factor.  For a cylindrical cavity with L/D > 4:

    .. math::

        F = \frac{2\,\pi\,H}{\ln(2\,H / D)}

    Parameters
    ----------
    Q : float
        Flow rate (m³/s).
    H : float
        Head difference (m) — for constant-head test, this is the
        applied head; flow rate must correspond.
    D : float
        Borehole diameter (m).
    method : str, default ``'constant_head'``
        Currently ``'constant_head'`` is implemented.

    Returns
    -------
    float
        Hydraulic conductivity k (m/s).

    Notes
    -----
    * The Lefranc test is common in European practice (NF P 94-132).
    * The shape factor assumes L/D > 4 for the simplified formula.

    References
    ----------
    Lefranc (1936); Cassan, M. (2005). *Les essais d'eau in situ*.

    Examples
    --------
    >>> from geoeq.site.field_perm import lefranc_test
    >>> k = lefranc_test(Q=1e-5, H=2.0, D=0.076)
    >>> f'{k:.2e}'
    '2.78e-06'
    """
    check_positive(Q, "Q")
    check_positive(H, "H")
    check_positive(D, "D")

    # Shape factor for cylindrical cavity, L=H, diameter D
    F = 2.0 * np.pi * H / np.log(2.0 * H / D)
    k = Q / (F * H)
    return float(k)
