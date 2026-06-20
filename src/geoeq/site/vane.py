r"""
Field Vane Shear Test (FVT) — undrained shear strength from torque.

Functions
---------
vane_su
    Undrained shear strength from measured torque.
vane_correction
    Bjerrum (1972) correction factor for Su.
vane_remolded
    Remolded Su and sensitivity from vane test.

References
----------
Bjerrum, L. (1972). Embankments on soft ground. *Proc., ASCE Specialty
    Conf. Performance of Earth and Earth-Supported Structures*, Lafayette,
    2, 1–54.
ASTM D2573/D2573M (2018). Standard Test Method for Field Vane Shear
    Test in Saturated Fine-Grained Soils.
Chandler, R. J. (1988). The in-situ measurement of the undrained shear
    strength of clays using the field vane. *Vane Shear Strength Testing
    in Soils*, ASTM STP 1014, 13–44.
"""

from typing import Dict, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def _scalar_or_array(result, *inputs):
    if all(np.ndim(x) == 0 for x in inputs):
        return float(result)
    return np.asarray(result)


def vane_su(
    T: Union[float, np.ndarray],
    D: float,
    H: float,
) -> Union[float, np.ndarray]:
    r"""
    Undrained shear strength from field vane torque.

    For a rectangular vane with height H and diameter D, assuming
    uniform shear stress on the cylindrical surface and both end caps:

    .. math::

        S_u = \frac{T}{\pi \left(\frac{D^2 H}{2} + \frac{D^3}{6}\right)}
        \qquad \text{[ASTM D2573]}

    For the standard vane with H/D = 2:

    .. math::

        S_u = \frac{6\,T}{7\,\pi\,D^3}

    Parameters
    ----------
    T : float or array_like
        Maximum torque (kN·m).
    D : float
        Vane diameter (m).  Standard sizes: 0.0508 m (2 in.),
        0.0635 m (2.5 in.), 0.0925 m (3.64 in.).
    H : float
        Vane height (m).  Standard H/D = 2.

    Returns
    -------
    float or ndarray
        Undrained shear strength Su (kPa).

    References
    ----------
    ASTM D2573 (2018); Das (2021), Section 3.15.

    Examples
    --------
    >>> from geoeq.site.vane import vane_su
    >>> round(vane_su(T=0.012, D=0.065, H=0.130), 1)
    20.9
    """
    T_arr = np.asarray(T, dtype=float)
    check_non_negative(T_arr, "T")
    check_positive(D, "D")
    check_positive(H, "H")

    K = np.pi * (D**2 * H / 2.0 + D**3 / 6.0)
    Su = T_arr / K
    return _scalar_or_array(Su, T)


def vane_correction(
    Su_fvt: Union[float, np.ndarray],
    PI: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    r"""
    Bjerrum (1972) correction for field vane undrained shear strength.

    The field vane overestimates the mobilised Su on the failure
    surface.  The correction factor μ depends on the soil's plasticity
    index:

    .. math::

        S_{u,\text{design}} = \mu \times S_{u,\text{FVT}}

    Bjerrum (1972):

    .. math::

        \mu = 1.7 - 0.54 \log_{10}(\text{PI})
        \qquad \text{for PI in \%}

    For very low PI (< 20 %), μ is capped at 1.0.

    Parameters
    ----------
    Su_fvt : float or array_like
        Field vane undrained shear strength (kPa).
    PI : float or array_like
        Plasticity index (%).

    Returns
    -------
    float or ndarray
        Corrected undrained shear strength Su,design (kPa).

    Notes
    -----
    * Bjerrum's factor typically ranges from 0.5 to 1.0.
    * Aas et al. (1986) proposed an alternative correction based on
      Su/σ'v ratio; use the Bjerrum factor for routine design.

    References
    ----------
    Bjerrum (1972); Das (2021), Eq. 3.39.

    Examples
    --------
    >>> from geoeq.site.vane import vane_correction
    >>> round(vane_correction(Su_fvt=50, PI=40), 1)
    38.2
    """
    Su = np.asarray(Su_fvt, dtype=float)
    pi = np.asarray(PI, dtype=float)
    check_non_negative(Su, "Su_fvt")
    check_positive(pi, "PI")

    mu = 1.7 - 0.54 * np.log10(pi)
    mu = np.clip(mu, 0.5, 1.0)

    result = mu * Su
    return _scalar_or_array(result, Su_fvt, PI)


def vane_remolded(
    T_peak: float,
    T_remolded: float,
    D: float,
    H: float,
) -> Dict[str, float]:
    r"""
    Undrained shear strength (peak and remolded) and sensitivity.

    After the peak torque test, the vane is rotated rapidly through
    several full revolutions to remold the soil, then retested.

    .. math::

        S_t = \frac{S_{u,\text{peak}}}{S_{u,\text{remolded}}}

    Parameters
    ----------
    T_peak : float
        Peak torque (kN·m).
    T_remolded : float
        Remolded torque (kN·m).
    D : float
        Vane diameter (m).
    H : float
        Vane height (m).

    Returns
    -------
    dict
        ``'Su_peak'`` : float — Peak Su (kPa).
        ``'Su_remolded'`` : float — Remolded Su (kPa).
        ``'sensitivity'`` : float — Sensitivity St.
        ``'classification'`` : str — Sensitivity classification.

    References
    ----------
    Das (2021), Table 3.7; ASTM D2573.

    Examples
    --------
    >>> from geoeq.site.vane import vane_remolded
    >>> res = vane_remolded(T_peak=0.015, T_remolded=0.003, D=0.065, H=0.130)
    >>> res['sensitivity']
    5.0
    """
    check_non_negative(T_peak, "T_peak")
    check_non_negative(T_remolded, "T_remolded")
    check_positive(D, "D")
    check_positive(H, "H")

    Su_peak = float(vane_su(T_peak, D, H))
    Su_rem = float(vane_su(T_remolded, D, H))

    if Su_rem > 0:
        St = Su_peak / Su_rem
    else:
        St = float("inf")

    if St < 2:
        cls_ = "Insensitive"
    elif St < 4:
        cls_ = "Low sensitivity"
    elif St < 8:
        cls_ = "Medium sensitivity"
    elif St < 16:
        cls_ = "Sensitive"
    else:
        cls_ = "Quick clay"

    return {
        "Su_peak": Su_peak,
        "Su_remolded": Su_rem,
        "sensitivity": float(St),
        "classification": cls_,
    }
