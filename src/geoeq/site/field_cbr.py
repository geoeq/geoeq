r"""
Field CBR and Dynamic Cone Penetrometer (DCP) correlations.

Functions
---------
dcp_cbr
    CBR from DCP penetration rate.
field_cbr_test
    Field (in-situ) CBR from load–penetration data.

References
----------
Webster, S. L., Grau, R. H. & Williams, T. P. (1992). *Description and
    Application of Dual Mass Dynamic Cone Penetrometer*. Report
    GL-92-3, US Army Waterways Experiment Station.
TRL (Transport Research Laboratory) (1993). *A guide to the structural
    design of bitumen-surfaced roads in tropical and sub-tropical
    countries*. Overseas Road Note 31, 4th ed.
ASTM D6951/D6951M (2018). Standard Test Method for Use of the Dynamic
    Cone Penetrometer in Shallow Pavement Applications.
"""

from typing import Dict, Union
import numpy as np
from geoeq.core.validation import check_positive


def dcp_cbr(
    dcpi: Union[float, np.ndarray],
    method: str = "webster",
) -> Union[float, np.ndarray]:
    r"""
    Estimate CBR from Dynamic Cone Penetrometer Index (DCPI).

    The DCP measures penetration per blow (mm/blow).  Several
    empirical correlations relate DCPI to CBR:

    **Webster et al. (1992)** — US Army Corps of Engineers:

    .. math::

        \log_{10}(\text{CBR}) = 1.12 - 0.39 \,\log_{10}(\text{DCPI})

    i.e.  :math:`\text{CBR} = 10^{1.12 - 0.39 \log(\text{DCPI})}`

    **TRL (1993)** — Transport Research Laboratory:

    .. math::

        \log_{10}(\text{CBR}) = 2.48 - 1.057 \,\log_{10}(\text{DCPI})

    **Kleyn (1975)** — South African method:

    .. math::

        \log_{10}(\text{CBR}) = 2.62 - 1.27 \,\log_{10}(\text{DCPI})

    Parameters
    ----------
    dcpi : float or array_like
        DCP Index — penetration per blow (mm/blow).
    method : str, default ``'webster'``
        Correlation method: ``'webster'``, ``'trl'``, or ``'kleyn'``.

    Returns
    -------
    float or ndarray
        Estimated CBR (%).

    References
    ----------
    Webster et al. (1992); TRL (1993); Kleyn (1975).

    Examples
    --------
    >>> from geoeq.site.field_cbr import dcp_cbr
    >>> round(dcp_cbr(dcpi=10, method='webster'), 1)
    5.4
    >>> round(dcp_cbr(dcpi=10, method='trl'), 1)
    26.5
    """
    d = np.asarray(dcpi, dtype=float)
    check_positive(d, "dcpi")
    method_l = method.lower()

    if method_l == "webster":
        cbr = 10 ** (1.12 - 0.39 * np.log10(d))
    elif method_l == "trl":
        cbr = 10 ** (2.48 - 1.057 * np.log10(d))
    elif method_l == "kleyn":
        cbr = 10 ** (2.62 - 1.27 * np.log10(d))
    else:
        raise ValueError(
            f"Unknown method '{method}'. Choose: webster, trl, kleyn."
        )

    if all(np.ndim(x) == 0 for x in [dcpi]):
        return float(cbr)
    return np.asarray(cbr)


def field_cbr_test(
    penetration: Union[list, np.ndarray],
    load: Union[list, np.ndarray],
    area: float = 1935.5,
) -> Dict[str, float]:
    r"""
    Field (in-situ) CBR from load–penetration data.

    The procedure is identical to the laboratory CBR (ASTM D1883)
    but performed on the in-situ subgrade.

    .. math::

        \text{CBR} = \frac{\text{Test load at standard penetration}}
                          {\text{Standard load}} \times 100

    Standard reference loads (for 1935.5 mm² = 3 in² piston):

    - At 2.54 mm (0.1 in.): 13.24 kN
    - At 5.08 mm (0.2 in.): 19.96 kN

    Parameters
    ----------
    penetration : array_like
        Penetration values (mm).
    load : array_like
        Corresponding loads (kN).
    area : float, default 1935.5
        Piston area (mm²).  Standard 3 in² = 1935.5 mm².

    Returns
    -------
    dict
        ``'CBR'`` : float — Governing CBR (%).
        ``'CBR_2_54'`` : float — CBR at 2.54 mm.
        ``'CBR_5_08'`` : float — CBR at 5.08 mm.
        ``'load_2_54'`` : float — Load at 2.54 mm (kN).
        ``'load_5_08'`` : float — Load at 5.08 mm (kN).

    References
    ----------
    ASTM D4429 (2009). Standard Test Method for CBR (California
    Bearing Ratio) of Soils in Place.

    Examples
    --------
    >>> from geoeq.site.field_cbr import field_cbr_test
    >>> pen = [0, 0.64, 1.27, 1.91, 2.54, 3.81, 5.08, 7.62, 10.16, 12.70]
    >>> load = [0, 0.5, 1.2, 2.5, 4.0, 6.5, 9.0, 13.0, 16.0, 18.0]
    >>> res = field_cbr_test(pen, load)
    >>> res['CBR'] > 0
    True
    """
    pen = np.asarray(penetration, dtype=float)
    ld = np.asarray(load, dtype=float)
    if len(pen) != len(ld):
        raise ValueError("penetration and load must have the same length.")
    check_positive(area, "area")

    # Standard reference loads (kN) for 1935.5 mm² piston
    std_load_2_54 = 13.24
    std_load_5_08 = 19.96

    # Interpolate loads at standard penetrations
    load_2_54 = float(np.interp(2.54, pen, ld))
    load_5_08 = float(np.interp(5.08, pen, ld))

    CBR_2_54 = (load_2_54 / std_load_2_54) * 100.0
    CBR_5_08 = (load_5_08 / std_load_5_08) * 100.0

    CBR = max(CBR_2_54, CBR_5_08)

    return {
        "CBR": float(CBR),
        "CBR_2_54": float(CBR_2_54),
        "CBR_5_08": float(CBR_5_08),
        "load_2_54": load_2_54,
        "load_5_08": load_5_08,
    }
