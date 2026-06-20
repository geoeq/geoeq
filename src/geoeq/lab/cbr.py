r"""
California Bearing Ratio (CBR) test processing for GeoEq.

Functions
---------
cbr_test
    Process CBR test data → CBR value (ASTM D1883).
cbr_plot
    Plot load–penetration curve with standard crushed-rock overlay.

References
----------
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 5.
ASTM D1883 (2016). Standard Test Method for California Bearing Ratio.
"""

from typing import Dict, List, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


# Standard loads for crushed rock at 2.54 mm and 5.08 mm penetration (kN)
_STANDARD_LOAD_2_54 = 13.24  # 1000 psi × area
_STANDARD_LOAD_5_08 = 19.96  # 1500 psi × area


def cbr_test(
    penetration: Union[List[float], np.ndarray],
    load: Union[List[float], np.ndarray],
    area: float = 19.35,
) -> Dict[str, float]:
    r"""
    Process CBR test data — compute the California Bearing Ratio.

    .. math::

        CBR\,(\%) = \frac{P_{\text{test}}}{P_{\text{standard}}} \times 100

    where :math:`P_{\text{test}}` is the load at 2.54 mm or 5.08 mm
    penetration and :math:`P_{\text{standard}}` is the standard load
    for crushed rock at the same penetration.

    The CBR is reported as the **larger** of the values at 2.54 mm and
    5.08 mm penetration.

    Parameters
    ----------
    penetration : array_like
        Piston penetration values (mm).
    load : array_like
        Corresponding load values (kN).
    area : float, default 19.35
        Piston area (cm²). Standard CBR piston diameter = 49.63 mm.

    Returns
    -------
    dict
        ``'CBR'`` : float — California Bearing Ratio (%).
        ``'CBR_2_54'`` : float — CBR at 2.54 mm penetration (%).
        ``'CBR_5_08'`` : float — CBR at 5.08 mm penetration (%).
        ``'load_2_54'`` : float — test load at 2.54 mm (kN).
        ``'load_5_08'`` : float — test load at 5.08 mm (kN).

    References
    ----------
    Das (2021), Section 5.9; ASTM D1883.

    Examples
    --------
    >>> from geoeq.lab.cbr import cbr_test
    >>> pen = [0, 0.64, 1.27, 2.54, 3.81, 5.08, 7.62, 10.16, 12.70]
    >>> load = [0, 0.8, 1.9, 4.2, 6.1, 7.8, 10.5, 12.1, 13.2]
    >>> res = cbr_test(pen, load)
    >>> res['CBR'] > 0
    True
    """
    p = np.asarray(penetration, dtype=float)
    f = np.asarray(load, dtype=float)

    if len(p) != len(f):
        raise ValueError("penetration and load must have the same length.")
    if len(p) < 3:
        raise ValueError("Need at least 3 data points.")

    idx = np.argsort(p)
    p = p[idx]
    f = f[idx]

    from scipy.interpolate import interp1d
    interp = interp1d(p, f, kind="linear", fill_value="extrapolate")

    load_2_54 = float(interp(2.54))
    load_5_08 = float(interp(5.08))

    cbr_2_54 = (load_2_54 / _STANDARD_LOAD_2_54) * 100.0
    cbr_5_08 = (load_5_08 / _STANDARD_LOAD_5_08) * 100.0

    return {
        "CBR": max(cbr_2_54, cbr_5_08),
        "CBR_2_54": cbr_2_54,
        "CBR_5_08": cbr_5_08,
        "load_2_54": load_2_54,
        "load_5_08": load_5_08,
    }


def cbr_plot(
    penetration: Union[List[float], np.ndarray],
    load: Union[List[float], np.ndarray],
    area: float = 19.35,
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot CBR load–penetration curve with standard reference points.

    Parameters
    ----------
    penetration : array_like
        Penetration values (mm).
    load : array_like
        Load values (kN).
    area : float, default 19.35
        Piston area (cm²).
    ax : matplotlib.axes.Axes, optional
    save_as : str, optional

    Returns
    -------
    dict
        ``'CBR'``, ``'ax'``, etc.

    Examples
    --------
    >>> from geoeq.lab.cbr import cbr_plot
    >>> pen = [0, 0.64, 1.27, 2.54, 3.81, 5.08, 7.62, 10.16, 12.70]
    >>> load = [0, 0.8, 1.9, 4.2, 6.1, 7.8, 10.5, 12.1, 13.2]
    >>> res = cbr_plot(pen, load)
    """
    import matplotlib.pyplot as plt

    p = np.asarray(penetration, dtype=float)
    f = np.asarray(load, dtype=float)

    result = cbr_test(p, f, area)

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 6))
    else:
        fig = ax.get_figure()

    idx = np.argsort(p)
    ax.plot(p[idx], f[idx], "bo-", markersize=6, linewidth=1.5, label="Test data")

    # Reference lines at 2.54 and 5.08 mm
    for pen_val, load_val, cbr_val, label in [
        (2.54, result["load_2_54"], result["CBR_2_54"], "2.54 mm"),
        (5.08, result["load_5_08"], result["CBR_5_08"], "5.08 mm"),
    ]:
        ax.axvline(pen_val, color="red", linestyle="--", alpha=0.5)
        ax.plot(pen_val, load_val, "r*", markersize=12)
        ax.annotate(f"CBR={cbr_val:.1f}%",
                    xy=(pen_val, load_val),
                    xytext=(pen_val + 0.5, load_val + 0.5),
                    fontsize=9, color="red",
                    arrowprops=dict(arrowstyle="->", color="red", alpha=0.5))

    # Standard reference loads
    ax.plot(2.54, _STANDARD_LOAD_2_54, "gs", markersize=8, label="Standard (2.54 mm)")
    ax.plot(5.08, _STANDARD_LOAD_5_08, "g^", markersize=8, label="Standard (5.08 mm)")

    ax.set_xlabel("Penetration (mm)", fontweight="bold")
    ax.set_ylabel("Load (kN)", fontweight="bold")
    ax.set_title(f"CBR Test — CBR = {result['CBR']:.1f}%", fontweight="bold")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    result["ax"] = ax
    return result
