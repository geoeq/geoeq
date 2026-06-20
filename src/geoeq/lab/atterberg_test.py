r"""
Atterberg limits laboratory test processing for GeoEq.

Functions
---------
liquid_limit_test
    Process Casagrande cup data → LL from the flow curve (ASTM D4318).
flow_curve_plot
    Plot the flow curve (w vs N) and determine LL.

References
----------
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 4.
ASTM D4318 (2017). Standard Test Methods for Liquid Limit, Plastic Limit.
Casagrande, A. (1932). Research on the Atterberg limits of soils.
"""

from typing import Dict, List, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive


def liquid_limit_test(
    blow_count: Union[List[float], np.ndarray],
    water_content: Union[List[float], np.ndarray],
) -> Dict[str, float]:
    r"""
    Determine liquid limit from Casagrande cup test data.

    The liquid limit is the water content at 25 blows on the
    semi-logarithmic flow curve:

    .. math::

        w = a \cdot N^b

    where *N* is the blow count and *a*, *b* are fitted constants.
    LL is evaluated at N = 25.

    Parameters
    ----------
    blow_count : array_like
        Number of blows for each trial (typically 3–5 trials spanning
        15–35 blows).
    water_content : array_like
        Water content for each trial (%).

    Returns
    -------
    dict
        ``'LL'`` : float — liquid limit (%).
        ``'slope'`` : float — flow index I_f (slope of flow curve).
        ``'r_squared'`` : float — goodness of fit.

    References
    ----------
    Das (2021), Section 4.4; ASTM D4318.

    Examples
    --------
    >>> from geoeq.lab.atterberg_test import liquid_limit_test
    >>> res = liquid_limit_test([15, 20, 28, 34], [42.1, 40.8, 38.5, 36.9])
    >>> 38 < res['LL'] < 40
    True
    """
    N = np.asarray(blow_count, dtype=float)
    w = np.asarray(water_content, dtype=float)

    if len(N) < 2:
        raise ValueError("Need at least 2 data points.")
    if len(N) != len(w):
        raise ValueError("blow_count and water_content must have the same length.")
    for v in N:
        check_positive(v, "blow_count")
    for v in w:
        check_positive(v, "water_content")

    # Semi-log fit: w vs log10(N)
    log_N = np.log10(N)
    coeffs = np.polyfit(log_N, w, 1)
    slope = coeffs[0]  # flow index
    intercept = coeffs[1]

    LL = float(np.polyval(coeffs, np.log10(25.0)))

    # R²
    w_pred = np.polyval(coeffs, log_N)
    ss_res = np.sum((w - w_pred) ** 2)
    ss_tot = np.sum((w - np.mean(w)) ** 2)
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "LL": LL,
        "slope": float(slope),
        "r_squared": float(r_sq),
    }


def flow_curve_plot(
    blow_count: Union[List[float], np.ndarray],
    water_content: Union[List[float], np.ndarray],
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot the flow curve (water content vs blow count on semi-log scale).

    Parameters
    ----------
    blow_count : array_like
        Number of blows.
    water_content : array_like
        Water content (%).
    ax : matplotlib.axes.Axes, optional
    save_as : str, optional

    Returns
    -------
    dict
        ``'LL'``, ``'slope'``, ``'r_squared'``, ``'ax'``.

    Examples
    --------
    >>> from geoeq.lab.atterberg_test import flow_curve_plot
    >>> res = flow_curve_plot([15, 20, 28, 34], [42.1, 40.8, 38.5, 36.9])
    """
    import matplotlib.pyplot as plt

    N = np.asarray(blow_count, dtype=float)
    w = np.asarray(water_content, dtype=float)

    result = liquid_limit_test(N, w)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    ax.scatter(N, w, s=80, color="steelblue", edgecolors="navy", zorder=5,
               label="Test data")

    N_fit = np.linspace(10, 40, 100)
    log_N_fit = np.log10(N_fit)
    w_fit = result["slope"] * log_N_fit + (result["LL"] - result["slope"] * np.log10(25.0))
    ax.plot(N_fit, w_fit, "r-", linewidth=2, label="Flow curve")

    ax.axvline(25, color="green", linestyle="--", alpha=0.7, label="N = 25")
    ax.plot(25, result["LL"], "r*", markersize=15,
            label=f"LL = {result['LL']:.1f}%")

    ax.set_xscale("log")
    ax.set_xlabel("Number of Blows N", fontweight="bold")
    ax.set_ylabel("Water Content w (%)", fontweight="bold")
    ax.set_title("Flow Curve — Liquid Limit Determination", fontweight="bold")
    ax.set_xlim(10, 100)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    ax.text(0.97, 0.97,
            f"LL = {result['LL']:.1f}%\n$I_f$ = {result['slope']:.2f}\n$R^2$ = {result['r_squared']:.4f}",
            transform=ax.transAxes, fontsize=10, va="top", ha="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="gray"))

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    result["ax"] = ax
    return result
