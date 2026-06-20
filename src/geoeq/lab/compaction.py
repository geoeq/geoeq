r"""
Compaction test processing for GeoEq.

Functions
---------
proctor
    Fit Standard or Modified Proctor curve → γ_d_max, w_opt.
zav_line
    Zero Air Voids (ZAV) line for a given Gs.
relative_compaction
    Relative compaction RC = γ_d / γ_d_max.
proctor_plot
    Plot compaction curve with ZAV line.

References
----------
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 5.
ASTM D698 (2012). Standard Proctor Test.
ASTM D1557 (2012). Modified Proctor Test.
Holtz, R. D. et al. (2011). *An Introduction to Geotechnical Engineering*, 2nd ed., Ch. 6.
"""

from typing import Dict, List, Optional, Union, Tuple
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative
from geoeq.core.constants import GAMMA_WATER


def proctor(
    water_content: Union[List[float], np.ndarray],
    dry_density: Union[List[float], np.ndarray],
) -> Dict[str, float]:
    r"""
    Process Proctor compaction test data to find optimum moisture content
    and maximum dry unit weight.

    Fits a second-degree polynomial to the (w, γ_d) data and locates
    the peak.

    Parameters
    ----------
    water_content : array_like
        Water content for each compaction point (%).
    dry_density : array_like
        Dry unit weight for each compaction point (kN/m³).

    Returns
    -------
    dict
        ``'w_opt'`` : float — optimum moisture content (%).
        ``'gamma_d_max'`` : float — maximum dry unit weight (kN/m³).
        ``'coeffs'`` : ndarray — polynomial coefficients [a, b, c].

    References
    ----------
    Das (2021), Section 5.6; ASTM D698 / D1557.

    Examples
    --------
    >>> from geoeq.lab.compaction import proctor
    >>> w = [8, 10, 12, 14, 16, 18]
    >>> gd = [17.5, 18.2, 18.8, 19.0, 18.7, 18.1]
    >>> res = proctor(w, gd)
    >>> 13 < res['w_opt'] < 15
    True
    >>> 18.5 < res['gamma_d_max'] < 19.5
    True
    """
    w = np.asarray(water_content, dtype=float)
    gd = np.asarray(dry_density, dtype=float)

    if len(w) != len(gd):
        raise ValueError("water_content and dry_density must have the same length.")
    if len(w) < 3:
        raise ValueError("Need at least 3 data points for Proctor test.")
    for v in w:
        check_non_negative(v, "water_content")
    for v in gd:
        check_positive(v, "dry_density")

    coeffs = np.polyfit(w, gd, 2)
    a, b, c = coeffs

    if a >= 0:
        # Parabola opens upward — no maximum; return the highest measured point
        i_max = np.argmax(gd)
        return {
            "w_opt": float(w[i_max]),
            "gamma_d_max": float(gd[i_max]),
            "coeffs": coeffs,
        }

    w_opt = -b / (2.0 * a)
    gamma_d_max = float(np.polyval(coeffs, w_opt))

    return {
        "w_opt": float(w_opt),
        "gamma_d_max": gamma_d_max,
        "coeffs": coeffs,
    }


def zav_line(
    Gs: float,
    w_range: Union[List[float], np.ndarray, None] = None,
    gamma_w: float = GAMMA_WATER,
) -> Dict[str, np.ndarray]:
    r"""
    Compute the Zero Air Voids (ZAV) line — the theoretical maximum dry
    unit weight at S = 100 %.

    .. math::

        \gamma_{d,\text{zav}} = \frac{G_s \, \gamma_w}{1 + w \, G_s}
        \qquad \text{[Das Eq.\,5.12]}

    Parameters
    ----------
    Gs : float
        Specific gravity of soil solids (typically 2.60–2.80).
    w_range : array_like, optional
        Water content values (%) at which to compute γ_d. Default
        ``[4, 6, 8, …, 30]``.
    gamma_w : float, default 9.81
        Unit weight of water (kN/m³).

    Returns
    -------
    dict
        ``'water_content'`` : ndarray — w values (%).
        ``'dry_density'`` : ndarray — γ_d at S = 100 % (kN/m³).

    References
    ----------
    Das (2021), Eq. 5.12.

    Examples
    --------
    >>> from geoeq.lab.compaction import zav_line
    >>> res = zav_line(Gs=2.70)
    >>> res['dry_density'][0] > res['dry_density'][-1]
    True
    """
    check_positive(Gs, "Gs")
    if w_range is None:
        w_range = np.arange(4, 32, 2, dtype=float)
    else:
        w_range = np.asarray(w_range, dtype=float)

    w_dec = w_range / 100.0
    gd = (Gs * gamma_w) / (1.0 + w_dec * Gs)

    return {
        "water_content": w_range,
        "dry_density": gd,
    }


def saturation_line(
    Gs: float,
    S: float,
    w_range: Union[List[float], np.ndarray, None] = None,
    gamma_w: float = GAMMA_WATER,
) -> Dict[str, np.ndarray]:
    r"""
    Compute a constant-saturation line on the compaction curve.

    .. math::

        \gamma_d = \frac{G_s \, \gamma_w}{1 + \frac{w \, G_s}{S}}

    where *S* is a fraction (0–1).

    Parameters
    ----------
    Gs : float
        Specific gravity.
    S : float
        Degree of saturation as a fraction (e.g. 0.8 for 80 %).
    w_range : array_like, optional
        Water content values (%). Default ``[4, 6, …, 30]``.
    gamma_w : float, default 9.81

    Returns
    -------
    dict
        ``'water_content'`` : ndarray (%).
        ``'dry_density'`` : ndarray (kN/m³).

    Examples
    --------
    >>> from geoeq.lab.compaction import saturation_line
    >>> res = saturation_line(Gs=2.70, S=0.8)
    >>> all(res['dry_density'] > 0)
    True
    """
    check_positive(Gs, "Gs")
    check_positive(S, "S")
    if S > 1.0:
        raise ValueError(f"S must be ≤ 1.0, got {S}.")

    if w_range is None:
        w_range = np.arange(4, 32, 2, dtype=float)
    else:
        w_range = np.asarray(w_range, dtype=float)

    w_dec = w_range / 100.0
    gd = (Gs * gamma_w) / (1.0 + w_dec * Gs / S)

    return {
        "water_content": w_range,
        "dry_density": gd,
    }


def relative_compaction(
    gamma_d: Union[float, np.ndarray],
    gamma_d_max: float,
) -> Union[float, np.ndarray]:
    r"""
    Compute relative compaction.

    .. math::

        RC\,(\%) = \frac{\gamma_d}{\gamma_{d,\max}} \times 100

    Parameters
    ----------
    gamma_d : float or array_like
        Field dry unit weight (kN/m³).
    gamma_d_max : float
        Maximum dry unit weight from Proctor test (kN/m³).

    Returns
    -------
    float or ndarray
        Relative compaction (%).

    Examples
    --------
    >>> from geoeq.lab.compaction import relative_compaction
    >>> relative_compaction(17.5, 19.0)
    92.1...
    """
    gd = np.asarray(gamma_d, dtype=float)
    check_positive(gamma_d_max, "gamma_d_max")
    check_positive(gd, "gamma_d")

    rc = (gd / gamma_d_max) * 100.0

    if np.ndim(gamma_d) == 0:
        return float(rc)
    return rc


def proctor_plot(
    water_content: Union[List[float], np.ndarray],
    dry_density: Union[List[float], np.ndarray],
    Gs: float = 2.65,
    show_zav: bool = True,
    show_sat_lines: bool = True,
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot Proctor compaction curve with ZAV line and saturation contours.

    Parameters
    ----------
    water_content : array_like
        Water content (%).
    dry_density : array_like
        Dry unit weight (kN/m³).
    Gs : float, default 2.65
        Specific gravity for ZAV and saturation lines.
    show_zav : bool, default True
        Plot the zero air voids line.
    show_sat_lines : bool, default True
        Plot 60%, 80% saturation contours.
    ax : matplotlib.axes.Axes, optional
    save_as : str, optional

    Returns
    -------
    dict
        ``'w_opt'``, ``'gamma_d_max'``, ``'ax'``.

    Examples
    --------
    >>> from geoeq.lab.compaction import proctor_plot
    >>> w = [8, 10, 12, 14, 16, 18]
    >>> gd = [17.5, 18.2, 18.8, 19.0, 18.7, 18.1]
    >>> res = proctor_plot(w, gd, Gs=2.70)
    """
    import matplotlib.pyplot as plt

    w = np.asarray(water_content, dtype=float)
    gd = np.asarray(dry_density, dtype=float)

    result = proctor(w, gd)

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 6))
    else:
        fig = ax.get_figure()

    # Plot data points
    ax.scatter(w, gd, s=80, color="steelblue", edgecolors="navy", zorder=5,
               label="Test data")

    # Fit curve
    w_fine = np.linspace(float(w.min()) - 1, float(w.max()) + 1, 200)
    gd_fine = np.polyval(result["coeffs"], w_fine)
    ax.plot(w_fine, gd_fine, "b-", linewidth=2, label="Best-fit curve")

    # Optimum point
    ax.axvline(result["w_opt"], color="red", linestyle="--", alpha=0.5)
    ax.axhline(result["gamma_d_max"], color="red", linestyle="--", alpha=0.5)
    ax.plot(result["w_opt"], result["gamma_d_max"], "r*", markersize=15,
            label=f"Optimum: w={result['w_opt']:.1f}%, γ_d={result['gamma_d_max']:.2f}")

    # ZAV line
    w_zav = np.linspace(float(w.min()) - 2, float(w.max()) + 4, 100)
    if show_zav:
        zav = zav_line(Gs, w_zav)
        ax.plot(zav["water_content"], zav["dry_density"], "k-", linewidth=2,
                label=f"ZAV (Gs={Gs})")

    # Saturation lines
    if show_sat_lines:
        for S_pct in [60, 80]:
            sat = saturation_line(Gs, S_pct / 100.0, w_zav)
            ax.plot(sat["water_content"], sat["dry_density"], "--",
                    color="gray", linewidth=1, alpha=0.7)
            # Label at the end
            ax.annotate(f"S={S_pct}%",
                        xy=(float(sat["water_content"][-1]),
                            float(sat["dry_density"][-1])),
                        fontsize=8, color="gray")

    ax.set_xlabel("Water Content w (%)", fontweight="bold")
    ax.set_ylabel("Dry Unit Weight γ_d (kN/m³)", fontweight="bold")
    ax.set_title("Proctor Compaction Test", fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    # Reasonable y-limits
    margin = 0.5
    ax.set_ylim(float(gd.min()) - margin, float(max(gd.max(), result["gamma_d_max"])) + margin)

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    result["ax"] = ax
    return result
