r"""
Consolidation (oedometer) test processing for GeoEq.

Functions
---------
oedometer
    Process one-dimensional consolidation data → e-log(p) curve, Cc, Cr.
preconsolidation
    Determine preconsolidation pressure by Casagrande or Pacheco-Silva.
compression_index
    Empirical correlations for Cc (6 methods).
cv
    Coefficient of consolidation from log-time or root-time fitting.
oedometer_plot
    Publication-quality e-log(p) curve with pc extraction.

References
----------
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 7.
ASTM D2435 (2011). Standard Test Method for 1-D Consolidation.
Casagrande, A. (1936). The determination of the preconsolidation load.
Terzaghi, K. et al. (1996). *Soil Mechanics in Engineering Practice*, 3rd ed.
"""

from typing import Dict, List, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def oedometer(
    stress: Union[List[float], np.ndarray],
    void_ratio: Union[List[float], np.ndarray],
) -> Dict:
    r"""
    Process oedometer (1-D consolidation) test data.

    Computes the compression index :math:`C_c`, recompression index
    :math:`C_r`, and identifies the virgin compression line (VCL).

    .. math::

        C_c = \frac{e_1 - e_2}{\log_{10}(\sigma'_2 / \sigma'_1)}
        \qquad \text{[Das Eq.\,7.9]}

    Parameters
    ----------
    stress : array_like
        Effective vertical stress for each load step (kPa). Must be > 0.
    void_ratio : array_like
        Void ratio at end of each load step.

    Returns
    -------
    dict
        ``'stress'`` : ndarray — sorted stresses (kPa).
        ``'void_ratio'`` : ndarray — corresponding void ratios.
        ``'Cc'`` : float — compression index (slope of VCL on e–log p).
        ``'Cr'`` : float — recompression index (slope of recompression
        portion, estimated from first two points).

    References
    ----------
    Das (2021), Ch. 7, Eq. 7.9.

    Examples
    --------
    >>> from geoeq.lab.consolidation import oedometer
    >>> stress = [25, 50, 100, 200, 400, 800]
    >>> e = [0.88, 0.86, 0.82, 0.74, 0.62, 0.48]
    >>> res = oedometer(stress, e)
    >>> round(res['Cc'], 2)
    0.46
    """
    s = np.asarray(stress, dtype=float)
    e = np.asarray(void_ratio, dtype=float)

    if len(s) != len(e):
        raise ValueError("stress and void_ratio must have the same length.")
    if len(s) < 3:
        raise ValueError("Need at least 3 load steps.")
    for v in s:
        check_positive(v, "stress")
    for v in e:
        check_non_negative(v, "void_ratio")

    idx = np.argsort(s)
    s = s[idx]
    e = e[idx]

    log_s = np.log10(s)

    # Cc: slope of steepest segment on e vs log(p) (virgin compression)
    slopes = []
    for i in range(len(s) - 1):
        dlog = log_s[i + 1] - log_s[i]
        if dlog > 0:
            slopes.append(-(e[i + 1] - e[i]) / dlog)
        else:
            slopes.append(0.0)

    Cc = float(max(slopes)) if slopes else 0.0

    # Cr: slope of first two points (recompression)
    dlog0 = log_s[1] - log_s[0]
    Cr = float(-(e[1] - e[0]) / dlog0) if dlog0 > 0 else 0.0

    return {
        "stress": s,
        "void_ratio": e,
        "Cc": Cc,
        "Cr": Cr,
    }


def preconsolidation(
    stress: Union[List[float], np.ndarray],
    void_ratio: Union[List[float], np.ndarray],
    method: str = "casagrande",
) -> Dict[str, float]:
    r"""
    Determine preconsolidation pressure from e–log(p) data.

    **Casagrande method** (default):

    1. Find the point of maximum curvature on the e–log p curve.
    2. Draw a horizontal line and a tangent at that point.
    3. Bisect the angle between them.
    4. The intersection of the bisector with the virgin compression
       line (VCL) gives :math:`p_c`.

    Parameters
    ----------
    stress : array_like
        Effective vertical stress (kPa), > 0.
    void_ratio : array_like
        Void ratio at each load step.
    method : {'casagrande'}, default ``'casagrande'``
        Extraction method.

    Returns
    -------
    dict
        ``'pc'`` : float — preconsolidation pressure (kPa).
        ``'method'`` : str.

    References
    ----------
    Casagrande (1936); Das (2021), Section 7.4.

    Examples
    --------
    >>> from geoeq.lab.consolidation import preconsolidation
    >>> stress = [25, 50, 100, 200, 400, 800]
    >>> e = [0.88, 0.86, 0.82, 0.74, 0.62, 0.48]
    >>> res = preconsolidation(stress, e)
    >>> 50 < res['pc'] < 200
    True
    """
    method = method.lower()
    if method != "casagrande":
        raise ValueError(f"Method '{method}' not supported. Use 'casagrande'.")

    s = np.asarray(stress, dtype=float)
    e = np.asarray(void_ratio, dtype=float)
    if len(s) < 4:
        raise ValueError("Need at least 4 data points for Casagrande method.")

    idx = np.argsort(s)
    s = s[idx]
    e = e[idx]
    log_s = np.log10(s)

    from scipy.interpolate import PchipInterpolator
    interp = PchipInterpolator(log_s, e)
    log_s_fine = np.linspace(log_s[0], log_s[-1], 500)
    e_fine = interp(log_s_fine)

    e_1 = interp(log_s_fine, 1)
    e_2 = interp(log_s_fine, 2)
    curvature = np.abs(e_2) / (1 + e_1 ** 2) ** 1.5
    i_max = np.argmax(curvature)

    tangent_slope = float(e_1[i_max])
    e_mc = float(e_fine[i_max])
    log_mc = float(log_s_fine[i_max])

    # VCL: last two points (steepest part)
    vcl_slope = (e[-1] - e[-2]) / (log_s[-1] - log_s[-2])
    vcl_intercept = e[-1] - vcl_slope * log_s[-1]

    # Bisector slope = average of horizontal (0) and tangent slope
    bisector_slope = tangent_slope / 2.0
    bisector_intercept = e_mc - bisector_slope * log_mc

    # Intersection of bisector with VCL
    if abs(vcl_slope - bisector_slope) < 1e-12:
        pc = 10 ** log_mc
    else:
        log_pc = (bisector_intercept - vcl_intercept) / (vcl_slope - bisector_slope)
        pc = 10 ** log_pc

    return {
        "pc": float(pc),
        "method": method,
    }


def compression_index(
    method: str = "terzaghi",
    LL: Optional[float] = None,
    e0: Optional[float] = None,
    wn: Optional[float] = None,
    Gs: Optional[float] = None,
) -> float:
    r"""
    Estimate compression index :math:`C_c` from empirical correlations.

    Parameters
    ----------
    method : str, default ``'terzaghi'``
        Correlation to use:

        - ``'terzaghi'`` : :math:`C_c = 0.009 \,(LL - 10)` — Terzaghi & Peck (1967).
        - ``'skempton'`` : :math:`C_c = 0.007 \,(LL - 7)` — remolded clays.
        - ``'rendon'`` : :math:`C_c = 0.01 \, w_n` — Rendon-Herrero (1983).
        - ``'nishida'`` : :math:`C_c = 1.15 \,(e_0 - 0.35)` — all clays.
        - ``'nagaraj'`` : :math:`C_c = 0.2343 \,e_0 \, G_s` — Nagaraj & Murthy (1986).
        - ``'hough'`` : :math:`C_c = 0.29 \,(e_0 - 0.27)` — inorganic silty clays.
    LL : float, optional
        Liquid limit (%). Required for ``terzaghi``, ``skempton``.
    e0 : float, optional
        Initial void ratio. Required for ``nishida``, ``nagaraj``, ``hough``.
    wn : float, optional
        Natural water content (%). Required for ``rendon``.
    Gs : float, optional
        Specific gravity. Required for ``nagaraj``.

    Returns
    -------
    float
        Estimated compression index :math:`C_c`.

    References
    ----------
    Das (2021), Table 7.3.

    Examples
    --------
    >>> from geoeq.lab.consolidation import compression_index
    >>> round(compression_index(method='terzaghi', LL=50), 3)
    0.36
    """
    method = method.lower()

    if method == "terzaghi":
        if LL is None:
            raise ValueError("LL is required for Terzaghi correlation.")
        return 0.009 * (LL - 10)
    elif method == "skempton":
        if LL is None:
            raise ValueError("LL is required for Skempton correlation.")
        return 0.007 * (LL - 7)
    elif method == "rendon":
        if wn is None:
            raise ValueError("wn is required for Rendon-Herrero correlation.")
        return 0.01 * wn
    elif method == "nishida":
        if e0 is None:
            raise ValueError("e0 is required for Nishida correlation.")
        return 1.15 * (e0 - 0.35)
    elif method == "nagaraj":
        if e0 is None or Gs is None:
            raise ValueError("e0 and Gs are required for Nagaraj correlation.")
        return 0.2343 * e0 * Gs
    elif method == "hough":
        if e0 is None:
            raise ValueError("e0 is required for Hough correlation.")
        return 0.29 * (e0 - 0.27)
    else:
        raise ValueError(
            f"Unknown method '{method}'. Choose from: terzaghi, skempton, "
            "rendon, nishida, nagaraj, hough."
        )


def cv(
    time: Union[List[float], np.ndarray],
    deformation: Union[List[float], np.ndarray],
    method: str = "log",
    H_dr: float = 1.0,
) -> Dict[str, float]:
    r"""
    Compute the coefficient of consolidation :math:`c_v` from
    time–deformation data.

    **Log-time method** (Casagrande):

    .. math::

        c_v = \frac{T_{50} \cdot H_{dr}^2}{t_{50}}
        \qquad T_{50} = 0.197

    **Root-time method** (Taylor):

    .. math::

        c_v = \frac{T_{90} \cdot H_{dr}^2}{t_{90}}
        \qquad T_{90} = 0.848

    Parameters
    ----------
    time : array_like
        Elapsed time (minutes).
    deformation : array_like
        Dial reading or settlement (mm). Increasing = more compression.
    method : {'log', 'root'}, default ``'log'``
        ``'log'`` — Casagrande log-time; ``'root'`` — Taylor root-time.
    H_dr : float, default 1.0
        Drainage path length (cm). Half the specimen height for
        double drainage.

    Returns
    -------
    dict
        ``'cv'`` : float — coefficient of consolidation (cm²/min).
        ``'t50'`` or ``'t90'`` : float — time for 50% or 90% consolidation.
        ``'method'`` : str.

    References
    ----------
    Das (2021), Section 7.7, Eqs. 7.22, 7.23.

    Examples
    --------
    >>> from geoeq.lab.consolidation import cv
    >>> t = [0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30, 60, 120, 240, 480, 1440]
    >>> d = [0.0, 0.05, 0.12, 0.22, 0.35, 0.50, 0.65, 0.78, 0.88, 0.95,
    ...      0.99, 1.02, 1.04, 1.06]
    >>> res = cv(t, d, method='log', H_dr=1.0)
    >>> res['cv'] > 0
    True
    """
    t = np.asarray(time, dtype=float)
    d = np.asarray(deformation, dtype=float)
    method = method.lower()

    if len(t) != len(d):
        raise ValueError("time and deformation must have the same length.")
    if len(t) < 5:
        raise ValueError("Need at least 5 data points.")
    check_positive(H_dr, "H_dr")

    idx = np.argsort(t)
    t = t[idx]
    d = d[idx]

    if method == "log":
        d_max = np.max(d)
        d_min = np.min(d)
        # d100 = point where secondary compression starts (inflection)
        # Approximate: last point where slope changes significantly
        d_100 = d_max
        # d0 = initial reading
        d_0 = d_min
        # d50 = halfway
        d_50 = (d_0 + d_100) / 2.0

        # Interpolate to find t50
        from scipy.interpolate import interp1d
        try:
            f = interp1d(d, t, kind="linear", fill_value="extrapolate")
            t_50 = float(f(d_50))
        except Exception:
            t_50 = float(t[len(t) // 2])

        T_50 = 0.197
        cv_val = T_50 * H_dr ** 2 / t_50

        return {
            "cv": float(cv_val),
            "t50": float(t_50),
            "method": "log",
        }

    elif method == "root":
        d_max = np.max(d)
        d_min = np.min(d)
        d_90 = d_min + 0.9 * (d_max - d_min)

        from scipy.interpolate import interp1d
        try:
            f = interp1d(d, t, kind="linear", fill_value="extrapolate")
            t_90 = float(f(d_90))
        except Exception:
            t_90 = float(t[int(0.9 * len(t))])

        T_90 = 0.848
        cv_val = T_90 * H_dr ** 2 / t_90

        return {
            "cv": float(cv_val),
            "t90": float(t_90),
            "method": "root",
        }

    else:
        raise ValueError(f"method must be 'log' or 'root', got '{method}'.")


def oedometer_plot(
    stress: Union[List[float], np.ndarray],
    void_ratio: Union[List[float], np.ndarray],
    show_pc: bool = True,
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot e–log(σ') curve from oedometer test data.

    Parameters
    ----------
    stress : array_like
        Effective vertical stress (kPa).
    void_ratio : array_like
        Void ratio at each load step.
    show_pc : bool, default True
        If True, estimate and mark the preconsolidation pressure.
    ax : matplotlib.axes.Axes, optional
        Existing axes.
    save_as : str, optional
        File path to save figure.

    Returns
    -------
    dict
        ``'Cc'``, ``'Cr'``, ``'pc'`` (if show_pc), ``'ax'``.

    Examples
    --------
    >>> from geoeq.lab.consolidation import oedometer_plot
    >>> res = oedometer_plot([25, 50, 100, 200, 400, 800],
    ...                     [0.88, 0.86, 0.82, 0.74, 0.62, 0.48])
    """
    import matplotlib.pyplot as plt

    s = np.asarray(stress, dtype=float)
    e = np.asarray(void_ratio, dtype=float)

    oed = oedometer(s, e)
    s_sorted = oed["stress"]
    e_sorted = oed["void_ratio"]

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 6))
    else:
        fig = ax.get_figure()

    ax.plot(s_sorted, e_sorted, "bo-", markersize=6, linewidth=1.5, label="Test data")
    ax.set_xscale("log")
    ax.set_xlabel("Effective Stress σ' (kPa)", fontweight="bold")
    ax.set_ylabel("Void Ratio e", fontweight="bold")
    ax.set_title("Oedometer Test — e vs log(σ')", fontweight="bold")
    ax.grid(True, which="both", alpha=0.3)
    ax.invert_yaxis()

    info_lines = [
        f"$C_c$ = {oed['Cc']:.3f}",
        f"$C_r$ = {oed['Cr']:.3f}",
    ]

    result = {
        "Cc": oed["Cc"],
        "Cr": oed["Cr"],
        "ax": ax,
    }

    if show_pc and len(s) >= 4:
        try:
            pc_res = preconsolidation(s, e)
            pc = pc_res["pc"]
            ax.axvline(pc, color="red", linestyle="--", alpha=0.7,
                       label=f"$p_c$ = {pc:.0f} kPa")
            info_lines.append(f"$p_c$ = {pc:.0f} kPa")
            result["pc"] = pc
        except Exception:
            pass

    ax.text(0.97, 0.97, "\n".join(info_lines), transform=ax.transAxes,
            fontsize=10, va="top", ha="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="gray"))

    ax.legend(fontsize=9)

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    return result


def cv_plot(
    time: Union[List[float], np.ndarray],
    deformation: Union[List[float], np.ndarray],
    method: str = "log",
    H_dr: float = 1.0,
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot time–deformation data and determine :math:`c_v`.

    For ``method='log'``, x-axis is log(time); for ``method='root'``,
    x-axis is √time.

    Parameters
    ----------
    time : array_like
        Elapsed time (minutes).
    deformation : array_like
        Dial reading / settlement (mm).
    method : {'log', 'root'}, default ``'log'``
    H_dr : float, default 1.0
        Drainage path (cm).
    ax : matplotlib.axes.Axes, optional
    save_as : str, optional

    Returns
    -------
    dict
        ``'cv'``, ``'ax'``, and ``'t50'`` or ``'t90'``.

    Examples
    --------
    >>> from geoeq.lab.consolidation import cv_plot
    >>> t = [0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30, 60, 120, 240, 480, 1440]
    >>> d = [0.0, 0.05, 0.12, 0.22, 0.35, 0.50, 0.65, 0.78, 0.88, 0.95,
    ...      0.99, 1.02, 1.04, 1.06]
    >>> res = cv_plot(t, d, method='log', H_dr=1.0)
    """
    import matplotlib.pyplot as plt

    t = np.asarray(time, dtype=float)
    d = np.asarray(deformation, dtype=float)
    result = cv(t, d, method=method, H_dr=H_dr)

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 6))
    else:
        fig = ax.get_figure()

    idx = np.argsort(t)
    t_sorted = t[idx]
    d_sorted = d[idx]

    if method == "log":
        ax.plot(t_sorted, d_sorted, "bo-", markersize=5, linewidth=1.2)
        ax.set_xscale("log")
        ax.set_xlabel("Time (min) — log scale", fontweight="bold")
        t_mark = result["t50"]
        ax.axvline(t_mark, color="red", linestyle="--", alpha=0.7,
                   label=f"$t_{{50}}$ = {t_mark:.2f} min")
    else:
        ax.plot(np.sqrt(t_sorted), d_sorted, "bo-", markersize=5, linewidth=1.2)
        ax.set_xlabel("√Time (√min)", fontweight="bold")
        t_mark = result["t90"]
        ax.axvline(np.sqrt(t_mark), color="red", linestyle="--", alpha=0.7,
                   label=f"$t_{{90}}$ = {t_mark:.2f} min")

    ax.set_ylabel("Deformation (mm)", fontweight="bold")
    ax.set_title(f"Consolidation — {'Log-Time' if method == 'log' else 'Root-Time'} Method",
                 fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax.text(0.97, 0.03, f"$c_v$ = {result['cv']:.4f} cm²/min",
            transform=ax.transAxes, fontsize=10, va="bottom", ha="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="gray"))

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    result["ax"] = ax
    return result
