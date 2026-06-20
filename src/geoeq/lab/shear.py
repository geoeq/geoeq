r"""
Shear strength laboratory tests for GeoEq.

Functions
---------
direct_shear
    Process direct shear test data → c, φ (ASTM D3080).
triaxial
    Process triaxial test data (UU/CU/CD) → c, φ (ASTM D2850/D4767).
unconfined
    Unconfined compression test → Su (ASTM D2166).
mohr_circle
    Draw Mohr circles and failure envelope.

References
----------
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 8.
ASTM D3080 (2011). Standard Test Method for Direct Shear Test.
ASTM D2850 (2015). Standard Test Method for UU Triaxial Compression.
ASTM D4767 (2011). Standard Test Method for CU Triaxial Compression.
ASTM D2166 (2016). Standard Test Method for Unconfined Compressive Strength.
"""

from typing import Dict, List, Optional, Union, Tuple
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def _scalar_or_array(result, *inputs):
    if all(np.ndim(x) == 0 for x in inputs):
        return float(result)
    return np.asarray(result)


def direct_shear(
    normal_stress: Union[List[float], np.ndarray],
    shear_stress: Union[List[float], np.ndarray],
) -> Dict[str, float]:
    r"""
    Process direct shear test results to obtain shear strength parameters.

    Fits the Mohr–Coulomb failure criterion by least-squares linear
    regression through the (σ', τ_f) data:

    .. math::

        \tau_f = c' + \sigma' \tan\phi'  \qquad \text{[Das Eq.\,8.3]}

    Parameters
    ----------
    normal_stress : array_like
        Effective normal stress on the failure plane for each specimen (kPa).
        Minimum 3 values required.
    shear_stress : array_like
        Peak (or residual) shear stress at failure for each specimen (kPa).

    Returns
    -------
    dict
        ``'c'`` : float — cohesion intercept (kPa).
        ``'phi'`` : float — friction angle (degrees).
        ``'r_squared'`` : float — coefficient of determination.

    Examples
    --------
    >>> from geoeq.lab.shear import direct_shear
    >>> res = direct_shear([50, 100, 150], [38, 62, 86])
    >>> round(res['phi'], 1)
    25.6
    >>> round(res['c'], 1)
    13.3
    """
    sigma = np.asarray(normal_stress, dtype=float)
    tau = np.asarray(shear_stress, dtype=float)

    if len(sigma) < 3:
        raise ValueError("Need at least 3 data points for direct shear test.")
    if len(sigma) != len(tau):
        raise ValueError("normal_stress and shear_stress must have the same length.")
    for s in sigma:
        check_non_negative(s, "normal_stress")
    for t in tau:
        check_non_negative(t, "shear_stress")

    coeffs = np.polyfit(sigma, tau, 1)
    tan_phi = coeffs[0]
    c = coeffs[1]

    phi_deg = float(np.degrees(np.arctan(tan_phi)))

    ss_res = np.sum((tau - np.polyval(coeffs, sigma)) ** 2)
    ss_tot = np.sum((tau - np.mean(tau)) ** 2)
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "c": float(max(c, 0.0)),
        "phi": phi_deg,
        "r_squared": float(r_sq),
    }


def triaxial(
    sigma3: Union[List[float], np.ndarray],
    delta_sigma: Union[List[float], np.ndarray],
    kind: str = "CD",
) -> Dict[str, float]:
    r"""
    Process triaxial compression test data to obtain shear strength parameters.

    From each specimen the major principal stress at failure is:

    .. math::

        \sigma_1 = \sigma_3 + \Delta\sigma_f

    A Mohr–Coulomb envelope tangent to the circles yields *c* and *φ*.

    For UU tests (*kind* = ``"UU"``), undrained shear strength is returned
    as :math:`S_u = \Delta\sigma_f / 2`.

    Parameters
    ----------
    sigma3 : array_like
        Cell (confining) pressure for each specimen (kPa). Minimum 2 for
        UU, 3 for CD/CU.
    delta_sigma : array_like
        Deviator stress at failure Δσ_f = σ₁ − σ₃ for each specimen (kPa).
    kind : {'UU', 'CU', 'CD'}, default ``'CD'``
        Test type:
        - ``'UU'`` — Unconsolidated-undrained (total stress, φ = 0).
        - ``'CU'`` — Consolidated-undrained (total stress parameters).
        - ``'CD'`` — Consolidated-drained (effective stress parameters).

    Returns
    -------
    dict
        ``'c'`` : float — cohesion (kPa).
        ``'phi'`` : float — friction angle (degrees).
        ``'sigma1'`` : ndarray — major principal stress at failure (kPa).
        ``'sigma3'`` : ndarray — confining pressures (kPa).
        For UU: also ``'Su'`` — average undrained shear strength (kPa).

    References
    ----------
    Das (2021), Ch. 8, Sections 8.5–8.9.

    Examples
    --------
    >>> from geoeq.lab.shear import triaxial
    >>> res = triaxial([100, 200, 300], [220, 380, 540], kind='CD')
    >>> round(res['phi'], 1)
    26.6
    """
    s3 = np.asarray(sigma3, dtype=float)
    ds = np.asarray(delta_sigma, dtype=float)
    kind = kind.upper()

    if kind not in ("UU", "CU", "CD"):
        raise ValueError(f"kind must be 'UU', 'CU', or 'CD', got '{kind}'.")
    if len(s3) != len(ds):
        raise ValueError("sigma3 and delta_sigma must have the same length.")

    min_pts = 2 if kind == "UU" else 3
    if len(s3) < min_pts:
        raise ValueError(f"Need at least {min_pts} specimens for {kind} test.")

    for v in s3:
        check_non_negative(v, "sigma3")
    for v in ds:
        check_positive(v, "delta_sigma")

    s1 = s3 + ds

    if kind == "UU":
        Su = float(np.mean(ds / 2.0))
        return {
            "c": Su,
            "phi": 0.0,
            "Su": Su,
            "sigma1": s1,
            "sigma3": s3,
        }

    # For CD/CU: fit Mohr–Coulomb from p-q space
    # p = (σ₁ + σ₃)/2, q = (σ₁ - σ₃)/2
    p = (s1 + s3) / 2.0
    q = (s1 - s3) / 2.0

    # q = a + p * tan(alpha), where sin(phi) = tan(alpha), c = a / cos(phi)
    coeffs = np.polyfit(p, q, 1)
    tan_alpha = coeffs[0]
    a = coeffs[1]

    sin_phi = tan_alpha
    sin_phi = np.clip(sin_phi, -1.0, 1.0)
    phi_rad = np.arcsin(sin_phi)
    phi_deg = float(np.degrees(phi_rad))
    cos_phi = np.cos(phi_rad)
    c = float(a / cos_phi) if cos_phi > 1e-10 else float(a)

    return {
        "c": max(c, 0.0),
        "phi": phi_deg,
        "sigma1": s1,
        "sigma3": s3,
    }


def unconfined(
    qu: Union[float, List[float], np.ndarray],
) -> Dict[str, Union[float, np.ndarray]]:
    r"""
    Process unconfined compression test — compute undrained shear strength.

    .. math::

        S_u = \frac{q_u}{2}  \qquad \text{[Das Eq.\,8.9]}

    Parameters
    ----------
    qu : float or array_like
        Unconfined compressive strength (kPa).

    Returns
    -------
    dict
        ``'qu'`` : float or ndarray — unconfined compressive strength (kPa).
        ``'Su'`` : float or ndarray — undrained shear strength (kPa).
        ``'consistency'`` : str — consistency description.

    References
    ----------
    Das (2021), Table 8.4; ASTM D2166.

    Examples
    --------
    >>> from geoeq.lab.shear import unconfined
    >>> res = unconfined(120)
    >>> res['Su']
    60.0
    >>> res['consistency']
    'Stiff'
    """
    qu_arr = np.asarray(qu, dtype=float)
    check_positive(qu_arr, "qu")

    Su = qu_arr / 2.0

    def _classify(val):
        if val < 12.5:
            return "Very soft"
        elif val < 25:
            return "Soft"
        elif val < 50:
            return "Medium"
        elif val < 100:
            return "Stiff"
        elif val < 200:
            return "Very stiff"
        else:
            return "Hard"

    if np.ndim(qu) == 0:
        consistency = _classify(float(Su))
    else:
        consistency = [_classify(float(v)) for v in np.atleast_1d(Su)]

    return {
        "qu": _scalar_or_array(qu_arr, qu),
        "Su": _scalar_or_array(Su, qu),
        "consistency": consistency,
    }


def mohr_circle(
    sigma1: Union[float, List[float], np.ndarray],
    sigma3: Union[float, List[float], np.ndarray],
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Draw Mohr circles and compute the failure envelope.

    For each pair (σ₁, σ₃), draws a Mohr circle centred at
    :math:`\bigl(\tfrac{\sigma_1+\sigma_3}{2},\;0\bigr)` with radius
    :math:`\tfrac{\sigma_1-\sigma_3}{2}`.

    If ≥ 3 circles are given, a best-fit Mohr–Coulomb envelope is drawn.

    Parameters
    ----------
    sigma1 : float or array_like
        Major principal stress at failure (kPa).
    sigma3 : float or array_like
        Minor principal stress (confining pressure) at failure (kPa).
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on.
    save_as : str, optional
        File path to save figure (e.g. ``'mohr.png'``).

    Returns
    -------
    dict
        ``'c'`` : float — cohesion intercept (kPa), or ``None`` if < 3 circles.
        ``'phi'`` : float — friction angle (degrees), or ``None``.
        ``'ax'`` : matplotlib.axes.Axes.

    References
    ----------
    Das (2021), Ch. 8, Fig. 8.8.

    Examples
    --------
    >>> from geoeq.lab.shear import mohr_circle
    >>> res = mohr_circle([320, 580, 840], [100, 200, 300])
    >>> round(res['phi'], 1)
    26.6
    """
    import matplotlib.pyplot as plt

    s1 = np.atleast_1d(np.asarray(sigma1, dtype=float))
    s3 = np.atleast_1d(np.asarray(sigma3, dtype=float))

    if len(s1) != len(s3):
        raise ValueError("sigma1 and sigma3 must have the same length.")

    for v in s3:
        check_non_negative(v, "sigma3")
    for v in s1:
        check_positive(v, "sigma1")

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.get_figure()

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(s1), 10)))

    for i, (sig1, sig3) in enumerate(zip(s1, s3)):
        center = (sig1 + sig3) / 2.0
        radius = (sig1 - sig3) / 2.0
        theta = np.linspace(0, np.pi, 200)
        x = center + radius * np.cos(theta)
        y = radius * np.sin(theta)
        ax.plot(x, y, color=colors[i], linewidth=1.5,
                label=f"σ₃={sig3:.0f}, σ₁={sig1:.0f} kPa")

    c_val, phi_val = None, None
    if len(s1) >= 3:
        p = (s1 + s3) / 2.0
        q = (s1 - s3) / 2.0
        coeffs = np.polyfit(p, q, 1)
        sin_phi = np.clip(coeffs[0], -1.0, 1.0)
        phi_rad = np.arcsin(sin_phi)
        phi_val = float(np.degrees(phi_rad))
        cos_phi = np.cos(phi_rad)
        c_val = float(coeffs[1] / cos_phi) if cos_phi > 1e-10 else float(coeffs[1])
        c_val = max(c_val, 0.0)

        sigma_env = np.linspace(0, float(np.max(s1)) * 1.1, 200)
        tau_env = c_val + sigma_env * np.tan(phi_rad)
        ax.plot(sigma_env, tau_env, 'k--', linewidth=2,
                label=f"Envelope: c={c_val:.1f} kPa, φ={phi_val:.1f}°")

    ax.set_xlabel("Normal Stress σ (kPa)", fontweight="bold")
    ax.set_ylabel("Shear Stress τ (kPa)", fontweight="bold")
    ax.set_title("Mohr Circles & Failure Envelope", fontweight="bold")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    ax.axhline(0, color="black", linewidth=0.5)

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    return {
        "c": c_val,
        "phi": phi_val,
        "ax": ax,
    }


def direct_shear_plot(
    normal_stress: Union[List[float], np.ndarray],
    shear_stress: Union[List[float], np.ndarray],
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot direct shear test results with failure envelope.

    Plots the (σ', τ_f) data and the best-fit Mohr–Coulomb line.

    Parameters
    ----------
    normal_stress : array_like
        Effective normal stress on the failure plane (kPa).
    shear_stress : array_like
        Shear stress at failure (kPa).
    ax : matplotlib.axes.Axes, optional
        Existing axes.
    save_as : str, optional
        File path to save.

    Returns
    -------
    dict
        ``'c'``, ``'phi'``, ``'r_squared'``, ``'ax'``.

    Examples
    --------
    >>> from geoeq.lab.shear import direct_shear_plot
    >>> res = direct_shear_plot([50, 100, 150], [38, 62, 86])
    """
    import matplotlib.pyplot as plt

    sigma = np.asarray(normal_stress, dtype=float)
    tau = np.asarray(shear_stress, dtype=float)

    result = direct_shear(sigma, tau)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    ax.scatter(sigma, tau, s=80, color="steelblue", edgecolors="navy",
               zorder=5, label="Test data")

    sigma_fit = np.linspace(0, float(np.max(sigma)) * 1.3, 200)
    tau_fit = result["c"] + sigma_fit * np.tan(np.radians(result["phi"]))
    ax.plot(sigma_fit, tau_fit, "r-", linewidth=2,
            label=f"c = {result['c']:.1f} kPa, φ = {result['phi']:.1f}°")

    ax.set_xlabel("Normal Stress σ' (kPa)", fontweight="bold")
    ax.set_ylabel("Shear Stress τ (kPa)", fontweight="bold")
    ax.set_title("Direct Shear Test — Failure Envelope", fontweight="bold")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    result["ax"] = ax
    return result
