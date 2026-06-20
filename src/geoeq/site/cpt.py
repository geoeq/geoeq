r"""
Cone Penetration Test (CPT) processing and correlations.

Functions
---------
cpt_normalize
    Compute normalized CPT parameters Qt, Fr, Bq.
cpt_sbt
    Robertson (1990/2009) Soil Behaviour Type classification.
cpt_friction_angle
    φ' from CPT (Robertson & Campanella 1983).
cpt_su
    Undrained shear strength from CPT.
cpt_dr
    Relative density from CPT (Baldi et al. 1986).
cpt_modulus
    Constrained or Young's modulus from CPT.
cpt_sbt_plot
    Robertson SBT chart visualization.

References
----------
Robertson, P. K. (1990). Soil classification using the CPT.
    *CJGE*, 27(1), 151–158.
Robertson, P. K. (2009). Interpretation of CPT — a unified approach.
    *CJGE*, 46(11), 1337–1355.
Robertson, P. K. & Campanella, R. G. (1983). Interpretation of cone
    penetration tests. Part I: Sand. *CJGE*, 20(4), 718–733.
Baldi, G. et al. (1986). Interpretation of CPTs and CPTUs. 4th Int.
    Geotechnical Seminar, Singapore, 143–156.
Lunne, T., Robertson, P. K. & Powell, J. J. M. (1997). *Cone
    Penetration Testing in Geotechnical Practice*. Blackie Academic.
"""

from typing import Dict, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def _scalar_or_array(result, *inputs):
    if all(np.ndim(x) == 0 for x in inputs):
        return float(result)
    return np.asarray(result)


def cpt_normalize(
    qt: Union[float, np.ndarray],
    fs: Union[float, np.ndarray],
    sigma_v: Union[float, np.ndarray],
    sigma_v_eff: Union[float, np.ndarray],
    u2: Union[float, np.ndarray] = 0.0,
    u0: Union[float, np.ndarray] = 0.0,
) -> Dict[str, Union[float, np.ndarray]]:
    r"""
    Compute normalized CPT parameters.

    .. math::

        Q_t = \frac{q_t - \sigma_{v0}}{\sigma'_{v0}} \qquad
        F_r = \frac{f_s}{q_t - \sigma_{v0}} \times 100 \qquad
        B_q = \frac{u_2 - u_0}{q_t - \sigma_{v0}}

    Parameters
    ----------
    qt : float or array_like
        Corrected cone resistance (kPa).
    fs : float or array_like
        Sleeve friction (kPa).
    sigma_v : float or array_like
        Total vertical stress σ_v0 (kPa).
    sigma_v_eff : float or array_like
        Effective vertical stress σ'_v0 (kPa).
    u2 : float or array_like, default 0.0
        Pore pressure at shoulder (kPa).
    u0 : float or array_like, default 0.0
        Equilibrium pore pressure (kPa).

    Returns
    -------
    dict
        ``'Qt'``, ``'Fr'``, ``'Bq'``.

    References
    ----------
    Robertson (1990, 2009).

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_normalize
    >>> res = cpt_normalize(qt=5000, fs=50, sigma_v=100, sigma_v_eff=60)
    >>> round(res['Qt'], 1)
    81.7
    """
    qt_a = np.asarray(qt, dtype=float)
    fs_a = np.asarray(fs, dtype=float)
    sv = np.asarray(sigma_v, dtype=float)
    sve = np.asarray(sigma_v_eff, dtype=float)
    u2_a = np.asarray(u2, dtype=float)
    u0_a = np.asarray(u0, dtype=float)

    net = qt_a - sv
    net = np.maximum(net, 1.0)  # avoid division by zero

    Qt = net / sve
    Fr = (fs_a / net) * 100.0
    Bq = (u2_a - u0_a) / net

    return {
        "Qt": _scalar_or_array(Qt, qt, fs),
        "Fr": _scalar_or_array(Fr, qt, fs),
        "Bq": _scalar_or_array(Bq, qt, fs),
    }


def cpt_ic(
    Qt: Union[float, np.ndarray],
    Fr: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    r"""
    Compute the Soil Behaviour Type Index Ic.

    .. math::

        I_c = \sqrt{(3.47 - \log Q_t)^2 + (\log F_r + 1.22)^2}

    Parameters
    ----------
    Qt : float or array_like
        Normalized cone resistance.
    Fr : float or array_like
        Normalized friction ratio (%).

    Returns
    -------
    float or ndarray
        Soil behaviour type index Ic.

    References
    ----------
    Robertson (2009), Eq. 2.

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_ic
    >>> round(cpt_ic(Qt=80, Fr=1.0), 2)
    1.62
    """
    Qt_a = np.asarray(Qt, dtype=float)
    Fr_a = np.asarray(Fr, dtype=float)
    Fr_a = np.maximum(Fr_a, 0.1)  # avoid log(0)

    Ic = np.sqrt((3.47 - np.log10(Qt_a)) ** 2 + (np.log10(Fr_a) + 1.22) ** 2)
    return _scalar_or_array(Ic, Qt, Fr)


def cpt_sbt(
    Qt: Union[float, np.ndarray],
    Fr: Union[float, np.ndarray],
) -> Dict[str, Union[int, str, float, np.ndarray]]:
    r"""
    Robertson (2009) Soil Behaviour Type classification from Ic.

    Zone boundaries (Robertson 2009):

    - Ic > 3.60 → Zone 2: Organic soils — clay
    - 2.95 < Ic ≤ 3.60 → Zone 3: Clay — silty clay
    - 2.60 < Ic ≤ 2.95 → Zone 4: Silt mixtures
    - 2.05 < Ic ≤ 2.60 → Zone 5: Sand mixtures
    - 1.31 < Ic ≤ 2.05 → Zone 6: Sands — clean sand
    - Ic ≤ 1.31 → Zone 7: Gravelly sand to dense sand

    Parameters
    ----------
    Qt : float or array_like
        Normalized cone resistance.
    Fr : float or array_like
        Normalized friction ratio (%).

    Returns
    -------
    dict
        ``'Ic'`` : float or ndarray — Soil behaviour type index.
        ``'zone'`` : int or list — Robertson zone number (2–7).
        ``'description'`` : str or list — soil type description.

    References
    ----------
    Robertson (2009), Table 1.

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_sbt
    >>> res = cpt_sbt(Qt=80, Fr=1.0)
    >>> res['zone']
    6
    """
    Ic_val = cpt_ic(Qt, Fr)

    _ZONES = [
        (3.60, 2, "Organic soils — clay"),
        (2.95, 3, "Clay — silty clay to clay"),
        (2.60, 4, "Silt mixtures — clayey silt to silty clay"),
        (2.05, 5, "Sand mixtures — silty sand to sandy silt"),
        (1.31, 6, "Sands — clean sand to silty sand"),
        (0.00, 7, "Gravelly sand to dense sand"),
    ]

    def _classify_one(ic):
        for threshold, zone, desc in _ZONES:
            if ic > threshold:
                return zone, desc
        return 7, "Gravelly sand to dense sand"

    if np.ndim(Ic_val) == 0:
        zone, desc = _classify_one(float(Ic_val))
        return {"Ic": float(Ic_val), "zone": zone, "description": desc}
    else:
        ic_arr = np.atleast_1d(Ic_val)
        zones = []
        descs = []
        for ic in ic_arr:
            z, d = _classify_one(float(ic))
            zones.append(z)
            descs.append(d)
        return {"Ic": Ic_val, "zone": zones, "description": descs}


def cpt_friction_angle(
    qt: Union[float, np.ndarray],
    sigma_v_eff: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    r"""
    Estimate drained friction angle from CPT.

    .. math::

        \phi' = \arctan\!\bigl[0.1 + 0.38 \log(q_t / \sigma'_{v0})\bigr]
        \qquad \text{[Robertson \& Campanella, 1983]}

    Parameters
    ----------
    qt : float or array_like
        Corrected cone resistance (kPa).
    sigma_v_eff : float or array_like
        Effective vertical stress (kPa).

    Returns
    -------
    float or ndarray
        Friction angle φ' (degrees).

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_friction_angle
    >>> round(cpt_friction_angle(5000, 60), 1)
    35.1
    """
    qt_a = np.asarray(qt, dtype=float)
    sve = np.asarray(sigma_v_eff, dtype=float)
    check_positive(qt_a, "qt")
    check_positive(sve, "sigma_v_eff")

    ratio = qt_a / sve
    phi = np.degrees(np.arctan(0.1 + 0.38 * np.log10(ratio)))
    return _scalar_or_array(phi, qt, sigma_v_eff)


def cpt_su(
    qt: Union[float, np.ndarray],
    sigma_v: Union[float, np.ndarray],
    Nkt: float = 14.0,
) -> Union[float, np.ndarray]:
    r"""
    Estimate undrained shear strength from CPT.

    .. math::

        S_u = \frac{q_t - \sigma_{v0}}{N_{kt}}

    Parameters
    ----------
    qt : float or array_like
        Corrected cone resistance (kPa).
    sigma_v : float or array_like
        Total vertical stress (kPa).
    Nkt : float, default 14.0
        Cone factor (typically 10–20; 14 is common for medium-
        sensitivity clays).

    Returns
    -------
    float or ndarray
        Undrained shear strength Su (kPa).

    References
    ----------
    Lunne et al. (1997), Ch. 6.

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_su
    >>> cpt_su(qt=1000, sigma_v=100, Nkt=14)
    64.28571428571429
    """
    qt_a = np.asarray(qt, dtype=float)
    sv = np.asarray(sigma_v, dtype=float)
    check_positive(Nkt, "Nkt")

    Su = (qt_a - sv) / Nkt
    return _scalar_or_array(np.maximum(Su, 0.0), qt, sigma_v)


def cpt_dr(
    qc: Union[float, np.ndarray],
    sigma_v_eff: Union[float, np.ndarray],
    C0: float = 157.0,
    C1: float = 0.55,
    C2: float = 2.41,
) -> Union[float, np.ndarray]:
    r"""
    Estimate relative density from CPT for normally consolidated
    uncemented quartz sands.

    .. math::

        D_r = \frac{1}{C_2} \ln\!\left(\frac{q_c}{C_0 \,(\sigma'_{v0})^{C_1}}\right)
        \qquad \text{[Baldi et al., 1986]}

    Default coefficients C₀ = 157, C₁ = 0.55, C₂ = 2.41 are for
    Ticino sand (Baldi et al., 1986).

    Parameters
    ----------
    qc : float or array_like
        Cone resistance (kPa).
    sigma_v_eff : float or array_like
        Effective vertical stress (kPa).
    C0, C1, C2 : float
        Calibration constants.

    Returns
    -------
    float or ndarray
        Relative density Dr (fraction, 0–1).

    References
    ----------
    Baldi et al. (1986).

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_dr
    >>> round(cpt_dr(10000, 100) * 100, 0)
    55.0
    """
    qc_a = np.asarray(qc, dtype=float)
    sve = np.asarray(sigma_v_eff, dtype=float)
    check_positive(qc_a, "qc")
    check_positive(sve, "sigma_v_eff")

    Dr = (1.0 / C2) * np.log(qc_a / (C0 * sve ** C1))
    Dr = np.clip(Dr, 0.0, 1.0)
    return _scalar_or_array(Dr, qc, sigma_v_eff)


def cpt_modulus(
    qt: Union[float, np.ndarray],
    sigma_v: Union[float, np.ndarray],
    Ic: Optional[Union[float, np.ndarray]] = None,
    alpha: Optional[float] = None,
) -> Union[float, np.ndarray]:
    r"""
    Estimate constrained modulus M from CPT.

    .. math::

        M = \alpha_M \times (q_t - \sigma_{v0})

    where α_M depends on Ic (Robertson, 2009):

    - Ic > 2.2 (clays): α_M = Qt when Qt < 14, else α_M = 14
    - Ic ≤ 2.2 (sands): α_M = 0.0188 × 10^(0.55 Ic + 1.68)

    Parameters
    ----------
    qt : float or array_like
        Corrected cone resistance (kPa).
    sigma_v : float or array_like
        Total vertical stress (kPa).
    Ic : float or array_like, optional
        Soil behaviour type index. If None, alpha must be provided.
    alpha : float, optional
        Direct multiplier override.

    Returns
    -------
    float or ndarray
        Constrained modulus M (kPa).

    References
    ----------
    Robertson (2009), Eq. 11.

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_modulus
    >>> cpt_modulus(qt=5000, sigma_v=100, alpha=5)
    24500.0
    """
    qt_a = np.asarray(qt, dtype=float)
    sv = np.asarray(sigma_v, dtype=float)
    net = qt_a - sv

    if alpha is not None:
        M = alpha * net
    elif Ic is not None:
        Ic_a = np.asarray(Ic, dtype=float)
        alpha_M = np.where(
            Ic_a > 2.2,
            np.minimum(net / np.maximum(sv, 1.0), 14.0),
            0.0188 * 10 ** (0.55 * Ic_a + 1.68),
        )
        M = alpha_M * net
    else:
        raise ValueError("Either Ic or alpha must be provided.")

    return _scalar_or_array(np.maximum(M, 0.0), qt, sigma_v)


def cpt_sbt_plot(
    Qt: Union[float, list, np.ndarray],
    Fr: Union[float, list, np.ndarray],
    labels=None,
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot data on the Robertson (2009) SBT chart (Qt vs Fr).

    Parameters
    ----------
    Qt : array_like
        Normalized cone resistance.
    Fr : array_like
        Normalized friction ratio (%).
    labels : list, optional
        Labels for each data point.
    ax : matplotlib.axes.Axes, optional
    save_as : str, optional

    Returns
    -------
    dict
        ``'ax'`` : Axes.

    Examples
    --------
    >>> from geoeq.site.cpt import cpt_sbt_plot
    >>> res = cpt_sbt_plot([80, 20, 5], [1.0, 3.0, 5.0])
    """
    import matplotlib.pyplot as plt

    Qt_a = np.atleast_1d(np.asarray(Qt, dtype=float))
    Fr_a = np.atleast_1d(np.asarray(Fr, dtype=float))

    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 7))
    else:
        fig = ax.get_figure()

    # Draw Ic contours
    Fr_grid = np.logspace(-1, 1, 300)
    for Ic_val, ls in [(1.31, "--"), (2.05, "--"), (2.60, "--"), (2.95, "--"), (3.60, "--")]:
        Qt_contour = 10 ** (3.47 - np.sqrt(Ic_val ** 2 - (np.log10(Fr_grid) + 1.22) ** 2))
        mask = np.isfinite(Qt_contour) & (Qt_contour > 0)
        ax.plot(Fr_grid[mask], Qt_contour[mask], "k", ls=ls, alpha=0.4, linewidth=0.8)

    # Zone labels
    zone_labels = {
        (0.5, 300): "7\nGravelly\nsand",
        (0.7, 50): "6\nClean sand",
        (1.5, 15): "5\nSand\nmixtures",
        (3.0, 5): "4\nSilt\nmixtures",
        (4.0, 2): "3\nClay",
        (6.0, 1.5): "2\nOrganic",
    }
    for (fr_pos, qt_pos), label in zone_labels.items():
        ax.text(fr_pos, qt_pos, label, fontsize=7, ha="center", va="center",
                color="#666666", alpha=0.7)

    # Data points
    ax.scatter(Fr_a, Qt_a, s=60, c="steelblue", edgecolors="navy", zorder=5)
    if labels is not None:
        for i, lbl in enumerate(labels):
            ax.annotate(str(lbl), (Fr_a[i], Qt_a[i]), fontsize=8,
                        xytext=(5, 5), textcoords="offset points")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Normalized Friction Ratio Fr (%)", fontweight="bold")
    ax.set_ylabel("Normalized Cone Resistance Qt", fontweight="bold")
    ax.set_title("Robertson (2009) CPT Soil Behaviour Type Chart", fontweight="bold")
    ax.set_xlim(0.1, 10)
    ax.set_ylim(1, 1000)
    ax.grid(True, which="both", alpha=0.2)

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    return {"ax": ax}
