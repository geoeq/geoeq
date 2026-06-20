r"""
Static and dynamic pile load test interpretation for GeoEq.

This module provides safety-critical functions for interpreting pile load
tests.  All formulas follow published references and use SI units unless
noted otherwise.

Static Pile Load Test Interpretation
-------------------------------------
davisson          Davisson's Offset Limit Method (1972).
chin              Chin's Hyperbolic Method (1970).
de_beer           De Beer's Double-Log Method (1967).
hansen_80         Hansen's 80 %% Method (1963).
fhwa_5_percent    FHWA 5 %% Diameter Criterion for drilled shafts.

Dynamic Pile Testing
--------------------
case_method       Case Method — static resistance from PDA (Goble et al. 1975).
hiley             Hiley Driving Formula.
danish_formula    Danish Driving Formula.
enr               Engineering News Record Formula.

Pile Integrity
--------------
pile_impedance    Impedance Z = EA/c.
beta_integrity    Beta ratio for pile length verification.

References
----------
Davisson, M. T. (1972). "High capacity piles." *Proc., Lecture Series on
    Innovations in Foundation Construction*, ASCE Illinois Section.
Chin, F. K. (1970). "Estimation of the ultimate load of piles not carried
    to failure." *Proc., 2nd SE Asian Conf. Soil Eng.*, pp. 81--90.
De Beer, E. E. (1967). "Proefondervindelijke bijdrage tot de studie van het
    grensdraagvermogen van zand onder funderingen op staal." *Tijdschrift
    der Openbare Werken van Belgie*, 108 (6).
Hansen, J. B. (1963). "Discussion on hyperbolic stress-strain response."
    *ASCE J. Soil Mech. Found. Div.*, 89 (SM4).
Goble, G. G., Rausche, F., & Likins, G. E. (1975). "Bearing capacity of
    piles from dynamic measurements." Final Report, Dept. of Civil Eng.,
    Case Western Reserve Univ.
Fellenius, B. H. (2020). *Basics of Foundation Design*, Electronic ed.
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from geoeq.core.validation import check_positive, check_non_negative, check_range


# ---------------------------------------------------------------------------
# Static Pile Load Test Interpretation
# ---------------------------------------------------------------------------


def davisson(
    load: Union[List[float], np.ndarray],
    settlement: Union[List[float], np.ndarray],
    diameter: float,
    length: float,
    area: float,
    elastic_modulus: float,
) -> Dict[str, float]:
    r"""
    Davisson's Offset Limit Method (1972).

    The failure load is defined as the load at which the measured
    load-settlement curve intersects an *offset elastic compression line*.

    The offset line is:

    .. math::

        \Delta_{\text{offset}} = \frac{P \, L}{A \, E}
            + \left(0.15 + \frac{D}{120}\right) \times 25.4 \; \text{[mm]}

    where the term in parentheses is in **inches** and :math:`D` is the
    pile diameter in **inches**.  The constant 0.15 in. (3.81 mm) is the
    quake offset, and :math:`D/120` accounts for pile-tip displacement.

    Equivalently in SI (all mm):

    .. math::

        \Delta_{\text{offset}} = \frac{P \, L}{A \, E}
            + 3.81 + \frac{D}{4.724}

    Parameters
    ----------
    load : array_like
        Applied loads (kN), monotonically increasing.
    settlement : array_like
        Corresponding pile-head settlements (mm).
    diameter : float
        Pile diameter or equivalent diameter (mm).
    length : float
        Pile length (mm).
    area : float
        Cross-sectional area of the pile (mm^2).
    elastic_modulus : float
        Elastic modulus of the pile material (kN/mm^2 = GPa).
        Typical: steel 200 GPa, concrete 25-35 GPa.

    Returns
    -------
    dict
        ``'Qu_davisson'`` : float
            Davisson failure load (kN), or ``nan`` if no intersection.
        ``'settlement_at_failure'`` : float
            Settlement at the failure load (mm).
        ``'elastic_compression_offset'`` : float
            The fixed offset :math:`\Delta_0` (mm) = 3.81 + D/4.724.

    Notes
    -----
    * Originally formulated for driven piles up to 600 mm (24 in.) diameter.
    * For larger piles, consider the FHWA 5 %% criterion
      (:func:`fhwa_5_percent`).
    * The elastic compression slope is :math:`L / (A \, E)`.

    References
    ----------
    Davisson (1972); Fellenius (2020), Ch. 8.

    Examples
    --------
    >>> from geoeq.site.pile_load import davisson
    >>> Q = [0, 200, 400, 600, 800, 1000, 1200, 1400]
    >>> s = [0, 1.0, 2.5, 5.0, 8.5, 14.0, 22.0, 35.0]
    >>> res = davisson(Q, s, diameter=300, length=12000, area=70686, elastic_modulus=30)
    >>> 0 < res['Qu_davisson'] < 1500
    True
    """
    Q = np.asarray(load, dtype=float)
    s = np.asarray(settlement, dtype=float)
    if len(Q) != len(s):
        raise ValueError("load and settlement must have the same length.")
    if len(Q) < 3:
        raise ValueError("Need at least 3 data points.")
    check_positive(diameter, "diameter")
    check_positive(length, "length")
    check_positive(area, "area")
    check_positive(elastic_modulus, "elastic_modulus")

    # Sort by load
    idx = np.argsort(Q)
    Q = Q[idx]
    s = s[idx]

    # Convert diameter mm -> inches for offset formula, then back to mm
    D_in = diameter / 25.4
    offset_in = 0.15 + D_in / 120.0  # inches
    offset_mm = offset_in * 25.4      # mm

    # Elastic compression slope: delta_elastic = P * L / (A * E)
    # slope_elastic has units mm/kN
    slope_elastic = length / (area * elastic_modulus)

    # Offset line: s_offset(P) = P * slope_elastic + offset_mm
    s_offset = Q * slope_elastic + offset_mm

    # Find intersection: where measured settlement crosses the offset line
    # diff = s_measured - s_offset; intersection where diff changes sign
    diff = s - s_offset

    Qu = float("nan")
    s_fail = float("nan")

    for i in range(1, len(diff)):
        if diff[i - 1] <= 0 and diff[i] > 0:
            # Linear interpolation
            frac = diff[i - 1] / (diff[i - 1] - diff[i])
            Qu = float(Q[i - 1] + frac * (Q[i] - Q[i - 1]))
            s_fail = float(s[i - 1] + frac * (s[i] - s[i - 1]))
            break

    return {
        "Qu_davisson": Qu,
        "settlement_at_failure": s_fail,
        "elastic_compression_offset": offset_mm,
    }


def chin(
    load: Union[List[float], np.ndarray],
    settlement: Union[List[float], np.ndarray],
    start_index: Optional[int] = None,
    end_index: Optional[int] = None,
) -> Dict[str, float]:
    r"""
    Chin's Hyperbolic Method (1970).

    Plot :math:`\Delta / Q` (settlement / load) versus :math:`\Delta`
    (settlement).  The data approach a straight line at larger
    settlements.  The ultimate capacity is:

    .. math::

        Q_u = \frac{1}{C_1}

    where :math:`C_1` is the **slope** of the linear portion of
    :math:`\Delta / Q` versus :math:`\Delta`.

    Parameters
    ----------
    load : array_like
        Applied loads (kN).  Must be > 0 for the data points used.
    settlement : array_like
        Corresponding settlements (mm).
    start_index : int, optional
        Index at which the linear portion begins (0-based).
        If ``None``, uses the last 50 %% of data points.
    end_index : int, optional
        Index at which the linear portion ends (exclusive).

    Returns
    -------
    dict
        ``'Qu_chin'`` : float
            Chin ultimate load (kN).
        ``'slope'`` : float
            Slope C1 of the linear regression (mm/kN per mm = 1/kN).
        ``'intercept'`` : float
            Intercept C2 of the linear regression.
        ``'r_squared'`` : float
            Coefficient of determination of the linear fit.

    Notes
    -----
    * Chin's method typically **overestimates** the actual failure load
      by 20--40 %%.  Apply a reduction factor (commonly 0.80) for
      design.
    * Only data points where Q > 0 are used.

    References
    ----------
    Chin (1970); Fellenius (2020), Section 8.9.

    Examples
    --------
    >>> from geoeq.site.pile_load import chin
    >>> Q = [100, 200, 400, 600, 800, 1000, 1200]
    >>> s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
    >>> res = chin(Q, s)
    >>> res['Qu_chin'] > 0
    True
    """
    Q = np.asarray(load, dtype=float)
    s = np.asarray(settlement, dtype=float)
    if len(Q) != len(s):
        raise ValueError("load and settlement must have the same length.")

    # Filter Q > 0
    mask = Q > 0
    Q = Q[mask]
    s = s[mask]
    if len(Q) < 3:
        raise ValueError("Need at least 3 data points with Q > 0.")

    # Select the linear portion
    n = len(Q)
    i0 = start_index if start_index is not None else n // 2
    i1 = end_index if end_index is not None else n
    if i0 < 0 or i1 > n or i0 >= i1:
        raise ValueError("Invalid start_index / end_index range.")

    delta_over_Q = s / Q  # mm / kN  →  units of 1/kN * mm

    # Linear regression: delta/Q = C1 * delta + C2
    x = s[i0:i1]
    y = delta_over_Q[i0:i1]

    coeffs = np.polyfit(x, y, 1)
    C1 = coeffs[0]  # slope  (1/kN)
    C2 = coeffs[1]  # intercept

    # R-squared
    y_pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    Qu = 1.0 / C1 if C1 > 0 else float("nan")

    return {
        "Qu_chin": float(Qu),
        "slope": float(C1),
        "intercept": float(C2),
        "r_squared": float(r_sq),
    }


def de_beer(
    load: Union[List[float], np.ndarray],
    settlement: Union[List[float], np.ndarray],
) -> Dict[str, float]:
    r"""
    De Beer's Double-Logarithmic Method (1967).

    Plot :math:`\log_{10}(Q)` versus :math:`\log_{10}(\Delta)`.  The
    data form two approximately straight-line segments.  The failure
    load is the load at the **intersection** of these two segments.

    The method fits two straight lines to the first and second halves
    of the double-log data and finds their intersection.

    Parameters
    ----------
    load : array_like
        Applied loads (kN).  Must be > 0.
    settlement : array_like
        Corresponding settlements (mm).  Must be > 0.

    Returns
    -------
    dict
        ``'Qu_de_beer'`` : float
            Failure load (kN) at the intersection.
        ``'settlement_at_failure'`` : float
            Settlement at failure (mm).
        ``'log_load_at_failure'`` : float
            log10(Qu).
        ``'log_settlement_at_failure'`` : float
            log10(settlement) at failure.

    Notes
    -----
    * All data points must have Q > 0 and settlement > 0.
    * The method uses an automated split: it searches for the break
      point that minimises total sum-of-squared residuals from two
      piecewise-linear fits in log-log space.

    References
    ----------
    De Beer (1967); Fellenius (2020), Section 8.8.

    Examples
    --------
    >>> from geoeq.site.pile_load import de_beer
    >>> Q = [50, 100, 200, 400, 600, 800, 1000, 1200]
    >>> s = [0.3, 0.7, 1.5, 3.5, 7.0, 14.0, 25.0, 42.0]
    >>> res = de_beer(Q, s)
    >>> res['Qu_de_beer'] > 0
    True
    """
    Q = np.asarray(load, dtype=float)
    s = np.asarray(settlement, dtype=float)
    if len(Q) != len(s):
        raise ValueError("load and settlement must have the same length.")
    if np.any(Q <= 0) or np.any(s <= 0):
        raise ValueError("All load and settlement values must be > 0.")
    if len(Q) < 4:
        raise ValueError("Need at least 4 data points.")

    logQ = np.log10(Q)
    logS = np.log10(s)

    # Find optimal break point by minimising combined residuals
    n = len(Q)
    best_ssr = float("inf")
    best_k = 2

    for k in range(2, n - 1):
        # Fit line to first segment [0..k] and second segment [k..n]
        c1 = np.polyfit(logS[:k + 1], logQ[:k + 1], 1)
        c2 = np.polyfit(logS[k:], logQ[k:], 1)

        r1 = logQ[:k + 1] - np.polyval(c1, logS[:k + 1])
        r2 = logQ[k:] - np.polyval(c2, logS[k:])
        ssr = np.sum(r1 ** 2) + np.sum(r2 ** 2)

        if ssr < best_ssr:
            best_ssr = ssr
            best_k = k

    # Final fits at best break point
    c1 = np.polyfit(logS[:best_k + 1], logQ[:best_k + 1], 1)
    c2 = np.polyfit(logS[best_k:], logQ[best_k:], 1)

    # Intersection: c1[0]*x + c1[1] = c2[0]*x + c2[1]
    if abs(c1[0] - c2[0]) < 1e-12:
        # Parallel lines — no intersection
        return {
            "Qu_de_beer": float("nan"),
            "settlement_at_failure": float("nan"),
            "log_load_at_failure": float("nan"),
            "log_settlement_at_failure": float("nan"),
        }

    logS_int = (c2[1] - c1[1]) / (c1[0] - c2[0])
    logQ_int = c1[0] * logS_int + c1[1]

    return {
        "Qu_de_beer": float(10 ** logQ_int),
        "settlement_at_failure": float(10 ** logS_int),
        "log_load_at_failure": float(logQ_int),
        "log_settlement_at_failure": float(logS_int),
    }


def hansen_80(
    load: Union[List[float], np.ndarray],
    settlement: Union[List[float], np.ndarray],
    start_index: Optional[int] = None,
    end_index: Optional[int] = None,
) -> Dict[str, float]:
    r"""
    Hansen's 80 %% Method (1963).

    Plot :math:`\sqrt{\Delta / Q}` versus :math:`\Delta`.  At large
    settlements the relationship becomes linear:

    .. math::

        \sqrt{\frac{\Delta}{Q}} = a \, \Delta + b

    The ultimate capacity is:

    .. math::

        Q_u = \frac{1}{2 \sqrt{a \, b}}

    and the settlement at failure is:

    .. math::

        \Delta_u = \frac{b}{a}

    The name "80 %% method" comes from the fact that at the ultimate
    load the pile mobilises 80 %% of its capacity for any
    settlement beyond :math:`\Delta_u`.

    Parameters
    ----------
    load : array_like
        Applied loads (kN).  Must be > 0.
    settlement : array_like
        Corresponding settlements (mm).  Must be > 0.
    start_index : int, optional
        Start of linear portion (0-based).  Default: last 50 %%.
    end_index : int, optional
        End of linear portion (exclusive).

    Returns
    -------
    dict
        ``'Qu_hansen'`` : float
            Hansen ultimate load (kN).
        ``'settlement_at_failure'`` : float
            Settlement at failure (mm) = b/a.
        ``'a'`` : float
            Slope of the linear fit.
        ``'b'`` : float
            Intercept of the linear fit.
        ``'r_squared'`` : float
            Coefficient of determination.

    Notes
    -----
    * Like Chin's method, Hansen's method also tends to **overestimate**
      the actual failure load.
    * If a or b is negative the result is physically meaningless and
      ``nan`` is returned.

    References
    ----------
    Hansen (1963); Fellenius (2020), Section 8.10.

    Examples
    --------
    >>> from geoeq.site.pile_load import hansen_80
    >>> Q = [100, 200, 400, 600, 800, 1000, 1200]
    >>> s = [0.5, 1.2, 3.0, 5.5, 9.0, 14.0, 21.0]
    >>> res = hansen_80(Q, s)
    >>> res['Qu_hansen'] > 0
    True
    """
    Q = np.asarray(load, dtype=float)
    s = np.asarray(settlement, dtype=float)
    if len(Q) != len(s):
        raise ValueError("load and settlement must have the same length.")

    mask = (Q > 0) & (s > 0)
    Q = Q[mask]
    s = s[mask]
    if len(Q) < 3:
        raise ValueError("Need at least 3 data points with Q > 0 and s > 0.")

    n = len(Q)
    i0 = start_index if start_index is not None else n // 2
    i1 = end_index if end_index is not None else n
    if i0 < 0 or i1 > n or i0 >= i1:
        raise ValueError("Invalid start_index / end_index range.")

    y = np.sqrt(s / Q)  # sqrt(mm / kN)

    x_fit = s[i0:i1]
    y_fit = y[i0:i1]

    coeffs = np.polyfit(x_fit, y_fit, 1)
    a = coeffs[0]  # slope
    b = coeffs[1]  # intercept

    # R-squared
    y_pred = np.polyval(coeffs, x_fit)
    ss_res = np.sum((y_fit - y_pred) ** 2)
    ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    if a > 0 and b > 0:
        Qu = 1.0 / (2.0 * np.sqrt(a * b))
        s_fail = b / a
    else:
        Qu = float("nan")
        s_fail = float("nan")

    return {
        "Qu_hansen": float(Qu),
        "settlement_at_failure": float(s_fail),
        "a": float(a),
        "b": float(b),
        "r_squared": float(r_sq),
    }


def fhwa_5_percent(
    load: Union[List[float], np.ndarray],
    settlement: Union[List[float], np.ndarray],
    diameter: float,
) -> Dict[str, float]:
    r"""
    FHWA / Modified Davisson 5 %% Diameter Criterion for drilled shafts.

    The failure load is the load corresponding to a settlement equal to
    5 %% of the pile (shaft) diameter:

    .. math::

        \Delta_{\text{fail}} = 0.05 \, D

    This criterion is recommended by FHWA for drilled shafts and
    large-diameter bored piles (> 600 mm / 24 in.).

    Parameters
    ----------
    load : array_like
        Applied loads (kN).
    settlement : array_like
        Corresponding settlements (mm).
    diameter : float
        Shaft diameter (mm).

    Returns
    -------
    dict
        ``'Qu_fhwa'`` : float
            Failure load (kN) at 5 %% diameter settlement.
        ``'settlement_criterion'`` : float
            5 %% of diameter (mm).

    References
    ----------
    FHWA (2010). *Drilled Shafts: Construction Procedures and LRFD
    Design Methods*, FHWA-NHI-10-016, Section 13.3.

    Examples
    --------
    >>> from geoeq.site.pile_load import fhwa_5_percent
    >>> Q = [0, 500, 1000, 1500, 2000, 2500, 3000]
    >>> s = [0, 2.0, 5.0, 12.0, 25.0, 45.0, 75.0]
    >>> res = fhwa_5_percent(Q, s, diameter=600)
    >>> res['settlement_criterion']
    30.0
    """
    Q = np.asarray(load, dtype=float)
    s = np.asarray(settlement, dtype=float)
    check_positive(diameter, "diameter")
    if len(Q) != len(s):
        raise ValueError("load and settlement must have the same length.")

    s_crit = 0.05 * diameter  # mm

    idx = np.argsort(s)
    Q_sorted = Q[idx]
    s_sorted = s[idx]

    if s_crit > s_sorted[-1]:
        Qu = float("nan")
    elif s_crit <= s_sorted[0]:
        Qu = float(Q_sorted[0])
    else:
        from scipy.interpolate import interp1d
        interp = interp1d(s_sorted, Q_sorted, kind="linear")
        Qu = float(interp(s_crit))

    return {
        "Qu_fhwa": Qu,
        "settlement_criterion": float(s_crit),
    }


# ---------------------------------------------------------------------------
# Dynamic Pile Testing (PDA / CAPWAP)
# ---------------------------------------------------------------------------


def case_method(
    F1: float,
    F2: float,
    v1: float,
    v2: float,
    impedance: float,
    Jc: float = 0.5,
) -> Dict[str, float]:
    r"""
    Case Method for dynamic pile capacity (Goble et al. 1975).

    From force and velocity measurements at the pile head at times
    :math:`t_1` (impact peak) and :math:`t_2 = t_1 + 2L/c`:

    .. math::

        R_{\text{total}} = \frac{1}{2}\bigl[(F_1 + F_2) +
            Z\,(v_1 - v_2)\bigr]

    .. math::

        R_{\text{static}} = R_{\text{total}}\,(1 - J_c) +
            \frac{1}{2}\bigl[(F_2 + Z\,v_2)\bigr]\,J_c

    Or equivalently (closed-form used in practice):

    .. math::

        R_{\text{static}} = \frac{(1-J_c)}{2}(F_1 + Z\,v_1)
            + \frac{(1+J_c)}{2}(F_2 - Z\,v_2) \; \times (-1)

    The simplified FHWA form used here:

    .. math::

        R_{\text{static}} = \frac{1}{2}\bigl[(1-J_c)(F_1+Z\,v_1)
            + (1+J_c)(F_2-Z\,v_2)\bigr] \times \frac{1}{1}

    Note: the sign convention uses compressive force positive, velocity
    positive downward.  See Rausche et al. (1985) for derivation.

    Parameters
    ----------
    F1 : float
        Force at pile head at time t1 — first peak (kN).
    F2 : float
        Force at pile head at time t2 = t1 + 2L/c (kN).
    v1 : float
        Velocity at pile head at time t1 (m/s).
    v2 : float
        Velocity at pile head at time t2 (m/s).
    impedance : float
        Pile impedance Z = E*A/c (kN*s/m).
        Compute with :func:`pile_impedance`.
    Jc : float, default 0.5
        Case damping factor (dimensionless).

        Typical ranges by soil type:

        ======== ============
        Soil     Jc
        ======== ============
        Sand     0.10 -- 0.20
        Silt     0.20 -- 0.40
        Clay     0.40 -- 0.70
        Peat     0.70 -- 1.00
        ======== ============

    Returns
    -------
    dict
        ``'R_total'`` : float — Total dynamic resistance (kN).
        ``'R_static'`` : float — Estimated static resistance (kN).
        ``'Jc'`` : float — Damping factor used.

    References
    ----------
    Goble, Rausche & Likins (1975);
    Rausche, F., Goble, G. G., & Likins, G. E. (1985). "Dynamic
    determination of pile capacity." *ASCE J. Geotech. Eng.*, 111(3).

    Examples
    --------
    >>> from geoeq.site.pile_load import case_method
    >>> res = case_method(F1=2500, F2=800, v1=3.0, v2=0.5,
    ...                   impedance=4000, Jc=0.4)
    >>> res['R_static'] > 0
    True
    """
    check_positive(impedance, "impedance")
    check_range(Jc, "Jc", 0.0, 1.0)

    Z = impedance

    R_total = 0.5 * ((F1 + F2) + Z * (v1 - v2))

    # Standard FHWA Case Method (GRL formulation)
    R_static = 0.5 * ((1 - Jc) * (F1 + Z * v1) + (1 + Jc) * (F2 - Z * v2))

    return {
        "R_total": float(R_total),
        "R_static": float(R_static),
        "Jc": float(Jc),
    }


def hiley(
    hammer_weight: float,
    drop_height: float,
    efficiency: float,
    set_per_blow: float,
    elastic_compression: float,
) -> Dict[str, float]:
    r"""
    Hiley Driving Formula — ultimate pile capacity from driving data.

    .. math::

        Q_u = \frac{W \, h \, \eta}{s + c/2}

    Parameters
    ----------
    hammer_weight : float
        Weight of the hammer, W (kN).
    drop_height : float
        Drop height, h (m).
    efficiency : float
        Hammer efficiency, :math:`\eta` (dimensionless, 0 to 1).

        Typical values:

        ======================== ===========
        Hammer type              Efficiency
        ======================== ===========
        Drop hammer (free fall)  0.75 -- 1.0
        Single-acting air/steam  0.75 -- 0.85
        Double-acting air/steam  0.70 -- 0.80
        Diesel hammer            0.70 -- 0.90
        Hydraulic hammer         0.80 -- 0.95
        ======================== ===========

    set_per_blow : float
        Permanent set per blow, s (m).  Typically measured as the
        average of the last 10 blows.  s > 0.
    elastic_compression : float
        Total temporary (elastic) compression, c (m).
        Sum of pile, soil, and driving cap compressions.
        Typical range: 0.002 -- 0.025 m.

    Returns
    -------
    dict
        ``'Qu_hiley'`` : float — Ultimate pile capacity (kN).
        ``'energy_input'`` : float — Hammer energy W*h*eta (kN*m).

    References
    ----------
    Hiley, A. (1925). "A rational pile driving formula and its
    application in practice explained." *Engineering (London)*, 119.

    Examples
    --------
    >>> from geoeq.site.pile_load import hiley
    >>> res = hiley(hammer_weight=50, drop_height=1.5, efficiency=0.85,
    ...            set_per_blow=0.010, elastic_compression=0.005)
    >>> res['Qu_hiley'] > 0
    True
    """
    check_positive(hammer_weight, "hammer_weight")
    check_positive(drop_height, "drop_height")
    check_range(efficiency, "efficiency", 0.0, 1.0)
    check_positive(set_per_blow, "set_per_blow")
    check_non_negative(elastic_compression, "elastic_compression")

    W = hammer_weight
    h = drop_height
    eta = efficiency
    s = set_per_blow
    c = elastic_compression

    energy = W * h * eta
    Qu = energy / (s + c / 2.0)

    return {
        "Qu_hiley": float(Qu),
        "energy_input": float(energy),
    }


def danish_formula(
    hammer_weight: float,
    drop_height: float,
    efficiency: float,
    set_per_blow: float,
    pile_length: float,
    pile_area: float,
    pile_modulus: float,
) -> Dict[str, float]:
    r"""
    Danish Driving Formula — ultimate pile capacity.

    .. math::

        Q_u = \frac{e_h \, W \, h}
        {s + \sqrt{\frac{e_h \, W \, h \, L}{2 \, A \, E}}}

    where :math:`e_h = \eta` is the hammer efficiency.

    Parameters
    ----------
    hammer_weight : float
        Weight of the hammer, W (kN).
    drop_height : float
        Drop height, h (m).
    efficiency : float
        Hammer efficiency, :math:`e_h` (dimensionless, 0 to 1).
    set_per_blow : float
        Permanent set per blow, s (m).  s > 0.
    pile_length : float
        Pile length, L (m).
    pile_area : float
        Cross-sectional area, A (m^2).
    pile_modulus : float
        Elastic modulus, E (kN/m^2 = kPa).
        Typical: steel ~200e6 kPa, concrete ~25e6 kPa.

    Returns
    -------
    dict
        ``'Qu_danish'`` : float — Ultimate pile capacity (kN).
        ``'elastic_compression_term'`` : float
            The square-root term (m).

    References
    ----------
    Sorensen, T. & Hansen, B. (1957). "Pile driving formulae — an
    investigation based on dimensional analysis and statistical
    analysis." *Proc., 4th ICSMFE*, London, Vol. 2, pp. 61--65.

    Examples
    --------
    >>> from geoeq.site.pile_load import danish_formula
    >>> res = danish_formula(hammer_weight=50, drop_height=1.5,
    ...     efficiency=0.85, set_per_blow=0.010,
    ...     pile_length=15, pile_area=0.07, pile_modulus=25e6)
    >>> res['Qu_danish'] > 0
    True
    """
    check_positive(hammer_weight, "hammer_weight")
    check_positive(drop_height, "drop_height")
    check_range(efficiency, "efficiency", 0.0, 1.0)
    check_positive(set_per_blow, "set_per_blow")
    check_positive(pile_length, "pile_length")
    check_positive(pile_area, "pile_area")
    check_positive(pile_modulus, "pile_modulus")

    W = hammer_weight
    h = drop_height
    eh = efficiency
    s = set_per_blow
    L = pile_length
    A = pile_area
    E = pile_modulus

    energy = eh * W * h
    elastic_term = np.sqrt(energy * L / (2.0 * A * E))
    Qu = energy / (s + elastic_term)

    return {
        "Qu_danish": float(Qu),
        "elastic_compression_term": float(elastic_term),
    }


def enr(
    hammer_weight: float,
    drop_height: float,
    set_per_blow: float,
    rebound: float = 0.0254,
    factor_of_safety: float = 6.0,
) -> Dict[str, float]:
    r"""
    Engineering News Record (ENR) Formula — allowable pile capacity.

    .. math::

        Q_a = \frac{W \, h}{F_s \, (s + c)}

    Parameters
    ----------
    hammer_weight : float
        Weight of the hammer, W (kN).
    drop_height : float
        Drop height, h (m).
    set_per_blow : float
        Permanent set per blow, s (m).
    rebound : float, default 0.0254
        Elastic rebound / temporary compression, c (m).
        Original ENR uses c = 25.4 mm (1 inch) for drop hammers,
        c = 2.54 mm (0.1 inch) for steam/air hammers.
    factor_of_safety : float, default 6.0
        Built-in factor of safety.  The original ENR formula uses
        FS = 6.  The **Modified ENR** uses a variable FS.

    Returns
    -------
    dict
        ``'Qa_enr'`` : float — Allowable pile capacity (kN).
        ``'Qu_enr'`` : float — Implied ultimate capacity = Qa * FS (kN).

    Notes
    -----
    * The ENR formula is considered **unreliable** by modern standards.
      It has a very high coefficient of variation (> 0.5) and should
      only be used for preliminary estimates or where required by
      local building codes.
    * ASCE and FHWA discourage its use for design.

    References
    ----------
    Wellington, A. M. (1893). "Piles and pile driving."
    *Engineering News Record*.
    Das (2021), Section 11.16.

    Examples
    --------
    >>> from geoeq.site.pile_load import enr
    >>> res = enr(hammer_weight=30, drop_height=1.0, set_per_blow=0.010)
    >>> res['Qa_enr'] > 0
    True
    """
    check_positive(hammer_weight, "hammer_weight")
    check_positive(drop_height, "drop_height")
    check_positive(set_per_blow, "set_per_blow")
    check_non_negative(rebound, "rebound")
    check_positive(factor_of_safety, "factor_of_safety")

    W = hammer_weight
    h = drop_height
    s = set_per_blow
    c = rebound
    FS = factor_of_safety

    Qa = (W * h) / (FS * (s + c))
    Qu = Qa * FS

    return {
        "Qa_enr": float(Qa),
        "Qu_enr": float(Qu),
    }


# ---------------------------------------------------------------------------
# Pile Integrity
# ---------------------------------------------------------------------------


def pile_impedance(
    elastic_modulus: float,
    area: float,
    wave_speed: float,
) -> Dict[str, float]:
    r"""
    Pile impedance — fundamental property for wave-equation analysis.

    .. math::

        Z = \frac{E \, A}{c}

    Parameters
    ----------
    elastic_modulus : float
        Elastic modulus, E (kN/m^2 = kPa).
        Typical: steel ~200e6 kPa, concrete ~30e6 kPa.
    area : float
        Cross-sectional area, A (m^2).
    wave_speed : float
        Stress-wave speed, c (m/s).
        Typical: steel ~5120 m/s, concrete ~3600--4000 m/s.

    Returns
    -------
    dict
        ``'impedance'`` : float — Impedance Z (kN*s/m).

    Notes
    -----
    * Impedance changes along the pile indicate defects:

      - **Decrease** in Z → necking, crack, or void.
      - **Increase** in Z → bulging or inclusion of harder material.

    * Related: :math:`c = \sqrt{E / \rho}`, where :math:`\rho` is
      the mass density (kg/m^3).

    References
    ----------
    Rausche et al. (1985); ASTM D4945 (PDA testing).

    Examples
    --------
    >>> from geoeq.site.pile_load import pile_impedance
    >>> res = pile_impedance(elastic_modulus=30e6, area=0.07,
    ...                      wave_speed=3800)
    >>> round(res['impedance'], 1)
    552.6
    """
    check_positive(elastic_modulus, "elastic_modulus")
    check_positive(area, "area")
    check_positive(wave_speed, "wave_speed")

    Z = elastic_modulus * area / wave_speed

    return {"impedance": float(Z)}


def beta_integrity(
    measured_length: float,
    known_length: float,
) -> Dict[str, float]:
    r"""
    Beta integrity factor — pile length verification from sonic test.

    .. math::

        \beta = \frac{L_{\text{measured}}}{L_{\text{known}}}

    Parameters
    ----------
    measured_length : float
        Length obtained from low-strain integrity test (m).
    known_length : float
        As-built or design pile length (m).

    Returns
    -------
    dict
        ``'beta'`` : float — Ratio of measured to known length.
        ``'interpretation'`` : str — Qualitative assessment.

    Notes
    -----
    Interpretation guide:

    ============= =========================================
    Beta range    Interpretation
    ============= =========================================
    0.95 -- 1.05  Intact — no significant defect detected
    0.80 -- 0.95  Possible defect — further investigation
    < 0.80        Likely severe defect or shortened pile
    > 1.05        Possible bulging or wave speed error
    ============= =========================================

    References
    ----------
    ASTM D5882 (2016). Standard Test Method for Low Strain Integrity
    Testing of Deep Foundations.

    Examples
    --------
    >>> from geoeq.site.pile_load import beta_integrity
    >>> res = beta_integrity(measured_length=14.5, known_length=15.0)
    >>> round(res['beta'], 3)
    0.967
    """
    check_positive(measured_length, "measured_length")
    check_positive(known_length, "known_length")

    beta = measured_length / known_length

    if 0.95 <= beta <= 1.05:
        interp = "Intact — no significant defect detected."
    elif 0.80 <= beta < 0.95:
        interp = "Possible defect — further investigation recommended."
    elif beta < 0.80:
        interp = "Likely severe defect or shortened pile."
    else:
        interp = "Possible bulging or wave speed calibration error."

    return {
        "beta": float(beta),
        "interpretation": interp,
    }
