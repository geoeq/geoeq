r"""
Plate Load Test (PLT) interpretation.

Functions
---------
plt_bearing
    Bearing capacity from plate load test data.
plt_subgrade_modulus
    Modulus of subgrade reaction ks from PLT.
plt_settlement_correction
    Terzaghi & Peck settlement correction for foundation size.
plt_elastic_modulus
    Elastic modulus from plate load test.
plt_plot
    Pressure–settlement curve visualization.

References
----------
Terzaghi, K. & Peck, R. B. (1967). *Soil Mechanics in Engineering
    Practice*, 2nd ed. Wiley.
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.
ASTM D1194 (1994, withdrawn). Standard Test Method for Bearing Capacity
    of Soil for Static Load and Spread Footings (replaced by D1196).
ASTM D1196/D1196M (2012). Standard Test Method for Nonrepetitive Static
    Plate Load Tests of Soils and Flexible Pavement Components, for
    Highways and Airfields.
IS 1888 (1982). Method of Load Test on Soils.
"""

from typing import Dict, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def plt_bearing(
    pressure: Union[list, np.ndarray],
    settlement: Union[list, np.ndarray],
    failure_criterion: str = "log_log",
) -> Dict[str, float]:
    r"""
    Determine ultimate bearing capacity from plate load test data.

    Parameters
    ----------
    pressure : array_like
        Applied pressures (kPa).
    settlement : array_like
        Corresponding settlements (mm).
    failure_criterion : str, default ``'log_log'``
        - ``'log_log'`` : De Beer's double-log method — intersection
          of two straight-line segments in log–log space.
        - ``'tangent'`` : Tangent intersection method — failure at the
          intersection of the initial tangent and the final tangent.
        - ``'settlement'`` : Failure at settlement = B/10
          (where B = 0.30 m standard plate → 30 mm).

    Returns
    -------
    dict
        ``'qu'`` : float — Ultimate bearing capacity (kPa).
        ``'settlement_at_failure'`` : float — Settlement at failure (mm).

    Notes
    -----
    * Standard plate diameter is 300 mm (12 in.), but 450 mm and
      600 mm (18/24 in.) plates are also used.
    * For clay soils, the plate test gives the **bearing capacity of
      the plate**, which equals that of the actual foundation
      (independent of size for φ = 0).
    * For sandy soils, scale correction is required — see
      :func:`plt_settlement_correction`.

    References
    ----------
    Das (2021), Section 4.14; IS 1888 (1982).

    Examples
    --------
    >>> from geoeq.site.plt_test import plt_bearing
    >>> p = [0, 50, 100, 150, 200, 250, 300, 350, 400]
    >>> s = [0, 0.5, 1.2, 2.5, 5.0, 10.0, 20.0, 35.0, 55.0]
    >>> res = plt_bearing(p, s, failure_criterion='settlement')
    >>> res['qu'] > 0
    True
    """
    p = np.asarray(pressure, dtype=float)
    s = np.asarray(settlement, dtype=float)
    if len(p) != len(s):
        raise ValueError("pressure and settlement must have the same length.")
    if len(p) < 3:
        raise ValueError("Need at least 3 data points.")

    criterion = failure_criterion.lower().replace("-", "_").replace(" ", "_")

    if criterion == "settlement":
        # Failure at 30 mm (B/10 for standard 300 mm plate)
        s_fail = 30.0
        mask = s > 0
        if np.any(s >= s_fail):
            qu = float(np.interp(s_fail, s, p))
        else:
            qu = float("nan")
        return {"qu": qu, "settlement_at_failure": s_fail}

    elif criterion == "log_log":
        mask = (p > 0) & (s > 0)
        p_m = p[mask]
        s_m = s[mask]
        if len(p_m) < 4:
            raise ValueError("Need at least 4 positive data points for log-log method.")

        logp = np.log10(p_m)
        logs = np.log10(s_m)

        n = len(p_m)
        best_ssr = float("inf")
        best_k = 2
        for k in range(2, n - 1):
            c1 = np.polyfit(logp[:k + 1], logs[:k + 1], 1)
            c2 = np.polyfit(logp[k:], logs[k:], 1)
            r1 = logs[:k + 1] - np.polyval(c1, logp[:k + 1])
            r2 = logs[k:] - np.polyval(c2, logp[k:])
            ssr = np.sum(r1**2) + np.sum(r2**2)
            if ssr < best_ssr:
                best_ssr = ssr
                best_k = k

        c1 = np.polyfit(logp[:best_k + 1], logs[:best_k + 1], 1)
        c2 = np.polyfit(logp[best_k:], logs[best_k:], 1)

        if abs(c1[0] - c2[0]) < 1e-12:
            return {"qu": float("nan"), "settlement_at_failure": float("nan")}

        logp_int = (c2[1] - c1[1]) / (c1[0] - c2[0])
        logs_int = c1[0] * logp_int + c1[1]

        qu = float(10**logp_int)
        s_fail = float(10**logs_int)
        return {"qu": qu, "settlement_at_failure": s_fail}

    elif criterion == "tangent":
        # Initial and final tangent intersection
        n = len(p)
        # Initial tangent from first few points
        n_init = min(3, n // 2)
        c_init = np.polyfit(s[:n_init], p[:n_init], 1)
        # Final tangent from last few points
        n_fin = min(3, n // 2)
        c_fin = np.polyfit(s[-n_fin:], p[-n_fin:], 1)

        if abs(c_init[0] - c_fin[0]) < 1e-12:
            return {"qu": float("nan"), "settlement_at_failure": float("nan")}

        s_int = (c_fin[1] - c_init[1]) / (c_init[0] - c_fin[0])
        qu = float(c_init[0] * s_int + c_init[1])
        return {"qu": max(qu, 0.0), "settlement_at_failure": float(max(s_int, 0.0))}

    else:
        raise ValueError(
            f"Unknown criterion '{failure_criterion}'. "
            "Choose: log_log, tangent, settlement."
        )


def plt_subgrade_modulus(
    pressure: Union[float, np.ndarray],
    settlement: Union[float, np.ndarray],
    plate_size: float = 0.3,
) -> Union[float, Dict[str, float]]:
    r"""
    Modulus of subgrade reaction from plate load test.

    .. math::

        k_s = \frac{q}{\delta}

    where q is the applied pressure (kPa) and δ is the settlement (m).

    For a standard 300 mm plate, correction to foundation size:

    .. math::

        k_{s,\text{foundation}} = k_{s,\text{plate}} \times
        \left(\frac{B_p}{B_f}\right)
        \qquad \text{(clay, φ = 0)}

    .. math::

        k_{s,\text{foundation}} = k_{s,\text{plate}} \times
        \left(\frac{B_f + B_p}{2 \, B_f}\right)^{\!2}
        \qquad \text{(sand)}

    Parameters
    ----------
    pressure : float or array_like
        Applied pressure(s) (kPa).
    settlement : float or array_like
        Corresponding settlement(s) (mm).
    plate_size : float, default 0.3
        Plate diameter or width (m).

    Returns
    -------
    float or dict
        If scalar inputs: float — ks (kN/m³).
        If array inputs: dict with ``'ks'`` array and ``'mean_ks'``.

    Notes
    -----
    Settlement is converted from mm to m internally.

    References
    ----------
    Das (2021), Eq. 4.35; Terzaghi (1955).

    Examples
    --------
    >>> from geoeq.site.plt_test import plt_subgrade_modulus
    >>> plt_subgrade_modulus(pressure=200, settlement=5.0)
    40000.0
    """
    check_positive(plate_size, "plate_size")

    if np.ndim(pressure) == 0 and np.ndim(settlement) == 0:
        check_positive(pressure, "pressure")
        check_positive(settlement, "settlement")
        s_m = settlement / 1000.0  # mm → m
        return float(pressure / s_m)

    p = np.asarray(pressure, dtype=float)
    s = np.asarray(settlement, dtype=float)
    check_positive(p, "pressure")
    check_positive(s, "settlement")

    s_m = s / 1000.0
    ks = p / s_m
    return {"ks": ks, "mean_ks": float(np.mean(ks))}


def plt_settlement_correction(
    s_plate: float,
    B_plate: float,
    B_foundation: float,
    soil_type: str = "sand",
) -> float:
    r"""
    Terzaghi & Peck settlement correction from plate to foundation.

    For **sand** (Terzaghi & Peck, 1967):

    .. math::

        \frac{S_f}{S_p} = \left(\frac{B_f \,(B_p + 0.3)}
                                       {B_p \,(B_f + 0.3)}\right)^{\!2}

    For **clay** (φ = 0, bearing capacity independent of size):

    .. math::

        \frac{S_f}{S_p} = \frac{B_f}{B_p}

    Parameters
    ----------
    s_plate : float
        Settlement of the test plate (mm).
    B_plate : float
        Plate size (m).
    B_foundation : float
        Foundation size (m).
    soil_type : str, default ``'sand'``
        ``'sand'`` or ``'clay'``.

    Returns
    -------
    float
        Estimated foundation settlement (mm).

    References
    ----------
    Terzaghi & Peck (1967); Das (2021), Eq. 4.34.

    Examples
    --------
    >>> from geoeq.site.plt_test import plt_settlement_correction
    >>> round(plt_settlement_correction(s_plate=5.0, B_plate=0.3,
    ...       B_foundation=2.0, soil_type='sand'), 1)
    17.7
    """
    check_non_negative(s_plate, "s_plate")
    check_positive(B_plate, "B_plate")
    check_positive(B_foundation, "B_foundation")

    Bp = B_plate
    Bf = B_foundation

    soil = soil_type.lower()
    if soil == "sand":
        ratio = (Bf * (Bp + 0.3) / (Bp * (Bf + 0.3))) ** 2
    elif soil == "clay":
        ratio = Bf / Bp
    else:
        raise ValueError(f"Unknown soil_type '{soil_type}'. Choose: sand, clay.")

    return float(s_plate * ratio)


def plt_elastic_modulus(
    pressure: float,
    settlement: float,
    plate_diameter: float = 0.3,
    poisson_ratio: float = 0.3,
    Is: float = 0.79,
) -> float:
    r"""
    Elastic modulus from plate load test (Timoshenko & Goodier).

    .. math::

        E = \frac{q \, B \, (1 - \nu^2) \, I_s}{\delta}

    Parameters
    ----------
    pressure : float
        Applied pressure (kPa).
    settlement : float
        Settlement (mm).
    plate_diameter : float, default 0.3
        Plate diameter (m).
    poisson_ratio : float, default 0.3
        Poisson's ratio.
    Is : float, default 0.79
        Influence factor. 0.79 for rigid circular plate,
        π/4 for flexible circular plate.

    Returns
    -------
    float
        Elastic modulus E (kPa).

    References
    ----------
    Timoshenko & Goodier (1970); Das (2021), Eq. 4.33.

    Examples
    --------
    >>> from geoeq.site.plt_test import plt_elastic_modulus
    >>> round(plt_elastic_modulus(200, 5.0, plate_diameter=0.3), 0)
    8658.0
    """
    check_positive(pressure, "pressure")
    check_positive(settlement, "settlement")
    check_positive(plate_diameter, "plate_diameter")

    s_m = settlement / 1000.0  # mm → m
    E = pressure * plate_diameter * (1.0 - poisson_ratio**2) * Is / s_m
    return float(E)


def plt_plot(
    pressure: Union[list, np.ndarray],
    settlement: Union[list, np.ndarray],
    ax=None,
    save_as: Optional[str] = None,
) -> Dict:
    r"""
    Plot plate load test pressure–settlement curve.

    Parameters
    ----------
    pressure : array_like
        Applied pressures (kPa).
    settlement : array_like
        Corresponding settlements (mm).
    ax : matplotlib.axes.Axes, optional
    save_as : str, optional

    Returns
    -------
    dict
        ``'ax'`` : Axes.

    Examples
    --------
    >>> from geoeq.site.plt_test import plt_plot
    >>> res = plt_plot([0, 100, 200, 300], [0, 2, 5, 12])
    """
    import matplotlib.pyplot as plt

    p = np.asarray(pressure, dtype=float)
    s = np.asarray(settlement, dtype=float)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    else:
        fig = ax.get_figure()

    ax.plot(p, s, "o-", color="steelblue", linewidth=1.5, markersize=6,
            markerfacecolor="white", markeredgewidth=1.5)
    ax.set_xlabel("Applied Pressure (kPa)", fontweight="bold")
    ax.set_ylabel("Settlement (mm)", fontweight="bold")
    ax.set_title("Plate Load Test — Pressure vs Settlement", fontweight="bold")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)

    if save_as:
        plt.savefig(save_as, bbox_inches="tight", dpi=300)

    return {"ax": ax}
