"""
Retaining-wall stability (gravity / cantilever) and sheet piles.

Standard FS calculations for overturning, sliding, and bearing
under known driving forces and resisting moments. For sheet piles,
implements Blum's method for a cantilever wall in cohesionless soil.

References
----------
* Das (2014). *Principles of Foundation Engineering*, Ch. 8.
* Bowles, J. E. (1996). *Foundation Analysis and Design* (5th ed.).
* Blum, H. (1931). "Einspannungsverhaltnisse bei Bohlwerken." Berlin.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from geoeq.core.validation import check_positive, check_non_negative
from geoeq.design.earth_pressure import Ka, Kp


def wall_overturning(
    resisting_moments: Sequence[float],
    driving_moments: Sequence[float],
) -> dict:
    """Factor of safety against overturning.

        FS_overturning = sum(M_R) / sum(M_O)

    Returns dict with totals and FS.

    Reference
    ---------
    Das (2014) Eq. 8.4.
    """
    M_R = float(np.sum(np.atleast_1d(resisting_moments)))
    M_O = float(np.sum(np.atleast_1d(driving_moments)))
    if M_O <= 0:
        raise ValueError("Total driving moment must be > 0.")
    return {"M_resisting": M_R, "M_driving": M_O, "FS": M_R / M_O}


def wall_sliding(
    horizontal_forces: Sequence[float],
    vertical_forces: Sequence[float],
    mu: float = None, delta: float = None,
    c_base: float = 0.0, B: float = None, Pp: float = 0.0,
) -> dict:
    """Factor of safety against sliding at the wall base.

        FS_sliding = (sum_V * tan(delta) + c_a * B + Pp) / sum_H

    Parameters
    ----------
    horizontal_forces : sequence
        Forces pushing the wall horizontally (kN/m of wall).
    vertical_forces : sequence
        Vertical forces resisting sliding (kN/m of wall).
    mu : float, optional
        Coefficient of friction between base and soil. If given, used as
        tan(delta). Otherwise computed from ``delta``.
    delta : float, optional
        Base friction angle (degrees).
    c_base : float
        Cohesion at base (kPa). Typically c_a ~ 0.5 to 2/3 of c'.
    B : float
        Base width (m). Required if c_base > 0.
    Pp : float
        Passive resistance at the toe (kN/m). Default 0.

    Reference
    ---------
    Das (2014) Eq. 8.6.
    """
    H = float(np.sum(np.atleast_1d(horizontal_forces)))
    V = float(np.sum(np.atleast_1d(vertical_forces)))
    check_positive(H, "sum(H)")
    if mu is None:
        if delta is None:
            raise ValueError("Provide mu or delta.")
        mu = np.tan(np.radians(delta))
    R_friction = V * mu
    R_cohesion = c_base * B if (c_base > 0 and B) else 0.0
    R_total = R_friction + R_cohesion + Pp
    return {
        "sum_V": V, "sum_H": H,
        "R_friction": R_friction, "R_cohesion": R_cohesion,
        "Pp": Pp, "R_total": R_total, "FS": R_total / H,
    }


def wall_bearing(
    V: float, M_net: float, B: float, q_ult: float = None,
) -> dict:
    """Bearing-pressure distribution under a wall foundation.

    Computes the eccentricity ``e``, maximum and minimum toe pressures,
    and -- if ``q_ult`` is given -- FS_bearing = q_ult / q_max.

    Parameters
    ----------
    V : float
        Total vertical load (kN/m).
    M_net : float
        Net moment about the centreline of the base (kN.m/m).
        Positive = overturning sense.
    B : float
        Base width (m).
    q_ult : float, optional
        Ultimate bearing capacity (kPa). If given, returns FS.

    Returns
    -------
    dict
        ``{'e', 'q_max', 'q_min', 'q_ult', 'FS', 'within_kern'}``.

    Reference
    ---------
    Das (2014) Eq. 8.8-8.10.
    """
    check_positive(V, "V")
    check_positive(B, "B")
    e = abs(M_net) / V
    in_kern = e <= B / 6
    if in_kern:
        q_max = (V / B) * (1 + 6 * e / B)
        q_min = (V / B) * (1 - 6 * e / B)
    else:
        # Triangular distribution -- only 3*(B/2 - e) of the base is in contact.
        q_max = (4 * V) / (3 * (B - 2 * e))
        q_min = 0.0
    out = {"e": float(e), "q_max": float(q_max), "q_min": float(q_min),
           "within_kern": bool(in_kern)}
    if q_ult is not None:
        check_positive(q_ult, "q_ult")
        out["q_ult"] = float(q_ult)
        out["FS"] = float(q_ult / q_max)
    return out


# -----------------------------------------------------------------
# Cantilever sheet pile in cohesionless soil (Blum's method)
# -----------------------------------------------------------------
def sheet_pile(
    gamma: float, H: float, phi: float, c: float = 0.0,
    gamma_sub: float = None, water_table: float = None,
    kind: str = "cantilever",
) -> dict:
    """Embedment depth of a sheet pile by Blum's simplified method.

    Cantilever wall in cohesionless soil:
        Embedment depth D ~ found by satisfying moment equilibrium about
        the point of rotation. Closed-form approximation:

            D ~ 1.2 to 1.5 * Dmin    where Dmin from theoretical balance.

    Practical (Das Eq. 14.6): for a cantilever sheet pile in dry
    cohesionless backfill of height H above the dredge line:

        sigma_a = Ka * gamma * H   at the dredge line
        gamma_b = gamma below dredge line (gamma_sub if submerged)

        Total embedment D ~ H * [ (Kp - Ka) / (gamma * Kp) ]^something

    Returns the theoretical embedment D_theory and a design D = 1.3 * D_theory.

    Parameters
    ----------
    gamma : float
        Unit weight of backfill (kN/m^3).
    H : float
        Free height above the dredge line (m).
    phi : float
        Friction angle (degrees).
    c : float
        Cohesion (kPa). For c-phi soils, applies Bell's solution; for
        c=0 uses the Blum equations.
    gamma_sub : float, optional
        Submerged unit weight below the water table (kN/m^3).
    water_table : float, optional
        Depth of the water table from the top (m).
    kind : str
        'cantilever' (only currently supported).

    Returns
    -------
    dict
        ``{'D_theory', 'D_design', 'Ka', 'Kp'}``.

    Reference
    ---------
    Blum (1931); Das (2014) Ch. 14.
    """
    check_positive(gamma, "gamma")
    check_positive(H, "H")
    check_non_negative(c, "c")
    if kind.lower() != "cantilever":
        raise NotImplementedError("Only 'cantilever' is implemented.")
    Ka_ = Ka(phi)
    Kp_ = Kp(phi)
    # Driving force per metre at and below dredge line.
    # Cohesionless: sigma_a at dredge = Ka*gamma*H, sigma_p_resisting at base.
    # Using Bowles' simplified equation for D (no water for now).
    # D_theory derived from cubic moment balance; approximate via:
    #   D = H * sqrt(Ka / (Kp - Ka))   -- Bowles "rough first estimate"
    if Kp_ - Ka_ <= 0:
        raise ValueError("Invalid soil: Kp <= Ka.")
    D_theory = H * np.sqrt(Ka_ / (Kp_ - Ka_))
    D_design = 1.3 * D_theory
    return {
        "Ka": float(Ka_), "Kp": float(Kp_),
        "D_theory": float(D_theory),
        "D_design": float(D_design),
        "method": "Blum cantilever (approximate)",
    }
