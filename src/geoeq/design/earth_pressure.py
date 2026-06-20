"""
Lateral earth pressure.

K0 (at-rest), Ka (active), Kp (passive) -- Rankine and Coulomb theories,
plus distributions including water, surcharge, and cohesion (tension
cracks).

References
----------
* Coulomb, C. A. (1776). "Essai sur une application des regles de
  maximis et minimis a quelques problemes de statique relatifs a
  l'architecture." *Memoires de Mathematique et de Physique*, Paris.
* Rankine, W. J. M. (1857). "On the stability of loose earth."
  *Phil. Trans. Royal Society*, 147, 9-27.
* Jaky, J. (1944). "The coefficient of earth pressure at rest."
  *J. Society Hungarian Architects and Engineers*, 22, 355-358.
* Mayne, P. W., Kulhawy, F. H. (1982). "K0-OCR relationships in soil."
  *J. Geotech. Eng.*, ASCE, 108(GT6), 851-872.
* Das (2014), Ch. 7-8.
"""

from __future__ import annotations

import numpy as np

from geoeq.core.validation import check_non_negative, check_range, check_positive


def K0(phi: float, OCR: float = 1.0, method: str = "jaky") -> float:
    """At-rest earth pressure coefficient.

    Parameters
    ----------
    phi : float
        Effective friction angle (degrees).
    OCR : float
        Over-consolidation ratio (1 for NC soils).
    method : str
        'jaky'  -> K0 = 1 - sin(phi) (NC; Jaky 1944).
        'mayne' -> K0 = (1 - sin phi) * OCR^sin(phi)  (Mayne & Kulhawy 1982).
        'alpan' -> Alpan (1967) for OC clays.

    Reference
    ---------
    Jaky (1944); Mayne & Kulhawy (1982).
    """
    check_range(phi, "phi", 0, 50)
    check_positive(OCR, "OCR")
    phi_rad = np.radians(phi)
    method = method.lower()
    if method == "jaky":
        return float(1 - np.sin(phi_rad))
    if method == "mayne":
        return float((1 - np.sin(phi_rad)) * OCR ** np.sin(phi_rad))
    if method == "alpan":
        return float(0.19 + 0.233 * np.log10(max(phi, 1)))  # rough
    raise ValueError("method must be 'jaky', 'mayne', or 'alpan'.")


def Ka(phi: float, delta: float = 0.0, alpha: float = 0.0,
       beta: float = 0.0, method: str = "rankine") -> float:
    """Active earth pressure coefficient.

    Parameters
    ----------
    phi : float
        Effective friction angle of backfill (degrees).
    delta : float
        Wall-soil friction angle (degrees). Default 0 (smooth wall).
    alpha : float
        Wall back-face inclination from vertical (degrees, positive = wall
        leans toward backfill). Default 0.
    beta : float
        Backfill slope angle (degrees, positive = ascending). Default 0.
    method : str
        'rankine' -- requires delta=alpha=0; uses Rankine (1857) formula.
        'coulomb' -- general (Coulomb 1776).

    Reference
    ---------
    Rankine (1857); Coulomb (1776); Das (2014) Ch. 7-8.
    """
    check_range(phi, "phi", 0, 50)
    phi_r = np.radians(phi)
    delta_r = np.radians(delta)
    alpha_r = np.radians(alpha)
    beta_r = np.radians(beta)
    method = method.lower()

    if method == "rankine":
        if abs(beta) < 1e-6:
            return float(np.tan(np.pi / 4 - phi_r / 2) ** 2)
        # Sloping backfill, vertical smooth wall (Das Eq. 7.21).
        cos_b = np.cos(beta_r)
        cos_p = np.cos(phi_r)
        num = cos_b - np.sqrt(cos_b ** 2 - cos_p ** 2)
        den = cos_b + np.sqrt(cos_b ** 2 - cos_p ** 2)
        return float(cos_b * num / den)

    if method == "coulomb":
        # Coulomb's general Ka (Das Eq. 8.5).
        num = np.cos(phi_r - alpha_r) ** 2
        sin_term = (np.sqrt(np.sin(phi_r + delta_r)
                            * np.sin(phi_r - beta_r)
                            / (np.cos(alpha_r - beta_r)
                               * np.cos(alpha_r + delta_r))))
        den = (np.cos(alpha_r) ** 2 * np.cos(alpha_r + delta_r)
               * (1 + sin_term) ** 2)
        return float(num / den)

    raise ValueError("method must be 'rankine' or 'coulomb'.")


def Kp(phi: float, delta: float = 0.0, alpha: float = 0.0,
       beta: float = 0.0, method: str = "rankine") -> float:
    """Passive earth pressure coefficient.

    See ``Ka`` for parameters. Same conventions; passive coefficient is the
    sign-reversed version of active.

    Reference
    ---------
    Rankine (1857); Coulomb (1776); Das (2014) Ch. 7-8.
    """
    check_range(phi, "phi", 0, 50)
    phi_r = np.radians(phi)
    delta_r = np.radians(delta)
    alpha_r = np.radians(alpha)
    beta_r = np.radians(beta)
    method = method.lower()
    if method == "rankine":
        if abs(beta) < 1e-6:
            return float(np.tan(np.pi / 4 + phi_r / 2) ** 2)
        cos_b = np.cos(beta_r)
        cos_p = np.cos(phi_r)
        num = cos_b + np.sqrt(cos_b ** 2 - cos_p ** 2)
        den = cos_b - np.sqrt(cos_b ** 2 - cos_p ** 2)
        return float(cos_b * num / den)
    if method == "coulomb":
        num = np.cos(phi_r + alpha_r) ** 2
        sin_term = (np.sqrt(np.sin(phi_r + delta_r)
                            * np.sin(phi_r + beta_r)
                            / (np.cos(alpha_r - beta_r)
                               * np.cos(alpha_r - delta_r))))
        den = (np.cos(alpha_r) ** 2 * np.cos(alpha_r - delta_r)
               * (1 - sin_term) ** 2)
        return float(num / den)
    raise ValueError("method must be 'rankine' or 'coulomb'.")


def earth_pressure(
    gamma: float, H: float, phi: float, c: float = 0.0,
    kind: str = "active", surcharge: float = 0.0,
    water_table: float = None, K_method: str = "rankine",
) -> dict:
    """Lateral pressure distribution on a vertical wall of height H.

    Earth pressure at depth z:
        sigma_h = K * gamma * z   (cohesionless, no surcharge)
                + K * surcharge
                - 2 c sqrt(K)     (cohesion subtraction in active case)

    Returns the resultant force P (kN/m) and its line of action.

    Parameters
    ----------
    gamma : float
        Bulk unit weight (kN/m^3).
    H : float
        Wall height (m).
    phi : float
        Friction angle (degrees).
    c : float
        Cohesion (kPa).
    kind : str
        'active' | 'passive' | 'at_rest'.
    surcharge : float
        Uniform surface surcharge q (kPa).
    water_table : float, optional
        Depth of water table below the top of the wall (m). If given,
        an additional hydrostatic pressure (gamma_w*(z - z_w)) is added
        and the effective unit weight below z_w is gamma - gamma_w.

    Returns
    -------
    dict
        ``{'K': ..., 'P_total': ..., 'P_water': ..., 'P_soil': ...,
        'h_point': ...}``.

    Reference
    ---------
    Das (2014), Ch. 7-8.
    """
    check_positive(H, "H")
    check_non_negative(c, "c")
    check_positive(gamma, "gamma")
    kind = kind.lower()

    if kind == "active":
        K = Ka(phi, method=K_method)
    elif kind == "passive":
        K = Kp(phi, method=K_method)
    elif kind in ("at_rest", "atrest", "rest"):
        K = K0(phi)
    else:
        raise ValueError("kind must be 'active', 'passive', or 'at_rest'.")

    # Cohesion subtraction (active) or addition (passive).
    sign = -1.0 if kind == "active" else (1.0 if kind == "passive" else 0.0)
    cohesion_term = sign * 2 * c * np.sqrt(K)  # subtracted from active

    # Soil pressure -- single layer for now.
    P_soil_top = K * surcharge + cohesion_term
    P_soil_bot = K * surcharge + K * gamma * H + cohesion_term
    P_soil = 0.5 * (P_soil_top + P_soil_bot) * H

    # Water (if any).
    if water_table is not None and water_table < H:
        from geoeq.core.constants import GAMMA_WATER
        h_w = H - water_table
        P_water = 0.5 * GAMMA_WATER * h_w ** 2
    else:
        P_water = 0.0
    P_total = P_soil + P_water

    # Line of action (from base): for a trapezoid, the centroid is at
    # H * (2*top + bot) / (3*(top+bot)).
    if (P_soil_top + P_soil_bot) > 0:
        h_point = H * (P_soil_top + 2 * P_soil_bot) / (
            3 * (P_soil_top + P_soil_bot))
    else:
        h_point = H / 3
    return {
        "K": float(K), "P_soil": float(P_soil),
        "P_water": float(P_water), "P_total": float(P_total),
        "h_point": float(h_point),
        "sigma_top": float(P_soil_top), "sigma_bot": float(P_soil_bot),
    }


def tension_crack_depth(c: float, gamma: float, Ka_value: float = None,
                        phi: float = 0.0) -> float:
    """Depth of tension crack in cohesive soil behind a wall.

        z_c = 2 c / (gamma * sqrt(Ka))

    Parameters
    ----------
    c : float
        Cohesion (kPa).
    gamma : float
        Unit weight (kN/m^3).
    Ka_value : float, optional
        Active earth pressure coefficient. If None, computed from ``phi``.
    phi : float
        Friction angle (degrees). Used only if Ka_value is None.

    Reference
    ---------
    Das (2014) Eq. 7.18.
    """
    check_positive(c, "c")
    check_positive(gamma, "gamma")
    if Ka_value is None:
        Ka_value = Ka(phi, method="rankine")
    return float(2 * c / (gamma * np.sqrt(Ka_value)))


def earth_pressure_plot(gamma: float, H: float, phi: float, c: float = 0.0,
                        kind: str = "active", save_as: str = None):
    """Visualize the lateral pressure distribution behind a wall.

    Returns the Matplotlib figure.
    """
    import matplotlib.pyplot as plt
    z = np.linspace(0, H, 100)
    kind = kind.lower()
    K = (Ka(phi) if kind == "active" else
         Kp(phi) if kind == "passive" else K0(phi))
    sigma = K * gamma * z - 2 * c * np.sqrt(K) * (kind == "active") + \
        2 * c * np.sqrt(K) * (kind == "passive")
    fig, ax = plt.subplots(figsize=(5, 8))
    ax.fill_betweenx(z, 0, sigma, alpha=0.3, color="#c0392b")
    ax.plot(sigma, z, color="#c0392b", linewidth=1.8)
    ax.invert_yaxis()
    ax.set_xlabel("Lateral pressure $\\sigma_h$ (kPa)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(f"{kind.title()} earth pressure (K = {K:.3f})")
    ax.grid(alpha=0.3)
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig
