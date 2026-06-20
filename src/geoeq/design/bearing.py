"""
Bearing capacity of shallow foundations.

Four classical methods:
    * Terzaghi (1943)  -- strips, no shape/depth factors.
    * Meyerhof (1963)  -- with shape, depth, inclination factors.
    * Hansen   (1970)  -- adds base inclination, ground inclination.
    * Vesic    (1973)  -- modifies Hansen's gamma factor.

General Meyerhof / Hansen / Vesic equation:

    q_u = c * Nc * sc * dc * ic
        + q * Nq * sq * dq * iq
        + 0.5 * gamma * B * N_gamma * s_gamma * d_gamma * i_gamma

Terzaghi (strip footing, no shape correction):

    q_u = c * Nc + q * Nq + 0.5 * gamma * B * N_gamma

References
----------
* Terzaghi, K. (1943). *Theoretical Soil Mechanics*. Wiley.
* Meyerhof, G. G. (1963). "Some recent research on the bearing capacity
  of foundations." *Canadian Geotechnical Journal*, 1(1), 16-26.
* Hansen, J. B. (1970). "A revised and extended formula for bearing
  capacity." *Bulletin No. 28*, Danish Geotechnical Institute.
* Vesic, A. S. (1973). "Analysis of ultimate loads of shallow
  foundations." *J. Soil Mech. Found. Eng.*, ASCE, 99(SM1), 45-73.
* Das, B. M. (2014). *Principles of Foundation Engineering*, Ch. 4.
"""

from __future__ import annotations

import numpy as np

from geoeq.core.validation import check_positive, check_non_negative, check_range


# -----------------------------------------------------------------
# Bearing-capacity factors Nc, Nq, N_gamma
# -----------------------------------------------------------------
def bearing_factors(phi: float, method: str = "meyerhof") -> dict:
    """Bearing-capacity factors Nc, Nq, N_gamma for a given friction angle.

    Parameters
    ----------
    phi : float
        Effective friction angle (degrees).
    method : str
        ``'terzaghi'`` | ``'meyerhof'`` | ``'hansen'`` | ``'vesic'``.

    Returns
    -------
    dict
        ``{'Nc': ..., 'Nq': ..., 'Ngamma': ...}``.

    Reference
    ---------
    Das (2014) Tables 4.1-4.3.
    """
    check_non_negative(phi, "phi")
    check_range(phi, "phi", 0, 50)
    method = method.lower()
    phi_rad = np.radians(phi)

    # Nq is identical across Meyerhof/Hansen/Vesic (Prandtl-Reissner).
    if phi == 0:
        Nq = 1.0
    else:
        Nq = np.exp(np.pi * np.tan(phi_rad)) * \
            np.tan(np.pi / 4 + phi_rad / 2) ** 2

    # Nc via classical relation Nc = (Nq - 1) / tan(phi); Nc = 5.14 at phi=0.
    if phi == 0:
        Nc = 5.14  # = pi + 2
    else:
        Nc = (Nq - 1) / np.tan(phi_rad)

    if method == "terzaghi":
        # Terzaghi's original Nq, Nc differ slightly (he used local shear);
        # use the Bowles-tabulated forms for general shear.
        Nq = np.exp((1.5 * np.pi - phi_rad) * np.tan(phi_rad)) / \
            (2 * np.cos(np.pi / 4 + phi_rad / 2) ** 2)
        if phi == 0:
            Nc = 5.7  # Terzaghi original
        else:
            Nc = (Nq - 1) / np.tan(phi_rad)
        # Terzaghi N_gamma (empirical fit to chart, Bowles 1996):
        Ngamma = (Nq - 1) * np.tan(1.4 * phi_rad) if phi > 0 else 0.0
    elif method == "meyerhof":
        Ngamma = (Nq - 1) * np.tan(1.4 * phi_rad)
    elif method == "hansen":
        Ngamma = 1.5 * (Nq - 1) * np.tan(phi_rad)
    elif method == "vesic":
        Ngamma = 2 * (Nq + 1) * np.tan(phi_rad)
    else:
        raise ValueError(
            "method must be 'terzaghi', 'meyerhof', 'hansen', or 'vesic'.")
    return {"Nc": float(Nc), "Nq": float(Nq), "Ngamma": float(Ngamma)}


# -----------------------------------------------------------------
# Shape, depth, inclination factors (Meyerhof / Hansen / Vesic)
# -----------------------------------------------------------------
def bearing_shape_factors(B: float, L: float, phi: float,
                          method: str = "meyerhof") -> dict:
    """Shape factors sc, sq, s_gamma for B x L footing.

    Reference
    ---------
    Das (2014) Table 4.3.
    """
    check_positive(B, "B")
    check_positive(L, "L")
    if B > L:
        B, L = L, B  # ensure B is the shorter side
    phi_rad = np.radians(phi)
    method = method.lower()
    Kp = np.tan(np.pi / 4 + phi_rad / 2) ** 2
    if method == "meyerhof":
        sc = 1 + 0.2 * Kp * B / L
        if phi >= 10:
            sq = s_gamma = 1 + 0.1 * Kp * B / L
        else:
            sq = s_gamma = 1.0
    elif method in ("hansen", "vesic"):
        Nq = np.exp(np.pi * np.tan(phi_rad)) * \
            np.tan(np.pi / 4 + phi_rad / 2) ** 2 if phi > 0 else 1.0
        Nc = (Nq - 1) / np.tan(phi_rad) if phi > 0 else 5.14
        sc = 1 + (Nq / Nc) * (B / L)
        sq = 1 + (B / L) * np.tan(phi_rad)
        s_gamma = max(1 - 0.4 * (B / L), 0.6)
    else:
        raise ValueError("method must be 'meyerhof', 'hansen', or 'vesic'.")
    return {"sc": float(sc), "sq": float(sq), "s_gamma": float(s_gamma)}


def bearing_depth_factors(Df: float, B: float, phi: float,
                          method: str = "meyerhof") -> dict:
    """Depth factors dc, dq, d_gamma.

    Reference
    ---------
    Das (2014) Table 4.3.
    """
    check_positive(B, "B")
    check_non_negative(Df, "Df")
    phi_rad = np.radians(phi)
    method = method.lower()
    Kp = np.tan(np.pi / 4 + phi_rad / 2) ** 2
    if method == "meyerhof":
        dc = 1 + 0.2 * np.sqrt(Kp) * Df / B
        if phi >= 10:
            dq = d_gamma = 1 + 0.1 * np.sqrt(Kp) * Df / B
        else:
            dq = d_gamma = 1.0
    elif method in ("hansen", "vesic"):
        # k factor (Hansen): Df/B if <= 1; atan(Df/B) [rad] if > 1.
        k = Df / B if Df / B <= 1 else np.arctan(Df / B)
        if phi == 0:
            dc = 1 + 0.4 * k
        else:
            dc = 1 + 0.4 * k
        dq = 1 + 2 * np.tan(phi_rad) * (1 - np.sin(phi_rad)) ** 2 * k
        d_gamma = 1.0
    else:
        raise ValueError("method must be 'meyerhof', 'hansen', or 'vesic'.")
    return {"dc": float(dc), "dq": float(dq), "d_gamma": float(d_gamma)}


def bearing_inclination_factors(beta: float, phi: float, c: float = 0.0,
                                method: str = "meyerhof") -> dict:
    """Inclination factors ic, iq, i_gamma for load inclined beta degrees
    from vertical.

    Reference
    ---------
    Das (2014) Table 4.3; Hansen (1970).
    """
    check_non_negative(beta, "beta")
    check_range(beta, "beta", 0, 90)
    method = method.lower()
    if method == "meyerhof":
        ic = iq = (1 - beta / 90) ** 2
        i_gamma = (1 - beta / phi) ** 2 if phi > 0 else 1.0
    else:
        ic = iq = (1 - beta / 90) ** 2
        i_gamma = (1 - beta / phi) ** 2 if phi > 0 else 1.0
    return {"ic": float(ic), "iq": float(iq), "i_gamma": float(i_gamma)}


# -----------------------------------------------------------------
# Ultimate bearing capacity
# -----------------------------------------------------------------
def bearing_capacity(
    c: float, gamma: float, Df: float, B: float, phi: float,
    L: float = None, method: str = "meyerhof",
    shape: bool = True, depth: bool = True,
    inclination: float = 0.0, gamma_above: float = None,
    water_table: str = "deep",
) -> dict:
    """Ultimate bearing capacity q_u of a shallow footing.

    Parameters
    ----------
    c : float
        Effective cohesion (kPa). Use Su for undrained (phi = 0) analysis.
    gamma : float
        Effective unit weight of soil *below* the footing (kN/m^3).
        Use gamma' = gamma_sat - gamma_w if the water table is at the
        base of the footing.
    Df : float
        Depth of embedment (m).
    B : float
        Footing width (m).
    phi : float
        Effective friction angle (degrees).
    L : float, optional
        Footing length (m). ``None`` -> strip footing (B/L -> 0).
    method : str
        'terzaghi' | 'meyerhof' | 'hansen' | 'vesic'.
    shape, depth : bool
        Whether to apply shape/depth corrections. Terzaghi ignores both.
    inclination : float
        Load inclination from vertical (deg). Default 0.
    gamma_above : float
        Effective unit weight above the footing (kN/m^3). Default = gamma.
    water_table : str
        'deep' (no correction) | 'at_base' (gamma below = gamma') |
        'above_base' (full submergence -- conservative).

    Returns
    -------
    dict with 'q_u' (kPa), 'Nc', 'Nq', 'Ngamma', and the applied factors.

    Reference
    ---------
    Das (2014), Ch. 4.
    """
    check_positive(B, "B")
    check_non_negative(Df, "Df")
    check_non_negative(c, "c")
    check_positive(gamma, "gamma")
    if gamma_above is None:
        gamma_above = gamma

    factors = bearing_factors(phi, method=method)
    Nc, Nq, Ng = factors["Nc"], factors["Nq"], factors["Ngamma"]

    if method == "terzaghi":
        # Apply Terzaghi's shape multipliers for square/circle (Das Eq. 4.6-4.7).
        if L is not None and np.isclose(L, B):
            # Square
            sc_, sg_ = 1.3, 0.8
        elif L is None or L > 10 * B:
            sc_, sg_ = 1.0, 1.0
        else:
            sc_, sg_ = 1.3, 0.8  # treat near-square as square; conservative
        q = gamma_above * Df
        q_u = c * Nc * sc_ + q * Nq + 0.5 * gamma * B * Ng * sg_
        return {
            "q_u": float(q_u), "Nc": Nc, "Nq": Nq, "Ngamma": Ng,
            "method": method,
        }

    # Strip default
    L_eff = L if L is not None else 1e6 * B

    sc = sq = sg = 1.0
    dc = dq = dg = 1.0
    ic = iq = ig = 1.0
    if shape:
        sf = bearing_shape_factors(B, L_eff, phi, method=method)
        sc, sq, sg = sf["sc"], sf["sq"], sf["s_gamma"]
    if depth:
        df_ = bearing_depth_factors(Df, B, phi, method=method)
        dc, dq, dg = df_["dc"], df_["dq"], df_["d_gamma"]
    if inclination > 0:
        i_ = bearing_inclination_factors(inclination, phi, c, method=method)
        ic, iq, ig = i_["ic"], i_["iq"], i_["i_gamma"]

    q = gamma_above * Df
    q_u = (c * Nc * sc * dc * ic
           + q * Nq * sq * dq * iq
           + 0.5 * gamma * B * Ng * sg * dg * ig)
    return {
        "q_u": float(q_u),
        "Nc": Nc, "Nq": Nq, "Ngamma": Ng,
        "sc": sc, "sq": sq, "s_gamma": sg,
        "dc": dc, "dq": dq, "d_gamma": dg,
        "method": method,
    }


def bearing_allowable(q_u: float, FS: float = 3.0) -> float:
    """Allowable bearing pressure q_all = q_u / FS.

    Standard FS=3 for shallow foundations under static loads
    (Das 2014, p. 224).
    """
    check_positive(FS, "FS")
    return q_u / FS
