"""
Pile design -- axial capacity, skin friction, group efficiency, settlement.

The pile *interpretation* methods (Davisson, Chin, Hansen 80%, etc.) for
processing measured load tests live in ``ge.site.pile_load``. This module
is the design side: predicting capacity from soil properties.

References
----------
* Meyerhof, G. G. (1976). "Bearing capacity and settlement of pile
  foundations." *J. Geotech. Eng. Div.*, ASCE, 102(GT3), 197-228.
* Tomlinson, M. J. (1957). "The adhesion of piles driven in clay soils."
  *4th ICSMFE*, Vol. 2, 66-71.
* Vijayvergiya, V. N., Focht, J. A. (1972). "A new way to predict capacity
  of piles in clay." OTC Paper 1718.
* Burland, J. B. (1973). "Shaft friction of piles in clay -- a simple
  fundamental approach." *Ground Engineering*, 6(3), 30-42.
* Vesic, A. S. (1977). *Design of Pile Foundations*. NCHRP Synthesis 42.
* Converse, F. J. & Labarre, E. (1947). "Group efficiency of piles."
* Das, B. M. (2014). *Principles of Foundation Engineering*, Ch. 9-10.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from geoeq.core.validation import check_positive, check_non_negative, check_range


# -----------------------------------------------------------------
# 1. End-bearing (point) resistance
# -----------------------------------------------------------------
def pile_end_bearing(
    phi: float = None, sigma_v_eff: float = None, c: float = None,
    method: str = "meyerhof", Su: float = None,
) -> dict:
    """Unit tip resistance q_p (kPa) at the pile base.

    Two regimes:

    * **Sand (drained):** q_p = sigma'_v * N_q*
        Meyerhof:  N_q* from chart -- approximation N_q* = exp(pi tan phi) * tan^2(45+phi/2) but
                   capped per Meyerhof (1976) Fig 9.11.
        Vesic:     N_q* = (1 + 2 K_0) / 3 * (tan phi)^? -- simplified here.
    * **Clay (undrained):** q_p = 9 * Su    (Skempton)

    Parameters
    ----------
    phi : float
        Friction angle (degrees). Required for sand.
    sigma_v_eff : float
        Effective vertical stress at pile tip (kPa). Required for sand.
    Su : float
        Undrained shear strength at tip (kPa). Required for clay.
    method : str
        'meyerhof' | 'vesic' | 'skempton' (clay -- alias).
    c : float
        Optional drained cohesion (kPa). Adds c * Nc* contribution.

    Returns
    -------
    dict
        ``{'q_p': kPa, 'Nq*': ..., 'method': ...}``.

    Reference
    ---------
    Meyerhof (1976); Vesic (1977); Das (2014) Ch. 9.
    """
    method = method.lower()
    if Su is not None:
        # Clay -- Skempton's classical.
        check_positive(Su, "Su")
        return {"q_p": float(9 * Su), "Nc_star": 9.0, "method": "skempton_clay"}

    check_range(phi, "phi", 0, 50)
    check_positive(sigma_v_eff, "sigma_v_eff")
    phi_r = np.radians(phi)
    if method == "meyerhof":
        Nq_star = np.exp(np.pi * np.tan(phi_r)) * \
            np.tan(np.pi / 4 + phi_r / 2) ** 2
        # Meyerhof's empirical cap.
        Nq_star = min(Nq_star, 500.0)
    elif method == "vesic":
        # Reduced/refined per Vesic (1977).
        Nq_star = np.exp((np.pi / 2 + phi_r) * np.tan(phi_r)) * \
            np.tan(np.pi / 4 + phi_r / 2)
    else:
        raise ValueError("method must be 'meyerhof', 'vesic', or pass Su for clay.")
    q_p = sigma_v_eff * Nq_star
    if c:
        # Add cohesion contribution (Das Eq. 9.21).
        Nc_star = (Nq_star - 1) / np.tan(phi_r) if phi > 0 else 9.0
        q_p += c * Nc_star
    return {"q_p": float(q_p), "Nq_star": float(Nq_star), "method": method}


# -----------------------------------------------------------------
# 2. Skin friction (alpha / beta / lambda)
# -----------------------------------------------------------------
def pile_skin_friction(
    Su: float = None, sigma_v_eff: float = None,
    method: str = "alpha", phi: float = None,
    alpha: float = None, beta: float = None,
    K: float = None, delta: float = None, lambda_: float = None,
    layer_thicknesses: Sequence[float] = None,
) -> dict:
    """Unit shaft friction f_s along a pile.

    Methods (drained vs total stress):

    * **alpha (Tomlinson 1957)** -- total stress, for clays:
        f_s = alpha * Su
    * **beta (Burland 1973)** -- effective stress, for clays or sands:
        f_s = beta * sigma'_v   where beta = K * tan(delta).
        Default K = K0 = 1 - sin(phi); delta = phi.
    * **lambda (Vijayvergiya & Focht 1972)** -- offshore long piles:
        f_s,avg = lambda * (sigma'_v_avg + 2 * Su_avg)

    Parameters
    ----------
    Su : float
        Undrained shear strength (kPa). Required for alpha.
    sigma_v_eff : float
        Effective vertical stress at depth (kPa). Required for beta/lambda.
    method : str
        'alpha' | 'beta' | 'lambda'.
    alpha, beta, K, delta, lambda_ : float
        Method-specific parameters. Sensible defaults are computed.
    layer_thicknesses : sequence, optional
        For ``lambda``, weighting of multi-layer averages.

    Reference
    ---------
    Tomlinson (1957); Burland (1973); Vijayvergiya & Focht (1972).
    """
    method = method.lower()
    if method == "alpha":
        check_positive(Su, "Su")
        if alpha is None:
            # API RP 2A (1993) tabulation; simplified.
            if Su <= 25:
                alpha = 1.0
            elif Su <= 75:
                alpha = 1.0 - 0.5 * (Su - 25) / 50  # linear 1.0 -> 0.5
            elif Su <= 175:
                alpha = 0.5 - 0.0014 * (Su - 75)    # linear 0.5 -> ~0.35
            else:
                alpha = 0.35
        return {"f_s": float(alpha * Su), "alpha": float(alpha),
                "method": "alpha"}
    if method == "beta":
        check_positive(sigma_v_eff, "sigma_v_eff")
        if beta is None:
            if K is None:
                if phi is None:
                    raise ValueError("beta method needs phi (or beta/K).")
                K = 1 - np.sin(np.radians(phi))
            if delta is None:
                if phi is None:
                    raise ValueError("beta method needs phi (or delta).")
                delta = phi
            beta = K * np.tan(np.radians(delta))
        return {"f_s": float(beta * sigma_v_eff), "beta": float(beta),
                "method": "beta"}
    if method == "lambda":
        check_positive(sigma_v_eff, "sigma_v_eff")
        check_positive(Su, "Su")
        if lambda_ is None:
            # Vijayvergiya & Focht (1972) chart, approximated.
            # lambda decreases with embedment; rough default 0.3.
            lambda_ = 0.3
        f_s = lambda_ * (sigma_v_eff + 2 * Su)
        return {"f_s": float(f_s), "lambda": float(lambda_),
                "method": "lambda"}
    raise ValueError("method must be 'alpha', 'beta', or 'lambda'.")


# -----------------------------------------------------------------
# 3. Total axial capacity
# -----------------------------------------------------------------
def pile_capacity(
    D: float, L: float,
    q_p: float, f_s: float, area_base: float = None,
    perimeter: float = None, FS: float = 3.0,
) -> dict:
    """Axial capacity Q_ult = Q_p + Q_s = q_p * A_p + f_s * As.

    Parameters
    ----------
    D : float
        Pile diameter (m).
    L : float
        Pile length (embedded depth, m).
    q_p : float
        Unit tip resistance (kPa) -- from ``pile_end_bearing``.
    f_s : float
        Unit shaft friction (kPa) -- from ``pile_skin_friction``.
        Pass the depth-averaged value, or pre-multiply by length.
    area_base : float
        Base area (m^2). Default pi*D^2/4 (closed-ended round pile).
    perimeter : float
        Shaft perimeter (m). Default pi*D.
    FS : float
        Factor of safety on Q_ult.

    Reference
    ---------
    Das (2014) Eq. 9.20.
    """
    check_positive(D, "D")
    check_positive(L, "L")
    if area_base is None:
        area_base = np.pi * D ** 2 / 4.0
    if perimeter is None:
        perimeter = np.pi * D
    Q_p = q_p * area_base
    Q_s = f_s * perimeter * L
    Q_ult = Q_p + Q_s
    Q_all = Q_ult / FS
    return {
        "Q_p": float(Q_p), "Q_s": float(Q_s),
        "Q_ult": float(Q_ult), "Q_all": float(Q_all),
        "FS": float(FS),
    }


# -----------------------------------------------------------------
# 4. Group efficiency
# -----------------------------------------------------------------
def pile_group_efficiency(
    n: int, m: int, D: float, s: float, method: str = "converse_labarre",
) -> dict:
    """Group efficiency factor eta for an n x m pile group.

    Methods:
    * **Converse-Labarre**:
        eta = 1 - [theta * ((n-1)*m + (m-1)*n)] / (90 * n * m)
        with theta = atan(D/s) in degrees.
    * **Feld**: each pile loses 1/16 capacity for each adjacent pile.

    Reference
    ---------
    Converse & Labarre (1947); Feld (1943); Das (2014) Eq. 9.74.
    """
    if n < 1 or m < 1:
        raise ValueError("n, m must be >= 1.")
    check_positive(D, "D")
    check_positive(s, "s")
    method = method.lower()
    if method in ("converse_labarre", "cl"):
        theta = np.degrees(np.arctan(D / s))
        eta = 1 - theta * ((n - 1) * m + (m - 1) * n) / (90 * n * m)
    elif method == "feld":
        # Per pile: 1 - n_adj * 1/16. Average over the group.
        # Corner: 3 adj; edge: 5; interior: 8.
        total = 0.0
        for i in range(n):
            for j in range(m):
                adj = 0
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        if (di, dj) == (0, 0):
                            continue
                        if 0 <= i + di < n and 0 <= j + dj < m:
                            adj += 1
                total += max(0, 1 - adj / 16)
        eta = total / (n * m)
    else:
        raise ValueError("method must be 'converse_labarre' or 'feld'.")
    return {"eta": float(eta), "method": method,
            "n_piles": n * m, "spacing_ratio": s / D}


# -----------------------------------------------------------------
# 5. Pile-head settlement (Vesic semi-empirical)
# -----------------------------------------------------------------
def pile_settlement(
    Q_w: float, Q_p: float, Q_s: float, D: float, L: float, Es: float,
    Ep: float = 25e6, mu_s: float = 0.3,
    Cp: float = 0.03, Cs: float = None, qp_ult: float = None,
) -> dict:
    """Vesic's three-component settlement of a single pile.

        s = s1 + s2 + s3

    Parameters
    ----------
    Q_w : float
        Working axial load (kN).
    Q_p : float
        Tip load at working condition (kN).
    Q_s : float
        Shaft load at working condition (kN).
    D, L : float
        Pile diameter and length (m).
    Es : float
        Soil elastic modulus (kPa).
    Ep : float
        Pile material modulus (kPa). Default 25 GPa (concrete).
    mu_s : float
        Soil Poisson's ratio.
    Cp : float
        Empirical tip-settlement coefficient (Das Table 9.5). Default 0.03.
    Cs : float, optional
        Shaft-settlement coefficient. Default 0.93 + 0.16 sqrt(L/D) * Cp.
    qp_ult : float
        Ultimate tip resistance (kPa). If given, used in s2 calc.

    Reference
    ---------
    Vesic (1977); Das (2014) Eq. 9.83-9.86.
    """
    check_positive(D, "D")
    check_positive(L, "L")
    check_positive(Es, "Es")
    check_positive(Ep, "Ep")
    # s1: elastic compression of pile shaft.
    # s1 = (Q_p + xi * Q_s) * L / (A_p * E_p)   with xi ~ 0.5 to 0.67 (avg).
    xi = 0.5
    A_p = np.pi * D ** 2 / 4
    s1 = (Q_p + xi * Q_s) * L / (A_p * Ep)
    # s2: tip settlement caused by Q_p.
    if qp_ult is None:
        qp_ult = (Q_p / A_p) * 3  # rough -- mobilization at FS=3
    s2 = (Q_p / A_p) * D * (1 - mu_s ** 2) / Es * Cp  # Vesic semi-empirical
    if Cs is None:
        Cs = (0.93 + 0.16 * np.sqrt(L / D)) * Cp
    # s3: shaft skin friction transferred load.
    s3 = (Q_s / (np.pi * D * L)) * D * (1 - mu_s ** 2) / Es * Cs
    s_total = s1 + s2 + s3
    return {
        "s1_elastic": float(s1), "s2_tip": float(s2), "s3_shaft": float(s3),
        "s_total": float(s_total),
    }
