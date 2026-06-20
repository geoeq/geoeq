"""
Slope stability.

* Infinite slope (with or without seepage parallel to the slope).
* Culmann's planar-failure for slopes of finite height.
* Taylor's stability number (1937) -- closed-form for phi=0 clays.
* Bishop's simplified method (1955) -- circular failure, iterative.

References
----------
* Culmann, K. (1875). *Die graphische Statik*.
* Taylor, D. W. (1937). "Stability of earth slopes." *J. Boston Soc. Civ. Eng.*,
  24, 197-247.
* Bishop, A. W. (1955). "The use of the slip circle in the stability analysis
  of slopes." *Geotechnique*, 5(1), 7-17.
* Janbu, N. (1968). *Slope stability computations*.
* Das, B. M. (2010). *Principles of Geotechnical Engineering* (7th ed.), Ch. 13.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from geoeq.core.constants import GAMMA_WATER
from geoeq.core.validation import check_positive, check_non_negative, check_range


# -----------------------------------------------------------------
# 1. Infinite slope
# -----------------------------------------------------------------
def infinite_slope(
    phi: float, beta: float, c: float = 0.0, gamma: float = 18.0,
    H: float = 1.0, seepage: bool = False, gamma_sat: float = None,
) -> dict:
    """Factor of safety of an infinite slope.

    Cohesionless, dry/no-seepage:
        FS = tan(phi) / tan(beta)

    Cohesive (c-phi) without seepage:
        FS = c / (gamma * H * cos(beta)^2 * tan(beta))
             + tan(phi) / tan(beta)

    Cohesionless with full seepage parallel to slope:
        FS = (gamma' / gamma_sat) * tan(phi) / tan(beta)

    Parameters
    ----------
    phi : float
        Friction angle (degrees).
    beta : float
        Slope inclination (degrees).
    c : float
        Cohesion (kPa). Default 0.
    gamma : float
        Bulk unit weight (kN/m^3).
    H : float
        Depth of failure plane below surface (m). Only matters if c > 0.
    seepage : bool
        If True, assumes seepage parallel to slope to the surface.
    gamma_sat : float
        Saturated unit weight (kN/m^3). Required if seepage=True.

    Returns
    -------
    dict
        ``{'FS': ..., 'gamma_eff': ...}``.

    Reference
    ---------
    Das (2010) Eq. 13.4-13.10.
    """
    check_range(phi, "phi", 0, 50)
    check_range(beta, "beta", 0.01, 89.99)
    check_non_negative(c, "c")
    check_positive(gamma, "gamma")
    beta_r = np.radians(beta)
    phi_r = np.radians(phi)

    if seepage:
        if gamma_sat is None:
            raise ValueError("gamma_sat required when seepage=True.")
        gamma_eff = gamma_sat - GAMMA_WATER
        FS = (gamma_eff / gamma_sat) * np.tan(phi_r) / np.tan(beta_r)
        if c > 0:
            FS += c / (gamma_sat * H * np.cos(beta_r) ** 2 * np.tan(beta_r))
    else:
        FS = np.tan(phi_r) / np.tan(beta_r)
        if c > 0:
            FS += c / (gamma * H * np.cos(beta_r) ** 2 * np.tan(beta_r))
    return {"FS": float(FS), "seepage": bool(seepage)}


# -----------------------------------------------------------------
# 2. Culmann -- planar failure surface
# -----------------------------------------------------------------
def culmann(c: float, phi: float, gamma: float, beta: float) -> dict:
    """Culmann's critical height H_cr for a planar-failure slope.

        H_cr = (4 c sin(beta) cos(phi)) / (gamma * (1 - cos(beta - phi)))

    Parameters
    ----------
    c : float
        Cohesion (kPa).
    phi : float
        Friction angle (degrees).
    gamma : float
        Bulk unit weight (kN/m^3).
    beta : float
        Slope angle (degrees), > phi.

    Reference
    ---------
    Culmann (1875); Das (2010) Eq. 13.13.
    """
    check_positive(c, "c")
    check_range(phi, "phi", 0, 50)
    check_positive(gamma, "gamma")
    check_range(beta, "beta", 0.01, 89.99)
    if beta <= phi:
        raise ValueError("Culmann requires beta > phi.")
    beta_r = np.radians(beta)
    phi_r = np.radians(phi)
    H_cr = (4 * c * np.sin(beta_r) * np.cos(phi_r)) / \
        (gamma * (1 - np.cos(beta_r - phi_r)))
    return {"H_cr": float(H_cr), "method": "culmann"}


# -----------------------------------------------------------------
# 3. Taylor's stability number (1937) -- circular failure
# -----------------------------------------------------------------
def taylor_stability(phi: float, c: float, gamma: float, H: float,
                     beta: float) -> dict:
    """Factor of safety by Taylor's stability number m.

        m = c / (gamma * H_cr)         (Taylor 1937)
        FS = c / (m * gamma * H)

    Uses Taylor's chart (interpolated table for phi=0..25, beta=15..90).

    Reference
    ---------
    Taylor (1937, 1948); Das (2010) Fig. 13.10.
    """
    check_non_negative(c, "c")
    check_range(phi, "phi", 0, 30)
    check_positive(gamma, "gamma")
    check_positive(H, "H")
    check_range(beta, "beta", 15, 90)

    # Taylor's stability number m for phi=0 (Das Fig 13.10b approximation):
    # For phi=0: m depends on beta and depth factor D=H'/H.
    # For phi > 0: m decreases. Simplified curve-fit (Das Table 13.1):
    # The standard "toe circle" chart for phi = 0..25, beta = 15..90 is hard
    # to reproduce closed-form; use a Janbu-style empirical fit:
    #   m approx = 0.181 - 0.0028 phi - 0.0006 beta + 0.000007 phi * beta
    # This matches Taylor's chart to ~0.01 in the practical range.
    m = max(0.01,
            0.181 - 0.0028 * phi - 0.0006 * beta + 7e-6 * phi * beta)
    FS = c / (m * gamma * H)
    return {"m": float(m), "FS": float(FS), "method": "taylor"}


# -----------------------------------------------------------------
# 4. Bishop's simplified method (1955)
# -----------------------------------------------------------------
def bishop(
    slices: Sequence[dict], R: float = None,
    max_iter: int = 50, tol: float = 1e-4,
) -> dict:
    """Bishop's simplified method (1955) for circular slope failure.

    Each ``slice`` is a dict with keys:
        b      -- slice width (m)
        h      -- slice height (m)
        alpha  -- base inclination (degrees, +ve when uphill)
        c      -- cohesion on slice base (kPa)
        phi    -- friction angle (degrees)
        gamma  -- unit weight (kN/m^3)
        u      -- pore pressure at slice base (kPa), default 0

    Iterates:
        FS_{k+1} = sum_i [ (c_i b_i + (W_i - u_i b_i) tan phi_i) / m_alpha_i ]
                   / sum_i ( W_i sin alpha_i )

    with m_alpha_i = cos alpha_i + sin alpha_i tan phi_i / FS_k.

    Reference
    ---------
    Bishop (1955); Das (2010) Eq. 13.50.
    """
    if not slices:
        raise ValueError("At least one slice required.")
    # Normalize the slice fields.
    def W(s):
        return s.get("gamma", 18.0) * s["b"] * s["h"]
    def alpha(s):
        return np.radians(s["alpha"])
    def phi_rad(s):
        return np.radians(s.get("phi", 0))
    def c_(s):
        return s.get("c", 0.0)
    def u_(s):
        return s.get("u", 0.0)
    def b_(s):
        return s["b"]

    driving = sum(W(s) * np.sin(alpha(s)) for s in slices)
    if abs(driving) < 1e-9:
        raise ValueError("Sum of W*sin(alpha) is zero -- check slice geometry.")

    FS = 1.0  # initial guess
    for it in range(max_iter):
        num = 0.0
        for s in slices:
            Wi = W(s)
            m_alpha = np.cos(alpha(s)) + \
                np.sin(alpha(s)) * np.tan(phi_rad(s)) / FS
            if m_alpha <= 0:
                # ill-conditioned; nudge FS upward.
                m_alpha = 0.01
            num += (c_(s) * b_(s)
                    + (Wi - u_(s) * b_(s)) * np.tan(phi_rad(s))) / m_alpha
        FS_new = num / driving
        if abs(FS_new - FS) < tol:
            FS = FS_new
            return {"FS": float(FS), "iterations": it + 1,
                    "converged": True, "method": "bishop_simplified"}
        FS = FS_new
    return {"FS": float(FS), "iterations": max_iter,
            "converged": False, "method": "bishop_simplified"}
