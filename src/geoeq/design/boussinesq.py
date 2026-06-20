"""
Stress distribution beneath loaded areas.

Implements the classical elastic solutions (Boussinesq, Westergaard,
Newmark) plus the 2:1 empirical approximation.

References
----------
* Boussinesq, J. (1885). *Application des potentiels a l'etude de l'equilibre
  et du mouvement des solides elastiques*.
* Westergaard, H. M. (1938). "A problem of elasticity suggested by a problem
  in soil mechanics." *Mechanics of Solids in Honor of S. Timoshenko*.
* Fadum, R. E. (1948). "Influence values for estimating stresses in elastic
  foundations." *Proc. 2nd ICSMFE*, Vol. 3, 77-84.
* Newmark, N. M. (1935). "Simplified computation of vertical stresses
  in elastic foundations." Univ. of Illinois Eng. Exp. Station Circular 24.
* Das (2010), Ch. 6.
"""

from __future__ import annotations

from typing import Union

import numpy as np

from geoeq.core.validation import check_positive, check_non_negative

ArrayLike = Union[float, np.ndarray]


# -----------------------------------------------------------------
# 1. Point load (Boussinesq 1885)
# -----------------------------------------------------------------
def boussinesq_point(P: float, z: float, r: float = 0.0) -> float:
    """Vertical stress under a point load on a semi-infinite elastic mass.

        delta_sigma_z = (3 P / (2 pi z^2)) * (1 / (1 + (r/z)^2))^(5/2)

    Parameters
    ----------
    P : float
        Point load (kN).
    z : float
        Depth below ground surface (m).
    r : float
        Radial distance from the line of action (m).

    Returns
    -------
    delta_sigma : float
        Vertical stress increment (kPa).

    Reference
    ---------
    Boussinesq (1885); Das (2010) Eq. 6.4.
    """
    check_positive(z, "z")
    check_non_negative(r, "r")
    ratio = 1.0 / (1.0 + (r / z) ** 2)
    return float(3 * P / (2 * np.pi * z ** 2) * ratio ** 2.5)


# -----------------------------------------------------------------
# 2. Line load
# -----------------------------------------------------------------
def boussinesq_line(q: float, z: float, x: float = 0.0) -> float:
    """Vertical stress under an infinite line load q (kN/m).

        delta_sigma_z = (2 q / pi z) * 1 / (1 + (x/z)^2)^2

    Reference
    ---------
    Das (2010) Eq. 6.7.
    """
    check_positive(z, "z")
    return float(2 * q / (np.pi * z) * 1 / (1 + (x / z) ** 2) ** 2)


# -----------------------------------------------------------------
# 3. Strip load (uniform pressure, infinite length)
# -----------------------------------------------------------------
def boussinesq_strip(q: float, B: float, z: float, x: float = 0.0) -> float:
    """Vertical stress under a uniformly loaded infinite strip of width B.

        delta_sigma_z = (q / pi) * [ alpha + sin(alpha) cos(alpha + 2 beta) ]

    where alpha and beta are the angles subtended at the point by the
    edges of the strip (Das Eq. 6.10).

    Parameters
    ----------
    q : float
        Surface pressure (kPa).
    B : float
        Strip width (m).
    z : float
        Depth (m).
    x : float
        Horizontal distance from the *centre* of the strip (m).

    Reference
    ---------
    Das (2010) Eq. 6.10.
    """
    check_positive(z, "z")
    check_positive(B, "B")
    # Geometry: edges at x - B/2 and x + B/2, both at depth z.
    x1 = x - B / 2.0
    x2 = x + B / 2.0
    alpha = np.arctan(x2 / z) - np.arctan(x1 / z)
    # beta is the angle to the closer edge from the vertical (Das defn).
    beta = np.arctan(x1 / z)
    return float(q / np.pi * (alpha + np.sin(alpha) * np.cos(alpha + 2 * beta)))


# -----------------------------------------------------------------
# 4. Circular load (uniform pressure over a circle of radius R)
# -----------------------------------------------------------------
def boussinesq_circular(q: float, R: float, z: float) -> float:
    """Vertical stress on the centreline beneath a uniformly loaded circle.

        delta_sigma_z = q * [ 1 - 1 / (1 + (R/z)^2)^(3/2) ]

    Reference
    ---------
    Das (2010) Eq. 6.15.
    """
    check_positive(z, "z")
    check_positive(R, "R")
    return float(q * (1 - 1 / (1 + (R / z) ** 2) ** 1.5))


# -----------------------------------------------------------------
# 5. Rectangular load (Fadum 1948) -- corner formula
# -----------------------------------------------------------------
def newmark_influence(m: float, n: float) -> float:
    """Fadum / Newmark influence value I for a rectangle B x L at depth z.

        m = B / z,   n = L / z

    Returns
    -------
    I : float
        Dimensionless influence factor (Das Eq. 6.30). The stress at the
        corner of a uniformly loaded B x L rectangle is q * I.

    Reference
    ---------
    Fadum (1948); Newmark (1935); Das Eq. 6.30, Table 6.5.
    """
    check_positive(m, "m")
    check_positive(n, "n")
    # Fadum's closed-form.
    s = m * m + n * n + 1.0
    term1 = (2 * m * n * np.sqrt(s)) / (s + m * m * n * n) * (s + 1) / s
    # Newmark formula (Fadum 1948):
    # I = (1/4 pi) * { 2 m n sqrt(m^2+n^2+1) / (m^2+n^2+1 + m^2 n^2)
    #                  * (m^2+n^2+2)/(m^2+n^2+1)
    #                 + atan( 2 m n sqrt(m^2+n^2+1) / (m^2+n^2+1 - m^2 n^2) ) }
    a = 2 * m * n * np.sqrt(s)
    b = s + m * m * n * n
    c = (s + 1) / s
    d = s - m * m * n * n
    atan_arg = a / d
    # When d < 0, the angle is in the second quadrant and we add pi.
    if d > 0:
        atan_term = np.arctan(atan_arg)
    else:
        atan_term = np.arctan(atan_arg) + np.pi
    I = (1.0 / (4 * np.pi)) * (a / b * c + atan_term)
    return float(I)


def boussinesq_rect(
    q: float, B: float, L: float, z: float, position: str = "corner",
) -> float:
    """Vertical stress under a uniformly loaded rectangle B x L.

    Uses Fadum (1948) corner influence value, with superposition for
    'centre' (4x mB/2 x L/2 sub-rectangles) and 'edge_mid' positions.

    Parameters
    ----------
    q : float
        Surface pressure (kPa).
    B, L : float
        Rectangle dimensions (m). Convention: L >= B.
    z : float
        Depth below surface (m).
    position : str
        'corner'  -- under one corner of the rectangle
        'centre'  -- under the centre (sum of 4 sub-rectangles)
        'edge'    -- under the midpoint of a long edge (2 sub-rectangles)

    Reference
    ---------
    Fadum (1948); Das (2010) Eq. 6.29-6.31.
    """
    check_positive(z, "z")
    check_positive(B, "B")
    check_positive(L, "L")
    position = position.lower()

    if position == "corner":
        I = newmark_influence(B / z, L / z)
        return float(q * I)
    elif position in ("centre", "center"):
        # Four rectangles each B/2 x L/2 sharing a corner at the centre.
        I = newmark_influence((B / 2) / z, (L / 2) / z)
        return float(q * 4 * I)
    elif position == "edge":
        # Two rectangles each B/2 x L sharing a corner at the edge mid.
        I = newmark_influence((B / 2) / z, L / z)
        return float(q * 2 * I)
    raise ValueError("position must be 'corner', 'centre', or 'edge'.")


# -----------------------------------------------------------------
# 6. Westergaard (1938) -- accounts for thinly stratified soils.
# -----------------------------------------------------------------
def westergaard_point(P: float, z: float, r: float = 0.0, mu: float = 0.0) -> float:
    """Westergaard vertical stress for a point load on a layered medium.

        delta_sigma_z = (P / (pi z^2)) * eta / (eta^2 + (r/z)^2)^(3/2)

    where eta = sqrt( (1 - 2mu) / (2 - 2mu) ).

    For thinly stratified soils (alternating stiff/soft) the Westergaard
    solution gives lower stresses than Boussinesq.

    Reference
    ---------
    Westergaard (1938); Das (2010) Eq. 6.13.
    """
    check_positive(z, "z")
    check_non_negative(r, "r")
    if not 0 <= mu < 0.5:
        raise ValueError("mu must be in [0, 0.5).")
    eta_sq = (1 - 2 * mu) / (2 - 2 * mu) if mu < 0.5 else 0.0
    eta = np.sqrt(eta_sq) if eta_sq > 0 else 1e-9
    return float(P / (np.pi * z ** 2) * eta / (eta_sq + (r / z) ** 2) ** 1.5)


# -----------------------------------------------------------------
# 7. 2:1 approximation (empirical, widely used in practice)
# -----------------------------------------------------------------
def stress_2to1(q: float, B: float, L: float, z: float) -> float:
    """2:1 vertical-spread approximation for a B x L footing.

        delta_sigma_z = q * B * L / ((B + z) * (L + z))

    Crude but engineering-useful for spread footings up to z ~ 2B.

    Reference
    ---------
    Das (2010) Eq. 6.34.
    """
    check_positive(z, "z")
    check_positive(B, "B")
    check_positive(L, "L")
    return float(q * B * L / ((B + z) * (L + z)))


# -----------------------------------------------------------------
# 8. Stress bulb (visualization helper)
# -----------------------------------------------------------------
def stress_bulb(P: float, z_range=None, r_range=None, n: int = 50):
    """Compute a stress-bulb mesh under a point load.

    Returns
    -------
    dict with keys 'z', 'r', 'sigma' -- 2D grids suitable for contourf.

    Reference
    ---------
    Das (2010) Fig. 6.6 (qualitative isobar pattern).
    """
    if z_range is None:
        z_range = (0.1, 10.0)
    if r_range is None:
        r_range = (-5.0, 5.0)
    z = np.linspace(z_range[0], z_range[1], n)
    r = np.linspace(r_range[0], r_range[1], n)
    R, Z = np.meshgrid(r, z)
    ratio = 1.0 / (1.0 + (R / Z) ** 2)
    sigma = 3 * P / (2 * np.pi * Z ** 2) * ratio ** 2.5
    return {"z": Z, "r": R, "sigma": sigma}
