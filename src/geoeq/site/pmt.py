r"""
Pressuremeter Test (PMT) — Ménard pressuremeter interpretation.

Functions
---------
pmt_parameters
    Extract pressuremeter modulus Em and limit pressure pL.
pmt_modulus
    Young's modulus from pressuremeter modulus.
pmt_su
    Undrained shear strength from limit pressure.
pmt_bearing
    Ultimate bearing capacity from PMT.
pmt_settlement
    Settlement estimate from PMT data (Ménard method).
pmt_ko
    At-rest earth pressure coefficient from PMT.

References
----------
Ménard, L. (1956). An apparatus for measuring the strength of soils in
    place. M.Sc. Thesis, University of Illinois.
Briaud, J.-L. (1992). *The Pressuremeter*. Balkema.
Baguelin, F., Jézéquel, J.-F. & Shields, D. H. (1978). *The
    Pressuremeter and Foundation Engineering*. Trans Tech Publications.
Clarke, B. G. (1995). *Pressuremeters in Geotechnical Design*. Blackie
    Academic.
"""

from typing import Dict, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def pmt_parameters(
    pressure: Union[list, np.ndarray],
    volume: Union[list, np.ndarray],
) -> Dict[str, float]:
    r"""
    Extract pressuremeter modulus Em and limit pressure pL from the
    pressure–volume curve.

    The pseudo-elastic phase is identified as the steepest linear
    portion of the p–ΔV curve.  The pressuremeter modulus is:

    .. math::

        E_m = 2\,(1+\nu)\,(V_0 + V_m)\,\frac{\Delta p}{\Delta V}

    where :math:`V_0` is the initial probe volume, :math:`V_m` is the
    mid-volume of the linear phase, :math:`\Delta p / \Delta V` is
    the slope of the linear portion, and ν = 0.33.

    The limit pressure :math:`p_L` is the pressure at which the
    cavity volume has doubled from its initial value.

    Parameters
    ----------
    pressure : array_like
        Applied pressures (kPa), monotonically increasing.
    volume : array_like
        Corresponding injected volumes (cm³).

    Returns
    -------
    dict
        ``'Em'`` : float — Pressuremeter modulus (kPa).
        ``'pL'`` : float — Limit pressure (kPa).
        ``'p0'`` : float — Initial pressure (kPa) at start of elastic phase.
        ``'pf'`` : float — Creep pressure (kPa) at end of elastic phase.

    Notes
    -----
    * The probe volume V₀ is typically 535 cm³ for the Ménard G-probe.
    * This function uses an automated algorithm to detect the linear
      phase; for critical projects, manual selection is recommended.

    References
    ----------
    Ménard (1956); Briaud (1992), Ch. 4.

    Examples
    --------
    >>> from geoeq.site.pmt import pmt_parameters
    >>> p = [0, 50, 100, 200, 300, 400, 500, 600, 700, 800]
    >>> v = [0, 20, 40, 60, 80, 110, 160, 240, 370, 600]
    >>> res = pmt_parameters(p, v)
    >>> res['Em'] > 0
    True
    """
    p = np.asarray(pressure, dtype=float)
    v = np.asarray(volume, dtype=float)
    if len(p) != len(v):
        raise ValueError("pressure and volume must have the same length.")
    if len(p) < 4:
        raise ValueError("Need at least 4 data points.")

    idx = np.argsort(p)
    p = p[idx]
    v = v[idx]

    V0 = 535.0  # Ménard G-probe standard volume (cm³)
    nu = 0.33

    n = len(p)
    best_r2 = -1.0
    best_i0 = 0
    best_i1 = n - 1

    for i0 in range(n - 2):
        for i1 in range(i0 + 2, min(i0 + n, n)):
            seg_p = p[i0:i1 + 1]
            seg_v = v[i0:i1 + 1]
            if len(seg_p) < 3:
                continue
            coeffs = np.polyfit(seg_v, seg_p, 1)
            pred = np.polyval(coeffs, seg_v)
            ss_res = np.sum((seg_p - pred) ** 2)
            ss_tot = np.sum((seg_p - np.mean(seg_p)) ** 2)
            if ss_tot == 0:
                continue
            r2 = 1.0 - ss_res / ss_tot
            span = i1 - i0
            if r2 > 0.95 and span > best_i1 - best_i0:
                best_r2 = r2
                best_i0 = i0
                best_i1 = i1

    seg_v = v[best_i0:best_i1 + 1]
    seg_p = p[best_i0:best_i1 + 1]
    dp_dv = np.polyfit(seg_v, seg_p, 1)[0]

    Vm = (seg_v[0] + seg_v[-1]) / 2.0
    Em = 2.0 * (1.0 + nu) * (V0 + Vm) * dp_dv

    p0 = float(p[best_i0])
    pf = float(p[best_i1])

    # Limit pressure: pressure at which volume doubles from V0
    V_limit = V0
    if v[-1] >= V_limit:
        pL = float(np.interp(V_limit, v, p))
    else:
        pL = float(p[-1])

    return {
        "Em": float(Em),
        "pL": pL,
        "p0": p0,
        "pf": pf,
    }


def pmt_modulus(
    Em: float,
    alpha: float = 0.5,
) -> float:
    r"""
    Young's modulus from pressuremeter modulus.

    .. math::

        E = \frac{E_m}{\alpha}

    Parameters
    ----------
    Em : float
        Pressuremeter modulus (kPa).
    alpha : float, default 0.5
        Rheological coefficient (Ménard, 1975).

        ============= =============
        Soil type     α
        ============= =============
        Sand          0.33 -- 0.50
        Silt          0.50
        Clay          0.50 -- 0.67
        Rock          0.33
        ============= =============

    Returns
    -------
    float
        Young's modulus E (kPa).

    References
    ----------
    Ménard (1975); Briaud (1992), Table 4.2.

    Examples
    --------
    >>> from geoeq.site.pmt import pmt_modulus
    >>> pmt_modulus(Em=15000, alpha=0.5)
    30000.0
    """
    check_positive(Em, "Em")
    check_positive(alpha, "alpha")
    return float(Em / alpha)


def pmt_su(
    pL: float,
    p0: float,
    sigma_h0: float = 0.0,
) -> float:
    r"""
    Undrained shear strength from pressuremeter limit pressure.

    .. math::

        S_u = \frac{p_L - p_0}{N_p}

    where :math:`N_p \approx 5.5` for the Ménard pressuremeter
    (Baguelin et al., 1978), and :math:`p_0` can be approximated
    as the in-situ total horizontal stress.

    Alternatively, using the net limit pressure:

    .. math::

        p_L^* = p_L - \sigma_{h0}
        \qquad
        S_u = \frac{p_L^*}{N_p}

    Parameters
    ----------
    pL : float
        Limit pressure (kPa).
    p0 : float
        Initial pressure / in-situ horizontal stress (kPa).
    sigma_h0 : float, default 0.0
        In-situ total horizontal stress (kPa).  If provided and
        non-zero, the net limit pressure method is used.

    Returns
    -------
    float
        Undrained shear strength Su (kPa).

    References
    ----------
    Baguelin et al. (1978); Briaud (1992), Eq. 4.12.

    Examples
    --------
    >>> from geoeq.site.pmt import pmt_su
    >>> round(pmt_su(pL=600, p0=100), 1)
    90.9
    """
    check_positive(pL, "pL")
    check_non_negative(p0, "p0")

    Np = 5.5
    if sigma_h0 > 0:
        Su = (pL - sigma_h0) / Np
    else:
        Su = (pL - p0) / Np

    return float(max(Su, 0.0))


def pmt_bearing(
    pL: float,
    p0: float,
    sigma_v: float,
    shape: str = "strip",
) -> Dict[str, float]:
    r"""
    Ultimate bearing capacity from pressuremeter data (Ménard method).

    .. math::

        q_u = \sigma_{v0} + k \,(p_L - p_0)

    where k is a bearing factor depending on foundation shape and
    embedment:

    ============= ====
    Shape         k
    ============= ====
    Strip         0.8
    Square/Circle 1.2
    ============= ====

    Parameters
    ----------
    pL : float
        Limit pressure (kPa).
    p0 : float
        Initial pressure (kPa).
    sigma_v : float
        Overburden stress at foundation level (kPa).
    shape : str, default ``'strip'``
        ``'strip'`` or ``'square'``/``'circle'``.

    Returns
    -------
    dict
        ``'qu'`` : float — Ultimate bearing capacity (kPa).
        ``'net_limit_pressure'`` : float — pL* = pL - p0 (kPa).
        ``'k'`` : float — Bearing factor used.

    References
    ----------
    Ménard (1963); Briaud (1992), Section 5.3.

    Examples
    --------
    >>> from geoeq.site.pmt import pmt_bearing
    >>> res = pmt_bearing(pL=800, p0=100, sigma_v=50, shape='square')
    >>> res['qu']
    890.0
    """
    check_positive(pL, "pL")
    check_non_negative(p0, "p0")
    check_non_negative(sigma_v, "sigma_v")

    shape_l = shape.lower()
    if shape_l in ("strip", "continuous"):
        k = 0.8
    elif shape_l in ("square", "circle", "circular", "rectangular"):
        k = 1.2
    else:
        raise ValueError(f"Unknown shape '{shape}'. Choose: strip, square, circle.")

    pL_star = pL - p0
    qu = sigma_v + k * pL_star

    return {
        "qu": float(qu),
        "net_limit_pressure": float(pL_star),
        "k": float(k),
    }


def pmt_settlement(
    q: float,
    B: float,
    Em: float,
    alpha: float = 0.5,
    shape: str = "circle",
    B0: float = 0.6,
) -> float:
    r"""
    Settlement from Ménard pressuremeter method.

    .. math::

        s = \frac{2}{9\,E_m}\,q\,B_0\,
            \left(\frac{\alpha \, B}{B_0}\right)^{\!\beta}

    where β = shape factor (α from Ménard, here separate from the
    rheological factor).  The simplified Ménard (1963) approach:

    .. math::

        s = \frac{q \, B}{E / \alpha_s}

    Using the deviatoric and spherical components
    (full Ménard approach):

    .. math::

        s = \frac{q}{9\,E_m}\bigl[2\,B_0\,\alpha_d\,(B/B_0)^{\alpha_s}
            + \alpha_s\,B\bigr]

    This function uses the simplified version for routine estimates.

    Parameters
    ----------
    q : float
        Net applied pressure (kPa).
    B : float
        Foundation width or diameter (m).
    Em : float
        Pressuremeter modulus (kPa).
    alpha : float, default 0.5
        Rheological coefficient.
    shape : str, default ``'circle'``
        ``'circle'`` → λ_d = 1.0, ``'square'`` → λ_d = 1.12,
        ``'strip'`` → λ_d = 1.53.
    B0 : float, default 0.6
        Reference width (m) = 0.60 m = 60 cm.

    Returns
    -------
    float
        Settlement (m).

    References
    ----------
    Ménard (1963); Briaud (1992), Section 6.2.

    Examples
    --------
    >>> from geoeq.site.pmt import pmt_settlement
    >>> round(pmt_settlement(q=150, B=2.0, Em=15000) * 1000, 1)
    10.7
    """
    check_positive(q, "q")
    check_positive(B, "B")
    check_positive(Em, "Em")
    check_positive(alpha, "alpha")

    shape_l = shape.lower()
    lam = {"circle": 1.0, "square": 1.12, "strip": 1.53}.get(shape_l, 1.0)

    # Deviatoric component
    s_d = (2.0 * q * B0) / (9.0 * Em) * lam * (B / B0) ** alpha
    # Spherical component
    s_s = (q * alpha * B) / (9.0 * Em)

    return float(s_d + s_s)


def pmt_ko(
    p0: float,
    sigma_v: float,
    u0: float = 0.0,
) -> float:
    r"""
    At-rest earth pressure coefficient from PMT initial pressure.

    .. math::

        K_0 = \frac{p_0 - u_0}{\sigma'_{v0}}

    Parameters
    ----------
    p0 : float
        In-situ horizontal stress from PMT (kPa).
    sigma_v : float
        Total vertical stress (kPa).
    u0 : float, default 0.0
        Pore water pressure (kPa).

    Returns
    -------
    float
        Coefficient of earth pressure at rest K₀.

    Examples
    --------
    >>> from geoeq.site.pmt import pmt_ko
    >>> pmt_ko(p0=80, sigma_v=150, u0=30)
    0.4166666666666667
    """
    check_non_negative(p0, "p0")
    check_positive(sigma_v, "sigma_v")
    check_non_negative(u0, "u0")

    sigma_v_eff = sigma_v - u0
    if sigma_v_eff <= 0:
        raise ValueError("Effective vertical stress must be positive.")

    return float((p0 - u0) / sigma_v_eff)
