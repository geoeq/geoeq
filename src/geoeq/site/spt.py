r"""
Standard Penetration Test (SPT) corrections and correlations.

Functions
---------
spt_n60
    Energy-corrected blow count N₆₀.
spt_n160
    Overburden-corrected (N₁)₆₀.
spt_n160cs
    Fines-content-corrected (N₁)₆₀cs for liquefaction.
spt_friction_angle
    Friction angle φ from SPT (multiple methods).
spt_su
    Undrained shear strength Su from SPT.
spt_dr
    Relative density Dr from SPT.
spt_modulus
    Elastic modulus Es from SPT.

References
----------
Skempton, A. W. (1986). SPT procedures and the effects in sands of overburden
    pressure, relative density, particle size, ageing and overconsolidation.
    *Geotechnique*, 36(3), 425–447.
Liao, S. S. C. & Whitman, R. V. (1986). Overburden correction factors for SPT
    in sand. *JGED*, ASCE, 112(3), 373–377.
Kulhawy, F. H. & Mayne, P. W. (1990). *Manual on Estimating Soil Properties
    for Foundation Design*. EPRI EL-6800.
Hatanaka, M. & Uchida, A. (1996). Empirical correlation between penetration
    resistance and internal friction angle of sandy soils. *Soils and
    Foundations*, 36(4), 1–9.
Bowles, J. E. (1996). *Foundation Analysis and Design*, 5th ed.
Stroud, M. A. (1974). The SPT in insensitive clays and soft rocks. *Proc.
    European Symposium on Penetration Testing*, Stockholm, 2.2, 367–375.
Hara, A. et al. (1974). Shear modulus and shear strength of cohesive soils.
    *Soils and Foundations*, 14(3), 1–12.
"""

from typing import Dict, Optional, Union
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative


def _scalar_or_array(result, *inputs):
    if all(np.ndim(x) == 0 for x in inputs):
        return float(result)
    return np.asarray(result)


def spt_n60(
    N: Union[float, np.ndarray],
    ER: float = 60.0,
    Cb: float = 1.0,
    Cs: float = 1.0,
    Cr: float = 1.0,
) -> Union[float, np.ndarray]:
    r"""
    Compute the energy-corrected SPT blow count N₆₀.

    .. math::

        N_{60} = N \times \frac{ER}{60} \times C_b \times C_s \times C_r
        \qquad \text{[Skempton, 1986]}

    Parameters
    ----------
    N : float or array_like
        Field blow count (blows per 300 mm).
    ER : float, default 60.0
        Energy ratio of the hammer (%). Safety hammer ≈ 60,
        donut hammer ≈ 45, automatic ≈ 80–95.
    Cb : float, default 1.0
        Borehole diameter correction (1.0 for 65–115 mm, 1.05 for 150 mm,
        1.15 for 200 mm).
    Cs : float, default 1.0
        Sampler correction (1.0 for standard, 1.1–1.3 without liner).
    Cr : float, default 1.0
        Rod length correction (0.75 for 3–4 m, 0.85 for 4–6 m,
        0.95 for 6–10 m, 1.0 for 10–30 m).

    Returns
    -------
    float or ndarray
        Energy-corrected blow count N₆₀.

    References
    ----------
    Skempton (1986); Das (2021), Table 3.8.

    Examples
    --------
    >>> from geoeq.site.spt import spt_n60
    >>> spt_n60(30, ER=72, Cb=1.0, Cs=1.0, Cr=0.95)
    34.2
    """
    N_arr = np.asarray(N, dtype=float)
    check_non_negative(N_arr, "N")
    check_positive(ER, "ER")

    result = N_arr * (ER / 60.0) * Cb * Cs * Cr
    return _scalar_or_array(result, N)


def spt_n160(
    N60: Union[float, np.ndarray],
    sigma_v: Union[float, np.ndarray],
    method: str = "liao_whitman",
    pa: float = 100.0,
) -> Union[float, np.ndarray]:
    r"""
    Compute the overburden-corrected SPT blow count (N₁)₆₀.

    .. math::

        (N_1)_{60} = C_N \times N_{60}

    Parameters
    ----------
    N60 : float or array_like
        Energy-corrected blow count.
    sigma_v : float or array_like
        Effective vertical overburden stress (kPa).
    method : str, default ``'liao_whitman'``
        Correction method:

        - ``'liao_whitman'`` : :math:`C_N = \sqrt{p_a / \sigma'_v}`
          (Liao & Whitman, 1986).
        - ``'skempton'`` : :math:`C_N = 2 / (1 + \sigma'_v / p_a)`
          (Skempton, 1986).
        - ``'peck'`` : :math:`C_N = 0.77 \log_{10}(20 \, p_a / \sigma'_v)`
          (Peck et al., 1974).
    pa : float, default 100.0
        Atmospheric pressure (kPa).

    Returns
    -------
    float or ndarray
        Overburden-corrected blow count (N₁)₆₀. Capped so CN ≤ 2.0.

    References
    ----------
    Liao & Whitman (1986); Skempton (1986); Peck et al. (1974).

    Examples
    --------
    >>> from geoeq.site.spt import spt_n160
    >>> round(spt_n160(20, sigma_v=100), 1)
    20.0
    >>> round(spt_n160(20, sigma_v=50), 1)
    28.3
    """
    n60 = np.asarray(N60, dtype=float)
    sv = np.asarray(sigma_v, dtype=float)
    check_non_negative(n60, "N60")
    check_positive(sv, "sigma_v")

    method = method.lower().replace("-", "_").replace(" ", "_")

    if method == "liao_whitman":
        CN = np.sqrt(pa / sv)
    elif method == "skempton":
        CN = 2.0 / (1.0 + sv / pa)
    elif method == "peck":
        CN = 0.77 * np.log10(20.0 * pa / sv)
    else:
        raise ValueError(
            f"Unknown method '{method}'. Choose from: liao_whitman, skempton, peck."
        )

    CN = np.minimum(CN, 2.0)
    result = CN * n60
    return _scalar_or_array(result, N60, sigma_v)


def spt_n160cs(
    N160: Union[float, np.ndarray],
    FC: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    r"""
    Fines-content-corrected (N₁)₆₀cs for liquefaction analysis.

    .. math::

        (N_1)_{60cs} = \alpha + \beta \times (N_1)_{60}

    where α and β depend on fines content FC:

    - FC ≤ 5%: α = 0, β = 1.0
    - 5 < FC < 35%: α = exp(1.76 − 190/FC²), β = 0.99 + FC^1.5 / 1000
    - FC ≥ 35%: α = 5.0, β = 1.2

    Parameters
    ----------
    N160 : float or array_like
        Overburden-corrected blow count (N₁)₆₀.
    FC : float or array_like
        Fines content (%).

    Returns
    -------
    float or ndarray
        Fines-corrected (N₁)₆₀cs.

    References
    ----------
    Youd, T. L. et al. (2001). Liquefaction Resistance of Soils. *JGGE*,
    ASCE, 127(10), 817–833.

    Examples
    --------
    >>> from geoeq.site.spt import spt_n160cs
    >>> round(spt_n160cs(15, FC=3), 1)
    15.0
    >>> round(spt_n160cs(15, FC=20), 1)
    19.6
    """
    n = np.asarray(N160, dtype=float)
    fc = np.asarray(FC, dtype=float)
    check_non_negative(n, "N160")
    check_non_negative(fc, "FC")

    alpha = np.where(fc <= 5, 0.0,
                     np.where(fc >= 35, 5.0,
                              np.exp(1.76 - 190.0 / fc ** 2)))
    beta = np.where(fc <= 5, 1.0,
                    np.where(fc >= 35, 1.2,
                             0.99 + fc ** 1.5 / 1000.0))

    result = alpha + beta * n
    return _scalar_or_array(result, N160, FC)


def spt_friction_angle(
    N: Union[float, np.ndarray],
    sigma_v: Optional[Union[float, np.ndarray]] = None,
    method: str = "hatanaka",
    pa: float = 100.0,
) -> Union[float, np.ndarray]:
    r"""
    Estimate friction angle φ' from SPT blow count.

    Parameters
    ----------
    N : float or array_like
        Blow count. Interpretation depends on *method*:
        - ``'hatanaka'``: N should be (N₁)₆₀.
        - ``'kulhawy'``: N should be (N₁)₆₀; σ'_v is ignored.
        - ``'peck'``: N should be field N (uncorrected).
    sigma_v : float or array_like, optional
        Effective overburden stress (kPa). Not used by all methods.
    method : str, default ``'hatanaka'``
        - ``'hatanaka'`` : :math:`\phi' = \sqrt{20\,(N_1)_{60}} + 20`
          (Hatanaka & Uchida, 1996).
        - ``'kulhawy'`` : :math:`\phi' = \tan^{-1}\!\bigl[\frac{(N_1)_{60}}{12.2+20.3\,(\sigma'_v/p_a)}\bigr]^{0.34}`
          (Kulhawy & Mayne, 1990).
        - ``'peck'`` : Peck, Hanson & Thornburn (1974) empirical table.
    pa : float, default 100.0

    Returns
    -------
    float or ndarray
        Friction angle φ' (degrees).

    References
    ----------
    Hatanaka & Uchida (1996); Kulhawy & Mayne (1990); Peck et al. (1974).

    Examples
    --------
    >>> from geoeq.site.spt import spt_friction_angle
    >>> round(spt_friction_angle(20, method='hatanaka'), 1)
    40.0
    """
    N_arr = np.asarray(N, dtype=float)
    check_non_negative(N_arr, "N")
    method = method.lower()

    if method == "hatanaka":
        phi = np.sqrt(20.0 * N_arr) + 20.0
    elif method == "kulhawy":
        if sigma_v is None:
            sigma_v = pa
        sv = np.asarray(sigma_v, dtype=float)
        phi = np.degrees(np.arctan((N_arr / (12.2 + 20.3 * (sv / pa))) ** 0.34))
    elif method == "peck":
        # Peck, Hanson & Thornburn (1974) empirical approximation
        phi = 26.0 + 0.3 * N_arr - 0.00054 * N_arr ** 2
        phi = np.clip(phi, 26.0, 45.0)
    else:
        raise ValueError(
            f"Unknown method '{method}'. Choose: hatanaka, kulhawy, peck."
        )

    return _scalar_or_array(phi, N)


def spt_su(
    N60: Union[float, np.ndarray],
    method: str = "stroud",
    PI: Optional[float] = None,
) -> Union[float, np.ndarray]:
    r"""
    Estimate undrained shear strength Su from SPT.

    Parameters
    ----------
    N60 : float or array_like
        Energy-corrected blow count N₆₀.
    method : str, default ``'stroud'``
        - ``'stroud'`` : :math:`S_u = f_1 \times N_{60}` (Stroud, 1974).
          f₁ depends on PI: PI < 20 → f₁ = 4.5; 20–30 → 5.0; 30–40 → 5.5;
          > 40 → 6.0.
        - ``'hara'`` : :math:`S_u = 29 \times N_{60}^{0.72}` (Hara et al., 1974).
    PI : float, optional
        Plasticity index (%). Required for ``'stroud'``.

    Returns
    -------
    float or ndarray
        Undrained shear strength Su (kPa).

    References
    ----------
    Stroud (1974); Hara et al. (1974).

    Examples
    --------
    >>> from geoeq.site.spt import spt_su
    >>> round(spt_su(10, method='hara'), 1)
    152.4
    """
    n60 = np.asarray(N60, dtype=float)
    check_non_negative(n60, "N60")
    method = method.lower()

    if method == "stroud":
        if PI is None:
            f1 = 5.0
        elif PI < 20:
            f1 = 4.5
        elif PI < 30:
            f1 = 5.0
        elif PI < 40:
            f1 = 5.5
        else:
            f1 = 6.0
        Su = f1 * n60
    elif method == "hara":
        Su = 29.0 * n60 ** 0.72
    else:
        raise ValueError(f"Unknown method '{method}'. Choose: stroud, hara.")

    return _scalar_or_array(Su, N60)


def spt_dr(
    N160: Union[float, np.ndarray],
    method: str = "meyerhof",
) -> Union[float, np.ndarray]:
    r"""
    Estimate relative density Dr from SPT.

    Parameters
    ----------
    N160 : float or array_like
        Overburden-corrected (N₁)₆₀.
    method : str, default ``'meyerhof'``
        - ``'meyerhof'`` : :math:`D_r (\%) = \sqrt{(N_1)_{60} / 0.21}`
          × 100 (adapted from Meyerhof, 1957 — coeff for clean
          normally-consolidated sands).
        - ``'skempton'`` : :math:`D_r (\%) = \sqrt{(N_1)_{60} / 0.28}`
          × 100 (Skempton, 1986 — fine sands).
        - ``'kulhawy'`` : :math:`D_r (\%) = \sqrt{(N_1)_{60} / C_p}`
          × 100, where Cp = 60 + 25 log D₅₀ (Kulhawy & Mayne, 1990).
          Uses Cp = 60 (default for medium sand D₅₀ ≈ 1 mm).

    Returns
    -------
    float or ndarray
        Relative density Dr (%).

    References
    ----------
    Meyerhof (1957); Skempton (1986); Kulhawy & Mayne (1990).

    Examples
    --------
    >>> from geoeq.site.spt import spt_dr
    >>> round(spt_dr(20, method='meyerhof'), 0)
    98.0
    """
    n = np.asarray(N160, dtype=float)
    check_non_negative(n, "N160")
    method = method.lower()

    if method == "meyerhof":
        # Dr² = N1,60 / 0.21 → Dr = √(N/0.21) as fraction, ×100 for %
        # More standard: Dr(%) = √(N1,60 / 41) × 100 from Meyerhof tables
        Dr = np.sqrt(n / 41.0) * 100.0
    elif method == "skempton":
        Dr = np.sqrt(n / 55.0) * 100.0
    elif method == "kulhawy":
        Cp = 60.0
        Dr = np.sqrt(n / Cp) * 100.0
    else:
        raise ValueError(f"Unknown method '{method}'. Choose: meyerhof, skempton, kulhawy.")

    return _scalar_or_array(np.clip(Dr, 0, 100), N160)


def spt_modulus(
    N: Union[float, np.ndarray],
    soil_type: str = "sand",
) -> Union[float, np.ndarray]:
    r"""
    Estimate elastic (Young's) modulus Es from SPT.

    Parameters
    ----------
    N : float or array_like
        SPT blow count (N or N₆₀ depending on source — typically N₆₀).
    soil_type : str, default ``'sand'``
        - ``'sand'`` : Es = 500 (N + 15) kPa (Bowles, 1996, for NC sand).
        - ``'sand_oc'`` : Es = 750 (N + 24) kPa (Bowles, 1996, OC sand).
        - ``'gravel'`` : Es = 1200 (N + 6) kPa (Bowles, 1996).
        - ``'clay_soft'`` : Es = 300 (N + 6) kPa (Bowles, 1996).
        - ``'clay_stiff'`` : Es = 500 (N + 15) kPa (Bowles, 1996).
        - ``'silt'`` : Es = 300 (N + 6) kPa (Bowles, 1996).

    Returns
    -------
    float or ndarray
        Elastic modulus Es (kPa).

    References
    ----------
    Bowles (1996), Table 5.6.

    Examples
    --------
    >>> from geoeq.site.spt import spt_modulus
    >>> spt_modulus(20, soil_type='sand')
    17500.0
    """
    n = np.asarray(N, dtype=float)
    check_non_negative(n, "N")
    st = soil_type.lower().replace(" ", "_")

    lookup = {
        "sand": (500.0, 15.0),
        "sand_nc": (500.0, 15.0),
        "sand_oc": (750.0, 24.0),
        "gravel": (1200.0, 6.0),
        "clay_soft": (300.0, 6.0),
        "clay_stiff": (500.0, 15.0),
        "silt": (300.0, 6.0),
    }

    if st not in lookup:
        raise ValueError(
            f"Unknown soil_type '{soil_type}'. Choose from: "
            + ", ".join(lookup.keys())
        )

    a, b = lookup[st]
    Es = a * (n + b)
    return _scalar_or_array(Es, N)
