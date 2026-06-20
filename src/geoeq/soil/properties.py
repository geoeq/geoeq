r"""
Soil phase-relationship and index-property formulas for GeoEq.

This module provides functions for computing basic soil properties
from index values such as void ratio, specific gravity, water content,
and degree of saturation.

All functions:
- Validate inputs and raise ``ValueError`` for physically impossible values.
- Accept both scalars and numpy arrays (use ``np.asarray`` internally).
- Return plain floats for scalar inputs, numpy arrays for array inputs.

Units
-----
- Unit weight: kN/m³  (default)
- Density: kg/m³
- Void ratio, porosity, saturation, water content: dimensionless
- Specific gravity: dimensionless
- Atterberg limits and indices: percent (%) unless noted

References
----------
Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 3–4.
Holtz, R. D., Kovacs, W. D., & Sheahan, T. C. (2011). *An Introduction to
Geotechnical Engineering*, 2nd ed., Ch. 2.
Skempton, A. W. (1953). The colloidal "activity" of clays. *Proc. 3rd ICSMFE*, 1, 57–61.
"""

import numpy as np

from geoeq.core.constants import GAMMA_WATER, DENSITY_WATER
from geoeq.core.validation import (
    check_positive,
    check_non_negative,
    check_range,
    check_fraction,
)


def _scalar_or_array(result, *inputs):
    """Return a float if all inputs were scalar, else a numpy array."""
    if all(np.ndim(x) == 0 for x in inputs if x is not None):
        return float(result)
    return np.asarray(result, dtype=float)


# -------------------------------------------------------------------
# Phase relations  (Das Ch. 3)
# -------------------------------------------------------------------

def void_ratio(
    n=None,
    w=None,
    Gs=None,
    S=None,
    Vv=None,
    Vs=None,
):
    r"""Compute void ratio from one of several input combinations.

    The void ratio is the ratio of the volume of voids to the volume of
    solids in a soil mass.

    .. math::

        e = \frac{V_v}{V_s}                         \quad\text{[Das Eq. 3.1]}

        e = \frac{n}{1 - n}                          \quad\text{[Das Eq. 3.5]}

        e = \frac{w \cdot G_s}{S}                    \quad\text{[Das Eq. 3.10]}

    Parameters
    ----------
    n : float or array_like, optional
        Porosity, dimensionless, in (0, 1).
    w : float or array_like, optional
        Water content, dimensionless (decimal, not %).
    Gs : float or array_like, optional
        Specific gravity of solids, typically 2.60–2.80.
    S : float or array_like, optional
        Degree of saturation, dimensionless, in (0, 1].
    Vv : float or array_like, optional
        Volume of voids (any consistent unit).
    Vs : float or array_like, optional
        Volume of solids (same unit as *Vv*).

    Returns
    -------
    float or ndarray
        Void ratio (dimensionless, ≥ 0).

    Raises
    ------
    ValueError
        If inputs are physically invalid or insufficient.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Eqs. 3.1, 3.5, 3.10.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.void_ratio(n=0.5)
    1.0

    >>> sp.void_ratio(w=0.25, Gs=2.65, S=1.0)
    0.6625
    """
    if Vv is not None and Vs is not None:
        Vv = np.asarray(Vv, dtype=float)
        Vs = np.asarray(Vs, dtype=float)
        check_non_negative(Vv, "Vv")
        check_positive(Vs, "Vs")
        return _scalar_or_array(Vv / Vs, Vv, Vs)
    elif n is not None:
        n = np.asarray(n, dtype=float)
        check_range(n, "porosity", 0.0, 1.0, inclusive=False)
        return _scalar_or_array(n / (1.0 - n), n)
    elif w is not None and Gs is not None and S is not None:
        w = np.asarray(w, dtype=float)
        Gs = np.asarray(Gs, dtype=float)
        S = np.asarray(S, dtype=float)
        check_non_negative(w, "water_content")
        check_positive(Gs, "Gs")
        check_positive(S, "S")
        check_fraction(S, "degree_of_saturation")
        return _scalar_or_array(w * Gs / S, w, Gs, S)
    else:
        raise ValueError(
            "Insufficient inputs for void_ratio. Provide one of: "
            "(Vv, Vs), (n), or (w, Gs, S)."
        )


def porosity(e=None, Vv=None, V=None):
    r"""Compute porosity from void ratio or volumes.

    .. math::

        n = \frac{V_v}{V} = \frac{e}{1 + e}         \quad\text{[Das Eq. 3.3]}

    Parameters
    ----------
    e : float or array_like, optional
        Void ratio, dimensionless (≥ 0).
    Vv : float or array_like, optional
        Volume of voids (any consistent unit).
    V : float or array_like, optional
        Total volume (same unit as *Vv*).

    Returns
    -------
    float or ndarray
        Porosity, dimensionless, in [0, 1).

    Raises
    ------
    ValueError
        If inputs are physically invalid or insufficient.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Eq. 3.3.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.porosity(e=1.0)
    0.5

    >>> sp.porosity(Vv=0.4, V=1.0)
    0.4
    """
    if Vv is not None and V is not None:
        Vv = np.asarray(Vv, dtype=float)
        V = np.asarray(V, dtype=float)
        check_non_negative(Vv, "Vv")
        check_positive(V, "V")
        return _scalar_or_array(Vv / V, Vv, V)
    elif e is not None:
        e = np.asarray(e, dtype=float)
        check_non_negative(e, "void_ratio")
        return _scalar_or_array(e / (1.0 + e), e)
    else:
        raise ValueError(
            "Insufficient inputs for porosity. Provide (e) or (Vv, V)."
        )


def saturation(w=None, Gs=None, e=None, Vw=None, Vv=None):
    r"""Compute the degree of saturation.

    .. math::

        S = \frac{V_w}{V_v} = \frac{w \cdot G_s}{e}  \quad\text{[Das Eq. 3.6, 3.10]}

    Parameters
    ----------
    w : float or array_like, optional
        Water content, dimensionless (decimal).
    Gs : float or array_like, optional
        Specific gravity of solids.
    e : float or array_like, optional
        Void ratio.
    Vw : float or array_like, optional
        Volume of water.
    Vv : float or array_like, optional
        Volume of voids.

    Returns
    -------
    float or ndarray
        Degree of saturation, dimensionless, in [0, 1].

    Raises
    ------
    ValueError
        If computed saturation exceeds 1.0 (physically impossible) or
        inputs are insufficient.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Eqs. 3.6, 3.10.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.saturation(w=0.25, Gs=2.65, e=0.80)
    0.828125

    >>> sp.saturation(Vw=0.25, Vv=0.50)
    0.5
    """
    if Vw is not None and Vv is not None:
        Vw = np.asarray(Vw, dtype=float)
        Vv = np.asarray(Vv, dtype=float)
        check_non_negative(Vw, "Vw")
        check_positive(Vv, "Vv")
        val = Vw / Vv
        if np.any(val > 1.0 + 1e-9):
            raise ValueError(
                f"Computed saturation > 1.0; check inputs (Vw={Vw}, Vv={Vv})."
            )
        val = np.minimum(val, 1.0)
        return _scalar_or_array(val, Vw, Vv)
    elif w is not None and Gs is not None and e is not None:
        w = np.asarray(w, dtype=float)
        Gs = np.asarray(Gs, dtype=float)
        e = np.asarray(e, dtype=float)
        check_non_negative(w, "water_content")
        check_positive(Gs, "Gs")
        check_positive(e, "void_ratio")
        val = w * Gs / e
        if np.any(val > 1.0 + 1e-9):
            raise ValueError(
                f"Computed saturation > 1.0; check inputs (w={w}, Gs={Gs}, e={e})."
            )
        val = np.minimum(val, 1.0)
        return _scalar_or_array(val, w, Gs, e)
    else:
        raise ValueError(
            "Insufficient inputs for saturation. Provide (Vw, Vv) or (w, Gs, e)."
        )


def water_content(
    S=None, Gs=None, e=None,
    Mw=None, Ms=None,
    Ww=None, Ws=None,
):
    r"""Compute gravimetric water content.

    .. math::

        w = \frac{M_w}{M_s} = \frac{S \cdot e}{G_s}  \quad\text{[Das Eq. 3.7, 3.10]}

    Parameters
    ----------
    S : float or array_like, optional
        Degree of saturation, in [0, 1].
    Gs : float or array_like, optional
        Specific gravity of solids.
    e : float or array_like, optional
        Void ratio.
    Mw : float or array_like, optional
        Mass of water (any consistent unit).
    Ms : float or array_like, optional
        Mass of solids (same unit as *Mw*).
    Ww : float or array_like, optional
        Weight of water (any consistent unit).
    Ws : float or array_like, optional
        Weight of solids (same unit as *Ww*).

    Returns
    -------
    float or ndarray
        Gravimetric water content (dimensionless, ≥ 0).

    Raises
    ------
    ValueError
        If inputs are physically invalid or insufficient.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Eqs. 3.7, 3.10.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.water_content(Mw=20.0, Ms=100.0)
    0.2

    >>> sp.water_content(S=0.80, Gs=2.65, e=0.72)
    0.21735849056603773
    """
    if Mw is not None and Ms is not None:
        Mw = np.asarray(Mw, dtype=float)
        Ms = np.asarray(Ms, dtype=float)
        check_non_negative(Mw, "Mw")
        check_positive(Ms, "Ms")
        return _scalar_or_array(Mw / Ms, Mw, Ms)
    elif Ww is not None and Ws is not None:
        Ww = np.asarray(Ww, dtype=float)
        Ws = np.asarray(Ws, dtype=float)
        check_non_negative(Ww, "Ww")
        check_positive(Ws, "Ws")
        return _scalar_or_array(Ww / Ws, Ww, Ws)
    elif S is not None and Gs is not None and e is not None:
        S = np.asarray(S, dtype=float)
        Gs = np.asarray(Gs, dtype=float)
        e = np.asarray(e, dtype=float)
        check_fraction(S, "degree_of_saturation")
        check_positive(Gs, "Gs")
        check_positive(e, "void_ratio")
        return _scalar_or_array(S * e / Gs, S, Gs, e)
    else:
        raise ValueError(
            "Insufficient inputs for water_content. "
            "Provide (Mw, Ms), (Ww, Ws), or (S, Gs, e)."
        )


def specific_gravity(Ms=None, Vs=None, gamma_w=GAMMA_WATER):
    r"""Compute specific gravity of soil solids.

    .. math::

        G_s = \frac{M_s}{V_s \cdot \rho_w}
            = \frac{\gamma_s}{\gamma_w}              \quad\text{[Das Ch. 3]}

    Parameters
    ----------
    Ms : float or array_like
        Mass of solids (grams).
    Vs : float or array_like
        Volume of solids (cm³).
    gamma_w : float, optional
        Unit weight of water. Not used in mass/volume calculation but
        kept for API consistency. Default 9.81 kN/m³.

    Returns
    -------
    float or ndarray
        Specific gravity (dimensionless, typically 2.60–2.80).

    Raises
    ------
    ValueError
        If *Ms* or *Vs* are not positive.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Ch. 3.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.specific_gravity(Ms=265.0, Vs=100.0)
    2.65

    >>> sp.specific_gravity(Ms=270.0, Vs=100.0)
    2.7
    """
    if Ms is None or Vs is None:
        raise ValueError("Must provide both Ms and Vs for specific_gravity.")
    Ms = np.asarray(Ms, dtype=float)
    Vs = np.asarray(Vs, dtype=float)
    check_positive(Ms, "Ms")
    check_positive(Vs, "Vs")
    # Gs = Ms / (Vs * rho_w), with rho_w = 1.0 g/cm³
    Gs = Ms / (Vs * 1.0)
    return _scalar_or_array(Gs, Ms, Vs)


# -------------------------------------------------------------------
# Density and Unit Weight  (Das Ch. 3)
# -------------------------------------------------------------------

def density(
    Gs=None,
    e=None,
    S=None,
    mass=None,
    volume=None,
    kind="bulk",
    unit="kN/m3",
):
    r"""Compute unit weight or density from phase properties.

    .. math::

        \gamma_d   = \frac{G_s \cdot \gamma_w}{1 + e}
                                                    \quad\text{[Das Eq. 3.13]}

        \gamma     = \frac{(G_s + S e)\,\gamma_w}{1 + e}
                                                    \quad\text{[Das Eq. 3.14]}

        \gamma_{sat} = \frac{(G_s + e)\,\gamma_w}{1 + e}
                                                    \quad\text{[Das Eq. 3.15]}

        \gamma'    = \gamma_{sat} - \gamma_w
                   = \frac{(G_s - 1)\,\gamma_w}{1 + e}
                                                    \quad\text{[Das Eq. 3.16]}

    Parameters
    ----------
    Gs : float or array_like, optional
        Specific gravity of solids (typically 2.60–2.80, allow 1.0–4.0).
    e : float or array_like, optional
        Void ratio (typically 0.3–1.5, allow 0.0–15.0).
    S : float or array_like, optional
        Degree of saturation, in [0, 1]. Required for ``kind='bulk'``.
    mass : float or array_like, optional
        Total mass (kg). Alternative to (Gs, e).
    volume : float or array_like, optional
        Total volume (m³). Alternative to (Gs, e).
    kind : str, optional
        ``'dry'``, ``'saturated'``, ``'submerged'``, ``'bulk'``, or
        ``'all'``.  Default ``'bulk'``.
    unit : str, optional
        ``'kN/m3'``, ``'kg/m3'``, ``'pcf'``, or ``'g/cm3'``.
        Default ``'kN/m3'``.

    Returns
    -------
    float or ndarray or dict
        Unit weight / density value.  If ``kind='all'`` returns a dict
        with keys ``'dry'``, ``'saturated'``, ``'submerged'`` (and
        ``'bulk'`` if *S* is given).

    Raises
    ------
    ValueError
        If inputs are invalid or insufficient.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Eqs. 3.13–3.16.

    Examples
    --------
    >>> from geoeq.soil import properties as sp

    Dry unit weight (Das Example 3.2, Gs=2.67, e=0.72):

    >>> sp.density(Gs=2.67, e=0.72, kind="dry")
    15.237...

    Saturated unit weight:

    >>> sp.density(Gs=2.67, e=0.72, kind="saturated")
    19.33...

    All unit weights at once:

    >>> sp.density(Gs=2.65, e=0.70, S=0.8, kind="all")
    {'dry': ..., 'saturated': ..., 'submerged': ..., 'bulk': ...}
    """
    unit_water = {
        "kN/m3": GAMMA_WATER,
        "kg/m3": DENSITY_WATER,
        "pcf": 62.42796,
        "g/cm3": 1.0,
    }

    mass_vol_factor = {
        "kN/m3": GAMMA_WATER / DENSITY_WATER,
        "kg/m3": 1.0,
        "pcf": 62.42796 / DENSITY_WATER,
        "g/cm3": 1.0 / DENSITY_WATER,
    }

    if unit not in unit_water:
        raise ValueError(
            f"Unsupported unit '{unit}'. Supported: {list(unit_water.keys())}"
        )

    water_prop = unit_water[unit]

    if Gs is not None and e is not None:
        Gs = np.asarray(Gs, dtype=float)
        e = np.asarray(e, dtype=float)
        check_positive(Gs, "specific_gravity")
        check_non_negative(e, "void_ratio")

        calc_dry = Gs * water_prop / (1.0 + e)
        calc_sat = (Gs + e) * water_prop / (1.0 + e)
        calc_sub = (Gs - 1.0) * water_prop / (1.0 + e)

        if kind == "dry":
            return _scalar_or_array(calc_dry, Gs, e)
        elif kind == "saturated":
            return _scalar_or_array(calc_sat, Gs, e)
        elif kind == "submerged":
            return _scalar_or_array(calc_sub, Gs, e)
        elif kind == "bulk":
            if S is None:
                raise ValueError(
                    "Degree of saturation (S) required for kind='bulk'."
                )
            S = np.asarray(S, dtype=float)
            check_fraction(S, "degree_of_saturation")
            result = (Gs + S * e) * water_prop / (1.0 + e)
            return _scalar_or_array(result, Gs, e, S)
        elif kind == "all":
            res = {
                "dry": _scalar_or_array(calc_dry, Gs, e),
                "saturated": _scalar_or_array(calc_sat, Gs, e),
                "submerged": _scalar_or_array(calc_sub, Gs, e),
            }
            if S is not None:
                S = np.asarray(S, dtype=float)
                check_fraction(S, "degree_of_saturation")
                res["bulk"] = _scalar_or_array(
                    (Gs + S * e) * water_prop / (1.0 + e), Gs, e, S
                )
            return res
        else:
            raise ValueError(f"Unknown kind: '{kind}'")

    elif mass is not None and volume is not None:
        mass = np.asarray(mass, dtype=float)
        volume = np.asarray(volume, dtype=float)
        check_positive(mass, "mass")
        check_positive(volume, "volume")
        base_density_kg_m3 = mass / volume
        val = base_density_kg_m3 * mass_vol_factor[unit]
        if kind == "all":
            return {"bulk": _scalar_or_array(val, mass, volume)}
        return _scalar_or_array(val, mass, volume)

    else:
        raise ValueError(
            "Must provide either (Gs, e) or (mass, volume) to compute density."
        )


# -------------------------------------------------------------------
# Relative density  (Das Eq. 3.22)
# -------------------------------------------------------------------

def relative_density(
    e=None, e_max=None, e_min=None,
    rho=None, rho_max=None, rho_min=None,
    kind="void",
):
    r"""Compute relative density (density index).

    .. math::

        D_r = \frac{e_{\max} - e}{e_{\max} - e_{\min}}
                                                    \quad\text{[Das Eq. 3.22]}

        D_r = \frac{\rho_d - \rho_{d,\min}}
                   {\rho_{d,\max} - \rho_{d,\min}}
              \cdot \frac{\rho_{d,\max}}{\rho_d}
                                                    \quad\text{[Das Eq. 3.23]}

    Parameters
    ----------
    e : float or array_like, optional
        Current void ratio.
    e_max : float or array_like, optional
        Maximum void ratio (loosest state).
    e_min : float or array_like, optional
        Minimum void ratio (densest state).
    rho : float or array_like, optional
        Current dry density.
    rho_max : float or array_like, optional
        Maximum dry density (densest).
    rho_min : float or array_like, optional
        Minimum dry density (loosest).
    kind : str, optional
        ``'void'`` for void-ratio method or ``'density'`` for
        dry-density method.  Default ``'void'``.

    Returns
    -------
    float or ndarray
        Relative density, clamped to [0, 1].

    Raises
    ------
    ValueError
        If inputs are invalid or Dr is far outside [0, 1].

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Eqs. 3.22–3.23.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.relative_density(e=0.60, e_max=0.90, e_min=0.50)
    0.75

    >>> sp.relative_density(rho=1600, rho_max=1800, rho_min=1400, kind="density")
    0.5625
    """
    if kind == "void":
        if e is None or e_max is None or e_min is None:
            raise ValueError("Must provide e, e_max, e_min for kind='void'.")
        e = np.asarray(e, dtype=float)
        e_max = np.asarray(e_max, dtype=float)
        e_min = np.asarray(e_min, dtype=float)
        check_non_negative(e, "void_ratio")
        check_non_negative(e_max, "e_max")
        check_non_negative(e_min, "e_min")
        if np.any(e_max <= e_min):
            raise ValueError(
                f"e_max ({e_max}) must be greater than e_min ({e_min})."
            )
        Dr = (e_max - e) / (e_max - e_min)
        if np.any(Dr < -0.01) or np.any(Dr > 1.01):
            raise ValueError(
                f"Computed Dr is outside [0, 1]; check inputs."
            )
        Dr = np.clip(Dr, 0.0, 1.0)
        return _scalar_or_array(Dr, e, e_max, e_min)
    elif kind == "density":
        if rho is None or rho_max is None or rho_min is None:
            raise ValueError(
                "Must provide rho, rho_max, rho_min for kind='density'."
            )
        rho = np.asarray(rho, dtype=float)
        rho_max = np.asarray(rho_max, dtype=float)
        rho_min = np.asarray(rho_min, dtype=float)
        check_positive(rho, "rho")
        check_positive(rho_max, "rho_max")
        check_positive(rho_min, "rho_min")
        if np.any(rho_max <= rho_min):
            raise ValueError("rho_max must be greater than rho_min.")
        Dr = ((rho - rho_min) / (rho_max - rho_min)) * (rho_max / rho)
        if np.any(Dr < -0.01) or np.any(Dr > 1.01):
            raise ValueError(
                f"Computed Dr is outside [0, 1]; check inputs."
            )
        Dr = np.clip(Dr, 0.0, 1.0)
        return _scalar_or_array(Dr, rho, rho_max, rho_min)
    else:
        raise ValueError("kind must be 'void' or 'density'.")


# -------------------------------------------------------------------
# Atterberg limits and consistency  (Das Ch. 4)
# -------------------------------------------------------------------

def atterberg(LL=None, PL=None, w=None, kind="PI"):
    r"""Compute Atterberg-limit indices: PI, LI, or CI.

    .. math::

        PI = LL - PL                                 \quad\text{[Das Ch. 4]}

        LI = \frac{w - PL}{PI}                      \quad\text{[Das Ch. 4]}

        CI = \frac{LL - w}{PI}                       \quad\text{[Das Ch. 4]}

    Parameters
    ----------
    LL : float or array_like, optional
        Liquid limit (%).
    PL : float or array_like, optional
        Plastic limit (%).
    w : float or array_like, optional
        Natural water content (%, same scale as LL/PL).
    kind : str, optional
        ``'PI'``, ``'LI'``, ``'CI'``, or ``'all'``.  Default ``'PI'``.

    Returns
    -------
    float, ndarray, or dict
        Requested index or dict of all indices.

    Raises
    ------
    ValueError
        If PL > LL or insufficient inputs.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Ch. 4.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.atterberg(LL=45, PL=22, kind="PI")
    23.0

    >>> sp.atterberg(LL=48, PL=22, w=35, kind="all")
    {'PI': 26.0, 'LI': 0.5, 'CI': 0.5}
    """
    res = {}

    if LL is not None and PL is not None:
        LL = np.asarray(LL, dtype=float)
        PL = np.asarray(PL, dtype=float)
        check_non_negative(LL, "liquid_limit")
        check_non_negative(PL, "plastic_limit")
        if np.any(PL > LL):
            raise ValueError(
                f"Plastic limit ({PL}) cannot exceed liquid limit ({LL})."
            )
        PI = LL - PL
        res["PI"] = _scalar_or_array(PI, LL, PL)

    if w is not None and "PI" in res:
        w = np.asarray(w, dtype=float)
        check_non_negative(w, "water_content")
        PI = res["PI"]
        res["LI"] = _scalar_or_array((w - PL) / PI, w, PL)
        res["CI"] = _scalar_or_array((LL - w) / PI, w, LL)

    if kind == "all":
        if not res:
            raise ValueError(
                "Insufficient inputs to calculate Atterberg indices."
            )
        return res
    elif kind in res:
        return res[kind]
    else:
        raise ValueError(
            f"Cannot calculate '{kind}' with provided inputs "
            "(need LL, PL, and optionally w)."
        )


def activity(PI, clay_fraction):
    r"""Compute the activity of clay.

    .. math::

        A = \frac{PI}{CF}                            \quad\text{[Skempton, 1953]}

    where *CF* is the clay-size fraction (% finer than 2 µm, expressed
    as a percentage, e.g. 25 for 25 %).

    Parameters
    ----------
    PI : float or array_like
        Plasticity index (%).
    clay_fraction : float or array_like
        Clay-size fraction (%), i.e. % finer than 2 µm.

    Returns
    -------
    float or ndarray
        Activity (dimensionless).

    Raises
    ------
    ValueError
        If PI < 0 or clay_fraction ≤ 0.

    Notes
    -----
    Classification (Skempton 1953):

    - A < 0.75  → Inactive
    - 0.75 ≤ A ≤ 1.25  → Normal
    - A > 1.25  → Active

    References
    ----------
    Skempton, A. W. (1953). The colloidal "activity" of clays. *Proc. 3rd
    ICSMFE*, 1, 57–61.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.activity(PI=30, clay_fraction=40)
    0.75

    >>> sp.activity(PI=50, clay_fraction=20)
    2.5
    """
    PI = np.asarray(PI, dtype=float)
    clay_fraction = np.asarray(clay_fraction, dtype=float)
    check_non_negative(PI, "PI")
    check_positive(clay_fraction, "clay_fraction")
    A = PI / clay_fraction
    return _scalar_or_array(A, PI, clay_fraction)


def sensitivity(Su_undisturbed, Su_remolded):
    r"""Compute the sensitivity of clay.

    .. math::

        S_t = \frac{S_{u(\text{undisturbed})}}{S_{u(\text{remolded})}}
                                                    \quad\text{[Das Ch. 4]}

    Parameters
    ----------
    Su_undisturbed : float or array_like
        Undrained shear strength of undisturbed clay (kPa).
    Su_remolded : float or array_like
        Undrained shear strength of remolded clay (kPa).

    Returns
    -------
    float or ndarray
        Sensitivity (dimensionless, ≥ 1).

    Raises
    ------
    ValueError
        If either strength is not positive.

    Notes
    -----
    Classification (Das Ch. 4):

    - St = 1           → Insensitive
    - 1 < St ≤ 4       → Medium sensitive
    - 4 < St ≤ 8       → Sensitive
    - St > 8           → Extra sensitive / Quick clay

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Ch. 4.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.sensitivity(Su_undisturbed=100, Su_remolded=25)
    4.0

    >>> sp.sensitivity(Su_undisturbed=50, Su_remolded=50)
    1.0
    """
    Su_und = np.asarray(Su_undisturbed, dtype=float)
    Su_rem = np.asarray(Su_remolded, dtype=float)
    check_positive(Su_und, "Su_undisturbed")
    check_positive(Su_rem, "Su_remolded")
    St = Su_und / Su_rem
    return _scalar_or_array(St, Su_und, Su_rem)


def liquidity_index(w, PL, PI):
    r"""Compute the liquidity index (standalone convenience function).

    .. math::

        LI = \frac{w - PL}{PI}                      \quad\text{[Das Ch. 4]}

    Parameters
    ----------
    w : float or array_like
        Natural water content (%).
    PL : float or array_like
        Plastic limit (%).
    PI : float or array_like
        Plasticity index (%).

    Returns
    -------
    float or ndarray
        Liquidity index (dimensionless).

    Raises
    ------
    ValueError
        If PI ≤ 0 (non-plastic soil).

    Notes
    -----
    - LI < 0  → soil is in a semi-solid / solid state
    - 0 ≤ LI ≤ 1  → soil is in the plastic range
    - LI > 1  → soil is in the liquid state

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Ch. 4.

    Examples
    --------
    >>> from geoeq.soil import properties as sp
    >>> sp.liquidity_index(w=35, PL=22, PI=26)
    0.5

    >>> sp.liquidity_index(w=22, PL=22, PI=26)
    0.0
    """
    w = np.asarray(w, dtype=float)
    PL = np.asarray(PL, dtype=float)
    PI = np.asarray(PI, dtype=float)
    check_non_negative(w, "water_content")
    check_non_negative(PL, "plastic_limit")
    check_positive(PI, "PI")
    LI = (w - PL) / PI
    return _scalar_or_array(LI, w, PL, PI)
