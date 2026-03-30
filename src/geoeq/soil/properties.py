"""
Soil phase-relationship formulas for GeoEq.

This module provides functions for computing basic soil properties
from index values such as void ratio, specific gravity, water content,
and degree of saturation.

All functions validate inputs and raise ``ValueError`` for physically
impossible values.

Units
-----
- Unit weight: kN/m^3
- Density: kg/m^3
- Void ratio, porosity, saturation, water content: dimensionless
- Specific gravity: dimensionless
- Atterberg limits and indices: percent (%) unless noted

References
----------
Das, B. M. *Principles of Geotechnical Engineering*, 9th ed., Ch. 3.
Holtz, R. D., Kovacs, W. D., & Sheahan, T. C. *An Introduction to
Geotechnical Engineering*, 2nd ed., Ch. 2.
"""

from geoeq.core.constants import GAMMA_WATER, DENSITY_WATER
from geoeq.core.validation import (
    check_positive,
    check_non_negative,
    check_range,
    check_fraction,
)


# -------------------------------------------------------------------
# Void ratio / porosity
# -------------------------------------------------------------------


def void_ratio_from_porosity(n: float) -> float:
    """
    Compute void ratio from porosity.

    .. math:: e = \\frac{n}{1 - n}

    Parameters
    ----------
    n : float
        Porosity, strictly between 0 and 1.

    Returns
    -------
    float
        Void ratio.

    Raises
    ------
    ValueError
        If *n* is not in the open interval (0, 1).

    Examples
    --------
    >>> void_ratio_from_porosity(0.5)
    1.0
    """
    check_range(n, "porosity", 0.0, 1.0, inclusive=False)
    return n / (1.0 - n)


def porosity_from_void_ratio(e: float) -> float:
    """
    Compute porosity from void ratio.

    .. math:: n = \\frac{e}{1 + e}

    Parameters
    ----------
    e : float
        Void ratio (must be positive).

    Returns
    -------
    float
        Porosity.

    Examples
    --------
    >>> porosity_from_void_ratio(1.0)
    0.5
    """
    check_positive(e, "void_ratio")
    return e / (1.0 + e)


# -------------------------------------------------------------------
# Degree of saturation and water content
# -------------------------------------------------------------------


def degree_of_saturation(w: float, Gs: float, e: float) -> float:
    """
    Compute degree of saturation from water content, Gs, and void ratio.

    .. math:: S = \\frac{w \\, G_s}{e}

    Parameters
    ----------
    w : float
        Gravimetric water content (decimal, e.g. 0.25 for 25 %).
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.

    Returns
    -------
    float
        Degree of saturation (0 to 1).

    Raises
    ------
    ValueError
        If the computed saturation exceeds 1.0 (physically impossible;
        indicates inconsistent inputs).

    Examples
    --------
    >>> degree_of_saturation(w=0.25, Gs=2.65, e=0.80)
    0.828125
    """
    check_non_negative(w, "water_content")
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    S = w * Gs / e
    if S > 1.0 + 1e-9:
        raise ValueError(
            f"Computed saturation {S:.4f} > 1.0; check inputs "
            f"(w={w}, Gs={Gs}, e={e})."
        )
    return min(S, 1.0)


def water_content(S: float, Gs: float, e: float) -> float:
    """
    Compute gravimetric water content from saturation, Gs, and void ratio.

    .. math:: w = \\frac{S \\, e}{G_s}

    Parameters
    ----------
    S : float
        Degree of saturation (0 to 1).
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.

    Returns
    -------
    float
        Water content (decimal).

    Examples
    --------
    >>> water_content(S=0.80, Gs=2.65, e=0.72)
    0.2173...
    """
    check_fraction(S, "degree_of_saturation")
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    return S * e / Gs


# -------------------------------------------------------------------
# Unit weights
# -------------------------------------------------------------------


def dry_unit_weight(
    Gs: float, e: float, gamma_w: float = GAMMA_WATER
) -> float:
    """
    Dry unit weight of soil.

    .. math:: \\gamma_d = \\frac{G_s \\, \\gamma_w}{1 + e}

    Parameters
    ----------
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.
    gamma_w : float, optional
        Unit weight of water (kN/m^3).  Default 9.81.

    Returns
    -------
    float
        Dry unit weight (kN/m^3).

    Examples
    --------
    >>> round(dry_unit_weight(Gs=2.65, e=0.72), 2)
    15.11
    """
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    check_positive(gamma_w, "gamma_w")
    return Gs * gamma_w / (1.0 + e)


def saturated_unit_weight(
    Gs: float, e: float, gamma_w: float = GAMMA_WATER
) -> float:
    """
    Saturated unit weight (degree of saturation = 1).

    .. math:: \\gamma_{sat} = \\frac{(G_s + e) \\, \\gamma_w}{1 + e}

    Parameters
    ----------
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.
    gamma_w : float, optional
        Unit weight of water (kN/m^3).  Default 9.81.

    Returns
    -------
    float
        Saturated unit weight (kN/m^3).
    """
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    return (Gs + e) * gamma_w / (1.0 + e)


def bulk_unit_weight(
    Gs: float, e: float, S: float, gamma_w: float = GAMMA_WATER
) -> float:
    """
    Bulk (moist, total) unit weight.

    .. math:: \\gamma = \\frac{(G_s + S \\, e) \\, \\gamma_w}{1 + e}

    Parameters
    ----------
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.
    S : float
        Degree of saturation (0 to 1).
    gamma_w : float, optional
        Unit weight of water (kN/m^3).  Default 9.81.

    Returns
    -------
    float
        Bulk unit weight (kN/m^3).
    """
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    check_fraction(S, "degree_of_saturation")
    return (Gs + S * e) * gamma_w / (1.0 + e)


def submerged_unit_weight(
    Gs: float, e: float, gamma_w: float = GAMMA_WATER
) -> float:
    """
    Submerged (buoyant, effective) unit weight.

    .. math:: \\gamma' = \\frac{(G_s - 1) \\, \\gamma_w}{1 + e}

    Parameters
    ----------
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.
    gamma_w : float, optional
        Unit weight of water (kN/m^3).  Default 9.81.

    Returns
    -------
    float
        Submerged unit weight (kN/m^3).
    """
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    return (Gs - 1.0) * gamma_w / (1.0 + e)


# -------------------------------------------------------------------
# Density
# -------------------------------------------------------------------


def dry_density(
    Gs: float, e: float, rho_w: float = DENSITY_WATER
) -> float:
    """
    Dry density.

    .. math:: \\rho_d = \\frac{G_s \\, \\rho_w}{1 + e}

    Parameters
    ----------
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.
    rho_w : float, optional
        Density of water (kg/m^3).  Default 1000.

    Returns
    -------
    float
        Dry density (kg/m^3).
    """
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    return Gs * rho_w / (1.0 + e)


def bulk_density(
    Gs: float, e: float, S: float, rho_w: float = DENSITY_WATER
) -> float:
    """
    Bulk density.

    .. math:: \\rho = \\frac{(G_s + S \\, e) \\, \\rho_w}{1 + e}

    Parameters
    ----------
    Gs : float
        Specific gravity of solid particles.
    e : float
        Void ratio.
    S : float
        Degree of saturation (0 to 1).
    rho_w : float, optional
        Density of water (kg/m^3).  Default 1000.

    Returns
    -------
    float
        Bulk density (kg/m^3).
    """
    check_positive(Gs, "specific_gravity")
    check_positive(e, "void_ratio")
    check_fraction(S, "degree_of_saturation")
    return (Gs + S * e) * rho_w / (1.0 + e)


def density(
    Gs: float = None,
    e: float = None,
    S: float = None,
    mass: float = None,
    volume: float = None,
    kind: str = "bulk",
    unit: str = "kN/m3"
):
    """
    Unified function to compute density or unit weight based on index properties
    or mass and volume inputs.

    Parameters
    ----------
    Gs : float, optional
        Specific gravity of solid particles.
    e : float, optional
        Void ratio.
    S : float, optional
        Degree of saturation (0 to 1). Required if kind='bulk' using Gs and e.
    mass : float, optional
        Mass in kg (used if Gs and e are not provided).
    volume : float, optional
        Volume in m^3 (used if Gs and e are not provided).
    kind : str, default "bulk"
        Variant to compute: 'dry', 'saturated', 'submerged', 'bulk', or 'all'.
    unit : str, default "kN/m3"
        Output unit: 'kN/m3', 'kg/m3', 'pcf', or 'g/cm3'.

    Returns
    -------
    float or dict
        Computed unit weight/density or a dictionary of all variants if kind="all".
    """
    unit_water = {
        "kN/m3": GAMMA_WATER,   # 9.81
        "kg/m3": DENSITY_WATER, # 1000.0
        "pcf": 62.42796,        # roughly 1000 kg/m3 in lb/ft3
        "g/cm3": 1.0,           # 1000 kg/m3 in g/cm3
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
        check_positive(Gs, "specific_gravity")
        check_positive(e, "void_ratio")
        
        calc_dry = Gs * water_prop / (1.0 + e)
        calc_sat = (Gs + e) * water_prop / (1.0 + e)
        calc_sub = (Gs - 1.0) * water_prop / (1.0 + e)
        
        if kind == "dry":
            return calc_dry
        elif kind == "saturated":
            return calc_sat
        elif kind == "submerged":
            return calc_sub
        elif kind == "bulk":
            if S is None:
                raise ValueError("Degree of saturation (S) must be provided for 'bulk' if (Gs, e) are used.")
            check_fraction(S, "degree_of_saturation")
            return (Gs + S * e) * water_prop / (1.0 + e)
        elif kind == "all":
            res = {
                "dry": calc_dry,
                "saturated": calc_sat,
                "submerged": calc_sub,
            }
            if S is not None:
                try:
                    check_fraction(S, "degree_of_saturation")
                    res["bulk"] = (Gs + S * e) * water_prop / (1.0 + e)
                except ValueError:
                    pass
            return res
        else:
            raise ValueError(f"Unknown kind: '{kind}'")
            
    elif mass is not None and volume is not None:
        check_positive(mass, "mass")
        check_positive(volume, "volume")
        base_density_kg_m3 = mass / volume
        val = base_density_kg_m3 * mass_vol_factor[unit]
        
        if kind == "all":
            return {kind if kind != "all" else "bulk": val}
        return val
        
    else:
        raise ValueError("Must provide either (Gs, e) or (mass, volume) to compute density.")



# -------------------------------------------------------------------
# Void ratio from other quantities
# -------------------------------------------------------------------


def void_ratio_from_water_content(w: float, Gs: float, S: float) -> float:
    """
    Compute void ratio from water content, specific gravity, and saturation.

    .. math:: e = \\frac{w \\, G_s}{S}

    Parameters
    ----------
    w : float
        Water content (decimal).
    Gs : float
        Specific gravity of solid particles.
    S : float
        Degree of saturation (0 to 1, must be > 0).

    Returns
    -------
    float
        Void ratio.
    """
    check_non_negative(w, "water_content")
    check_positive(Gs, "specific_gravity")
    check_positive(S, "degree_of_saturation")
    check_fraction(S, "degree_of_saturation")
    return w * Gs / S


def void_ratio_from_dry_unit_weight(
    gamma_d: float, Gs: float, gamma_w: float = GAMMA_WATER
) -> float:
    """
    Compute void ratio from dry unit weight and specific gravity.

    .. math:: e = \\frac{G_s \\, \\gamma_w}{\\gamma_d} - 1

    Parameters
    ----------
    gamma_d : float
        Dry unit weight (kN/m^3).
    Gs : float
        Specific gravity of solid particles.
    gamma_w : float, optional
        Unit weight of water (kN/m^3).  Default 9.81.

    Returns
    -------
    float
        Void ratio.
    """
    check_positive(gamma_d, "dry_unit_weight")
    check_positive(Gs, "specific_gravity")
    return Gs * gamma_w / gamma_d - 1.0


# -------------------------------------------------------------------
# Relative density
# -------------------------------------------------------------------


def relative_density(e: float, e_max: float, e_min: float) -> float:
    """
    Relative density (density index) of a granular soil.

    .. math:: D_r = \\frac{e_{max} - e}{e_{max} - e_{min}}

    Parameters
    ----------
    e : float
        In-situ void ratio.
    e_max : float
        Maximum void ratio (loosest state).
    e_min : float
        Minimum void ratio (densest state).

    Returns
    -------
    float
        Relative density (0 to 1).

    Raises
    ------
    ValueError
        If ``e_max <= e_min`` or the result falls well outside [0, 1].

    Examples
    --------
    >>> relative_density(e=0.60, e_max=0.90, e_min=0.50)
    0.75
    """
    check_positive(e, "void_ratio")
    check_positive(e_max, "e_max")
    check_positive(e_min, "e_min")
    if e_max <= e_min:
        raise ValueError(
            f"e_max ({e_max}) must be greater than e_min ({e_min})."
        )
    Dr = (e_max - e) / (e_max - e_min)
    if Dr < -0.01 or Dr > 1.01:
        raise ValueError(
            f"Computed Dr={Dr:.4f} is outside [0, 1]; check inputs."
        )
    return max(0.0, min(1.0, Dr))


# -------------------------------------------------------------------
# Atterberg limits and consistency
# -------------------------------------------------------------------


def plasticity_index(LL: float, PL: float) -> float:
    """
    Plasticity index.

    .. math:: PI = LL - PL

    Parameters
    ----------
    LL : float
        Liquid limit (%).
    PL : float
        Plastic limit (%).

    Returns
    -------
    float
        Plasticity index (%).

    Raises
    ------
    ValueError
        If ``PL > LL``.

    Examples
    --------
    >>> plasticity_index(LL=48, PL=22)
    26
    """
    check_non_negative(LL, "liquid_limit")
    check_non_negative(PL, "plastic_limit")
    if PL > LL:
        raise ValueError(
            f"Plastic limit ({PL}) cannot exceed liquid limit ({LL})."
        )
    return LL - PL


def liquidity_index(w: float, PL: float, PI: float) -> float:
    """
    Liquidity index.

    .. math:: LI = \\frac{w - PL}{PI}

    Parameters
    ----------
    w : float
        Natural water content (%).
    PL : float
        Plastic limit (%).
    PI : float
        Plasticity index (%).

    Returns
    -------
    float
        Liquidity index (dimensionless).

    Examples
    --------
    >>> round(liquidity_index(w=35, PL=22, PI=26), 4)
    0.5
    """
    check_non_negative(w, "water_content")
    check_non_negative(PL, "plastic_limit")
    check_positive(PI, "plasticity_index")
    return (w - PL) / PI


def consistency_index(w: float, LL: float, PI: float) -> float:
    """
    Consistency index.

    .. math:: CI = \\frac{LL - w}{PI}

    Parameters
    ----------
    w : float
        Natural water content (%).
    LL : float
        Liquid limit (%).
    PI : float
        Plasticity index (%).

    Returns
    -------
    float
        Consistency index (dimensionless).

    Examples
    --------
    >>> consistency_index(w=35, LL=48, PI=26)
    0.5
    """
    check_non_negative(w, "water_content")
    check_non_negative(LL, "liquid_limit")
    check_positive(PI, "plasticity_index")
    return (LL - w) / PI


def atterberg(
    LL: float = None,
    PL: float = None,
    w: float = None,
    kind: str = "PI"
):
    """
    Unified function for Atterberg indices.
    
    Parameters
    ----------
    LL : float, optional
        Liquid limit (%).
    PL : float, optional
        Plastic limit (%).
    w : float, optional
        Natural water content (%).
    kind : str, default "PI"
        Variant to compute: 'PI', 'LI', 'CI', or 'all'.
        
    Returns
    -------
    float or dict
        Computed index or dictionary of all indices if kind="all".
    """
    res = {}
    if LL is not None and PL is not None:
        check_non_negative(LL, "liquid_limit")
        check_non_negative(PL, "plastic_limit")
        if PL > LL:
            raise ValueError(
                f"Plastic limit ({PL}) cannot exceed liquid limit ({LL})."
            )
        res["PI"] = LL - PL
        
    if w is not None and "PI" in res:
        check_non_negative(w, "water_content")
        res["LI"] = (w - PL) / res["PI"]
        res["CI"] = (LL - w) / res["PI"]
        
    if kind == "all":
        if not res:
            raise ValueError("Insufficient inputs to calculate Atterberg indices.")
        return res
    elif kind in res:
        return res[kind]
    else:
        raise ValueError(
            f"Cannot calculate '{kind}' with provided inputs (need LL, PL, and optionally w)."
        )

