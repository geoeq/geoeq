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
# Phase relations
# -------------------------------------------------------------------

def void_ratio(
    n: float = None,
    w: float = None,
    Gs: float = None,
    S: float = None,
    Vv: float = None,
    Vs: float = None
):
    """
    Unified function to compute void ratio.
    """
    if Vv is not None and Vs is not None:
        check_non_negative(Vv, "Vv")
        check_positive(Vs, "Vs")
        return Vv / Vs
    elif n is not None:
        check_range(n, "porosity", 0.0, 1.0, inclusive=False)
        return n / (1.0 - n)
    elif w is not None and Gs is not None and S is not None:
        check_non_negative(w, "water_content")
        check_positive(Gs, "Gs")
        check_positive(S, "S")
        check_fraction(S, "degree_of_saturation")
        return w * Gs / S
    else:
        raise ValueError("Insufficient inputs for void_ratio.")


def porosity(
    e: float = None,
    Vv: float = None,
    V: float = None
):
    """
    Unified function to compute porosity.
    """
    if Vv is not None and V is not None:
        check_non_negative(Vv, "Vv")
        check_positive(V, "V")
        return Vv / V
    elif e is not None:
        check_positive(e, "void_ratio")
        return e / (1.0 + e)
    else:
        raise ValueError("Insufficient inputs for porosity.")


def saturation(
    w: float = None,
    Gs: float = None,
    e: float = None,
    Vw: float = None,
    Vv: float = None
):
    """
    Unified function to compute degree of saturation.
    """
    if Vw is not None and Vv is not None:
        check_non_negative(Vw, "Vw")
        check_positive(Vv, "Vv")
        val = Vw / Vv
        if val > 1.0 + 1e-9:
            raise ValueError(f"Computed saturation {val:.4f} > 1.0; check inputs.")
        return min(val, 1.0)
    elif w is not None and Gs is not None and e is not None:
        check_non_negative(w, "water_content")
        check_positive(Gs, "Gs")
        check_positive(e, "e")
        val = w * Gs / e
        if val > 1.0 + 1e-9:
            raise ValueError(f"Computed saturation {val:.4f} > 1.0; check inputs.")
        return min(val, 1.0)
    else:
        raise ValueError("Insufficient inputs for saturation.")


def water_content(
    S: float = None,
    Gs: float = None,
    e: float = None,
    Mw: float = None,
    Ms: float = None,
    Ww: float = None,
    Ws: float = None
):
    """
    Unified function to compute gravimetric water content.
    """
    if Mw is not None and Ms is not None:
        check_non_negative(Mw, "Mw")
        check_positive(Ms, "Ms")
        return Mw / Ms
    elif Ww is not None and Ws is not None:
        check_non_negative(Ww, "Ww")
        check_positive(Ws, "Ws")
        return Ww / Ws
    elif S is not None and Gs is not None and e is not None:
        check_fraction(S, "degree_of_saturation")
        check_positive(Gs, "Gs")
        check_positive(e, "e")
        return S * e / Gs
    else:
        raise ValueError("Insufficient inputs for water_content.")


# -------------------------------------------------------------------
# Density and Unit Weight
# -------------------------------------------------------------------

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
        raise ValueError(f"Unsupported unit '{unit}'. Supported: {list(unit_water.keys())}")

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
# Relative density
# -------------------------------------------------------------------

def relative_density(
    e: float = None,
    e_max: float = None,
    e_min: float = None,
    rho: float = None,
    rho_max: float = None,
    rho_min: float = None,
    kind: str = "void"
):
    """
    Unified function for relative density.
    """
    if kind == "void":
        if e is None or e_max is None or e_min is None:
            raise ValueError("Must provide e, e_max, and e_min for kind='void'.")
        check_positive(e, "void_ratio")
        check_positive(e_max, "e_max")
        check_positive(e_min, "e_min")
        if e_max <= e_min:
            raise ValueError(f"e_max ({e_max}) must be greater than e_min ({e_min}).")
        Dr = (e_max - e) / (e_max - e_min)
        if Dr < -0.01 or Dr > 1.01:
            raise ValueError(f"Computed Dr={Dr:.4f} is outside [0, 1].")
        return max(0.0, min(1.0, Dr))
    elif kind == "density":
        if rho is None or rho_max is None or rho_min is None:
            raise ValueError("Must provide rho, rho_max, rho_min for kind='density'.")
        check_positive(rho, "rho")
        check_positive(rho_max, "rho_max")
        check_positive(rho_min, "rho_min")
        if rho_max <= rho_min:
            raise ValueError("rho_max must be greater than rho_min.")
        Dr = ((rho - rho_min) / (rho_max - rho_min)) * (rho_max / rho)
        if Dr < -0.01 or Dr > 1.01:
            raise ValueError(f"Computed Dr={Dr:.4f} is outside [0, 1].")
        return max(0.0, min(1.0, Dr))
    else:
        raise ValueError("kind must be 'void' or 'density'.")


# -------------------------------------------------------------------
# Atterberg limits and consistency
# -------------------------------------------------------------------

def atterberg(
    LL: float = None,
    PL: float = None,
    w: float = None,
    kind: str = "PI"
):
    """
    Unified function for Atterberg indices.
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
