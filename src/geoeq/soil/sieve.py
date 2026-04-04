"""
Sieve analysis and standard sieve mapping for GeoEq.
"""

from typing import Union, List, Dict, Optional
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative

# Standard sieve mappings (Designation -> Opening in mm)
SIEVE_STANDARDS = {
    "ASTM": {
        "3\"": 75.0,
        "1.5\"": 37.5,
        "3/4\"": 19.0,
        "3/8\"": 9.5,
        "#4": 4.75,
        "#8": 2.36,
        "#10": 2.0,
        "#16": 1.18,
        "#20": 0.85,
        "#30": 0.6,
        "#40": 0.425,
        "#50": 0.3,
        "#60": 0.25,
        "#100": 0.15,
        "#200": 0.075,
    },
    "BS": {
        "75mm": 75.0,
        "37.5mm": 37.5,
        "20mm": 20.0,
        "10mm": 10.0,
        "5mm": 5.0,
        "3.35mm": 3.35,
        "2mm": 2.0,
        "1.18mm": 1.18,
        "600um": 0.6,
        "425um": 0.425,
        "300um": 0.3,
        "212um": 0.212,
        "150um": 0.15,
        "63um": 0.063,
    },
    "IS": {
        "80mm": 80.0,
        "40mm": 40.0,
        "20mm": 20.0,
        "10mm": 10.0,
        "4.75mm": 4.75,
        "2.36mm": 2.36,
        "2mm": 2.0,
        "1.18mm": 1.18,
        "1mm": 1.0,
        "600um": 0.6,
        "425um": 0.425,
        "300um": 0.3,
        "150um": 0.15,
        "75um": 0.075,
    }
}

def get_opening(designation: Union[str, float], standard: str = "ASTM") -> float:
    """Helper to convert designation to mm."""
    if isinstance(designation, (int, float)):
        return float(designation)
    
    std = standard.upper()
    if std not in SIEVE_STANDARDS:
        raise ValueError(f"Unknown standard: {standard}. Use ASTM, BS, or IS.")
    
    if designation in SIEVE_STANDARDS[std]:
        return SIEVE_STANDARDS[std][designation]
    
    # Try cleaning the string if it's not found exactly
    clean = designation.strip().replace(" ", "")
    if clean in SIEVE_STANDARDS[std]:
        return SIEVE_STANDARDS[std][clean]
    
    try:
        return float(designation)
    except ValueError:
        raise ValueError(f"Unknown sieve designation '{designation}' for standard {standard}.")

def sieve_ana(
    opening: List[Union[str, float]], 
    mass_retained: List[float], 
    standard: str = "ASTM",
    total_mass: Optional[float] = None
) -> Dict[str, np.ndarray]:
    """
    Perform sieve analysis and return percent finer.
    
    Args:
        opening: List of sieve designations (e.g., "#4") or nominal openings (mm).
        mass_retained: Mass of soil retained on each sieve (g).
        standard: 'ASTM', 'BS', or 'IS'.
        total_mass: If None, calculated as sum(mass_retained).
        
    Returns:
        Dict containing arrays for opening, mass_retained, percent_retained, 
        cumulative_retained, and percent_finer.
    """
    if len(opening) != len(mass_retained):
        raise ValueError("Opening and mass_retained must have the same length.")
    
    # Convert designations to mm
    mm_openings = np.array([get_opening(o, standard) for o in opening])
    m_ret = np.array(mass_retained)
    
    # Sort by opening (descending)
    idx = np.argsort(mm_openings)[::-1]
    mm_openings = mm_openings[idx]
    m_ret = m_ret[idx]
    
    if total_mass is None:
        total_mass = np.sum(m_ret)
    
    check_positive(total_mass, "total_mass")
    
    pct_ret = (m_ret / total_mass) * 100.0
    cum_ret = np.cumsum(pct_ret)
    pct_finer = 100.0 - cum_ret
    
    # Ensure no negative values (numerical precision)
    pct_finer = np.maximum(pct_finer, 0.0)
    
    return {
        "diameter": mm_openings,
        "opening": mm_openings,
        "mass_retained": m_ret,
        "percent_retained": pct_ret,
        "cumulative_retained": cum_ret,
        "percent_finer": pct_finer,
        "total_mass": total_mass
    }
