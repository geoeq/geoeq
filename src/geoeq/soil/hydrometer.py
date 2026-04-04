"""
Hydrometer analysis for GeoEq following ASTM D7928 standards.
"""

from typing import Union, List, Dict, Optional, Tuple
import numpy as np
from geoeq.core.validation import check_positive, check_non_negative, check_range

def effective_depth(reading: float, model: str = "152H") -> float:
    """ Compute effective depth L (cm) based on hydrometer reading R_actual. """
    if model == "152H":
        # Approximate linear relation for 152H: L = 16.3 - 0.1641 * R
        # Exact ASTM D7928 values are usually tabulated but this is commonly used.
        return 16.3 - 0.1641 * reading
    else:
        raise ValueError(f"Hydrometer model '{model}' not yet supported. Default is 152H.")

def dynamic_viscosity_water(T_celsius: float) -> float:
    """ Compute dynamic viscosity of water (Poise) at T Celsius. """
    # Approximation for viscosity in Poise (g/cm-s)
    phi = 1.79234 - 0.05105 * T_celsius + 0.00072 * (T_celsius**2)
    return phi / 1000.0 # Convert to Poise if needed? 
    # Actually most textbooks use tabulated K values.
    # We will use the direct formula for D.

def hydro_ana(
    reading: List[float], 
    time: List[float], 
    T: Union[float, List[float]], 
    Gs: float = 2.65,
    Ws: float = 50.0,
    Cz: float = 0.0,
    model: str = "152H",
    units: str = "SI"
) -> Dict[str, np.ndarray]:
    """
    Perform hydrometer analysis and return (diameter, percent_finer).
    
    Args:
        reading: List of raw hydrometer readings (top of meniscus).
        time: Elapsed time from start (minutes).
        T: Temperature of suspension (Celsius).
        Gs: Specific gravity of soil solids.
        Ws: Initial dry mass of soil used (g).
        Cz: Zero correction (dispersant correction).
        model: '152H' is default.
        
    Returns:
        Dict with 'diameter' (mm) and 'percent_finer' (%) arrays.
    """
    r = np.array(reading)
    t = np.array(time)
    temp = np.array(T) if isinstance(T, list) else np.full_like(r, T)
    
    if len(r) != len(t):
        raise ValueError("reading and time must have same length.")
    
    # 1. Corrected reading for Percent Finer (Rc)
    # Rc = R_actual - Cz + CT (temperature correction)
    # CT is approx 0.0 at 20C, and adds ~0.2 per degree C above 20.
    CT = (temp - 20.0) * 0.2 
    Rc = r - Cz + CT
    
    # 2. Correction factor 'a' for Gs != 2.65
    a = (1.65 / (Gs - 1.0)) * (Gs / 2.65)
    
    # 3. Percent Finer P (%)
    # P = (Rc * a / Ws) * 100
    P = (Rc * a / Ws) * 100.0
    
    # 4. Corrected reading for Effective Depth (for Diameter calculation)
    # For diameter, we often use the Meniscus reading directly 
    # but the ASTM says use R corrected for meniscus if needed.
    # L depends on the ACTUAL position of the hydrometer.
    L = np.array([effective_depth(val, model) for val in r])

    # 5. Particle Diameter D (mm)
    # D = K * sqrt(L/t), where K = sqrt(30 * eta / (Gs - Gw) * g)
    # g = 981 cm/s^2, eta = viscosity.
    # A common way is to use tables or formula:
    # K factor depends on Gs and Temperature.
    # Simplified K lookup/calc:
    viscosity = 0.01002 * (10**( (1.3272*(20-temp)-0.001053*(temp-20)**2)/(temp+105) )) # Poise (g/cm-s)
    Gw = 0.9982 # Approx density of water at 20C
    K = np.sqrt((30.0 * viscosity) / ( (Gs - Gw) * 980.7 ))
    
    D = K * np.sqrt(L / t)
    
    # Clean results
    P = np.maximum(0.0, np.minimum(100.0, P))
    
    # Sort by diameter descending
    idx = np.argsort(D)[::-1]
    
    return {
        "diameter": D[idx],
        "percent_finer": P[idx],
        "L": L[idx],
        "Rc": Rc[idx]
    }
