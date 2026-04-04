"""
Grain size distribution interpolation and parameter getters.
"""

from typing import Union, List, Dict, Optional
import numpy as np
from scipy.interpolate import interp1d

def grain_interpolate(diameters: np.ndarray, percent_finer: np.ndarray, p_target: float) -> float:
    """ Log-linear interpolation for finding Dx. """
    # filter out zeros/negatives for log
    mask = (diameters > 0) & (percent_finer >= 0)
    d = diameters[mask]
    p = percent_finer[mask]
    
    if len(d) < 2:
        return np.nan
    
    # Sort by percent_finer (ascending for interp)
    idx = np.argsort(p)
    p_sorted = p[idx]
    d_sorted = d[idx]
    
    # Check if p_target is in range
    if p_target < p_sorted[0] or p_target > p_sorted[-1]:
        return np.nan
        
    # Interpolate in log-space (log(D) vs P)
    f = interp1d(p_sorted, np.log10(d_sorted), kind='linear')
    log_d = f(p_target)
    return 10.0**float(log_d)

def grain_d10(data: Dict[str, np.ndarray]) -> float:
    """ Get D10 (effective size). """
    return grain_interpolate(data["diameter"], data["percent_finer"], 10.0)

def grain_d30(data: Dict[str, np.ndarray]) -> float:
    """ Get D30. """
    return grain_interpolate(data["diameter"], data["percent_finer"], 30.0)

def grain_d60(data: Dict[str, np.ndarray]) -> float:
    """ Get D60. """
    return grain_interpolate(data["diameter"], data["percent_finer"], 60.0)

def grain_Cu(data: Dict[str, np.ndarray]) -> float:
    """ Get Uniformity Coefficient Cu = D60/D10. """
    d60 = grain_d60(data)
    d10 = grain_d10(data)
    if d10 > 0 and not np.isnan(d60) and not np.isnan(d10):
        return d60 / d10
    return np.nan

def grain_Cc(data: Dict[str, np.ndarray]) -> float:
    """ Get Coefficient of Curvature Cc = (D30^2)/(D60 * D10). """
    d60 = grain_d60(data)
    d30 = grain_d30(data)
    d10 = grain_d10(data)
    if d60 > 0 and d10 > 0 and not np.isnan(d60) and not np.isnan(d30) and not np.isnan(d10):
        return (d30**2) / (d60 * d10)
    return np.nan

def process_combined_data(sieve_data: Optional[Dict] = None, hydro_data: Optional[Dict] = None) -> Dict[str, np.ndarray]:
    """ Combine sieve and hydrometer results into one dataset. """
    d_list = []
    p_list = []
    
    if sieve_data:
        d_list.append(sieve_data["opening"])
        p_list.append(sieve_data["percent_finer"])
    
    if hydro_data:
        d_list.append(hydro_data["diameter"])
        p_list.append(hydro_data["percent_finer"])
        
    if not d_list:
        raise ValueError("Must provide at least sieve or hydrometer data.")
        
    all_d = np.concatenate(d_list)
    all_p = np.concatenate(p_list)
    
    # Sort and remove duplicates if needed? 
    # Usually we leave both and sort by D descending.
    idx = np.argsort(all_d)[::-1]
    
    return {
        "diameter": all_d[idx],
        "percent_finer": all_p[idx]
    }
