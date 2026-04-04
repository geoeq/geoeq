"""
Geotechnical grain size distribution plot using matplotlib.
"""

from typing import Union, List, Dict, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from geoeq.soil.grain_size import grain_d10, grain_d30, grain_d60, grain_Cu, grain_Cc

def grain_size_plot(
    data: Union[Dict[str, np.ndarray], Dict[str, Dict[str, np.ndarray]]], 
    smooth: bool = False,
    annotation: bool = False,
    D_para: bool = True, 
    Cu_para: bool = True, 
    Cc_para: bool = True,
    param_pos: Union[str, Tuple[float, float]] = "top right",
    save_as: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    **kwargs
) -> plt.Figure:
    """
    Professional grain size distribution plot with advanced smoothing.
    
    Args:
        data: Dict with 'diameter' and 'percent_finer', or Dict of Dicts for multi-source.
        smooth: If True, uses high-resolution PCHIP interpolation.
        annotation: If True, acknowledge Sieve and Hydrometer parts separately.
        D_para: Show markers/projection for D10, D30, D60.
        Cu_para: Display Cu on plot.
        Cc_para: Display Cc on plot.
        param_pos: Position of parameter box. E.g., 'top right', 'top left', or (x, y).
        save_as: Filename to save.
        ax: matplotlib axes.
        kwargs: matplotlib line properties.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7))
    else:
        fig = ax.get_figure()
    
    # Marker cycle for multi-datasets
    marker_cycle = ['s', '*', '^', 'D', 'o', 'v']
    
    # Handle single vs multi-dataset
    datasets = {}
    if "diameter" in data:
        datasets["Combined"] = data
        combined_data = data
    else:
        datasets = data
        all_d = np.concatenate([ds["diameter"] for ds in datasets.values()])
        all_p = np.concatenate([ds["percent_finer"] for ds in datasets.values()])
        combined_data = {"diameter": all_d, "percent_finer": all_p}
        
    # 1. Plot the continuous smooth line
    if smooth:
        d_all = combined_data["diameter"]
        p_all = combined_data["percent_finer"]
        idx_all = np.argsort(d_all)
        d_s_all = d_all[idx_all]
        p_s_all = p_all[idx_all]
        
        if len(d_s_all) > 2:
            log_d_all = np.log10(np.maximum(d_s_all, 1e-6))
            interp = PchipInterpolator(log_d_all, p_s_all)
            d_smooth = np.logspace(np.log10(d_s_all[0]), np.log10(d_s_all[-1]), 1000)
            p_smooth = interp(np.log10(d_smooth))
            
            line_kwargs = kwargs.copy()
            line_kwargs.pop('marker', None)
            ax.plot(d_smooth, p_smooth, **line_kwargs)
            
    # 2. Plot Markers and Legends for individual datasets
    for i, (label, dset) in enumerate(datasets.items()):
        d = dset["diameter"]
        p = dset["percent_finer"]
        idx = np.argsort(d)
        
        # Determine marker for this dataset
        m = kwargs.get('marker', marker_cycle[i % len(marker_cycle)])
        
        if not smooth:
            ax.plot(d[idx], p[idx], label=label if annotation else None, marker=m, **kwargs)
        else:
            # If smoothing, plot only markers for each part
            marker_kwargs = kwargs.copy()
            marker_kwargs['linestyle'] = 'None'
            marker_kwargs['marker'] = m
            ax.plot(d[idx], p[idx], label=label if annotation else None, **marker_kwargs)
    
    # Global Plot Styling
    ax.invert_xaxis()
    ax.set_xscale('log')
    ax.set_xlabel("Particle Diameter (mm)", fontweight="bold")
    ax.set_ylabel("Percent Passing (%)", fontweight="bold")
    ax.set_ylim(-2, 108)
    ax.set_xlim(100, 0.001)
    
    # Grid and Shading
    ax.grid(True, which="major", linestyle="-", alpha=0.4, color='gray')
    ax.grid(True, which="minor", linestyle="--", alpha=0.2, color='gray')
    ax.axvspan(100, 4.75, color='#e0e0e0', alpha=0.2)  # Gravel
    ax.axvspan(4.75, 0.075, color='#fdf5e6', alpha=0.2) # Sand
    ax.axvspan(0.075, 0.001, color='#e6f3ff', alpha=0.2) # Fines
    ax.axvline(4.75, color='black', alpha=0.3, linewidth=1, linestyle='-')
    ax.axvline(0.075, color='black', alpha=0.3, linewidth=1, linestyle='-')

    # Dx Parameters (Red Dotted / Professional)
    if D_para:
        d10 = grain_d10(combined_data)
        d30 = grain_d30(combined_data)
        d60 = grain_d60(combined_data)
        for val, target in zip([d60, d30, d10], [60, 30, 10]):
            if not np.isnan(val):
                ax.hlines(target, xmin=105.0, xmax=val, colors='red', linestyles=':', linewidth=1.1)
                ax.vlines(val, ymin=-2.0, ymax=target, colors='red', linestyles=':', linewidth=1.1)
                ax.plot(val, target, marker='o', markersize=5, color='red', alpha=0.8)
    
    # Text Box Positioning
    text_parts = []
    if D_para:
        text_parts.append(f"$D_{{60}}$ : {grain_d60(combined_data):.3f} mm")
        text_parts.append(f"$D_{{30}}$ : {grain_d30(combined_data):.3f} mm")
        text_parts.append(f"$D_{{10}}$ : {grain_d10(combined_data):.3f} mm")
    if Cu_para:
        cu = grain_Cu(combined_data)
        text_parts.append(f"$C_u$  : {cu:.2f}")
    if Cc_para:
        cc = grain_Cc(combined_data)
        text_parts.append(f"$C_c$  : {cc:.2f}")
    
    if text_parts:
        # Default positioning logic
        tx, ty, tha, tva = 0.97, 0.95, 'right', 'top' # top right
        if isinstance(param_pos, str):
            pos_map = {
                "top right": (0.97, 0.95, 'right', 'top'),
                "top left": (0.03, 0.95, 'left', 'top'),
                "bottom right": (0.97, 0.05, 'right', 'bottom'),
                "bottom left": (0.03, 0.05, 'left', 'bottom'),
            }
            if param_pos in pos_map:
                tx, ty, tha, tva = pos_map[param_pos]
        elif isinstance(param_pos, (tuple, list)) and len(param_pos) == 2:
            tx, ty = param_pos
            
        ax.text(tx, ty, "\n".join(text_parts), transform=ax.transAxes, 
                fontsize=9, verticalalignment=tva, horizontalalignment=tha,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
    
    # Class labels (Moved slightly to avoid box overlap)
    ax.text(25, 104, "GRAVEL", ha="center", va="center", fontsize=9, fontweight="bold", color="#666666")
    ax.text(0.6, 104, "SAND", ha="center", va="center", fontsize=9, fontweight="bold", color="#666666")
    ax.text(0.012, 104, "FINES", ha="left", va="center", fontsize=9, fontweight="bold", color="#666666")

    if annotation:
        ax.legend(loc="lower left", fontsize=8, framealpha=0.9, edgecolor='gray')

    if save_as:
        plt.savefig(save_as, bbox_inches='tight', dpi=300)
        
    return fig
