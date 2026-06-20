"""
Plot helpers for ge.dynamics:

* ``gmax_curves_plot``   -- G/Gmax and damping vs strain (Darendeli).
* ``liquefaction_chart`` -- CSR vs (N1)60,cs with Youd-2001 and IB-2008
                             triggering curves.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from geoeq.dynamics.modulus import modulus_reduction, damping_ratio
from geoeq.dynamics.liquefaction import liquefaction_crr


# -----------------------------------------------------------------
# G/Gmax + damping vs strain
# -----------------------------------------------------------------
def gmax_curves_plot(
    PI_values=(0, 15, 30, 50, 100),
    sigma_m_eff: float = 100.0, OCR: float = 1.0,
    strain_range=None, method: str = "darendeli",
    ax=None, save_as: str = None,
):
    """Plot G/Gmax (left axis) and damping D (right axis) versus
    shear strain for a family of plasticity indices.

    Reproduces the Vucetic-Dobry style family of curves using
    Darendeli (2001)'s closed-form expressions.
    """
    import matplotlib.pyplot as plt
    if strain_range is None:
        strain_range = np.logspace(-6, -2, 80)  # 0.0001 % .. 1 %
    if ax is None:
        fig, ax1 = plt.subplots(figsize=(8, 5))
    else:
        ax1 = ax
        fig = ax1.figure
    ax2 = ax1.twinx()

    cmap = plt.cm.viridis(np.linspace(0.2, 0.9, len(PI_values)))
    for PI, colour in zip(PI_values, cmap):
        G_ratio = [modulus_reduction(s, sigma_m_eff, PI, OCR, method)
                   for s in strain_range]
        D = [damping_ratio(s, sigma_m_eff, PI, OCR, method=method) * 100
             for s in strain_range]
        ax1.plot(strain_range * 100, G_ratio, color=colour, linewidth=2,
                 label=f"PI = {PI}")
        ax2.plot(strain_range * 100, D, color=colour, linewidth=1.2,
                 linestyle="--", alpha=0.7)

    ax1.set_xscale("log")
    ax1.set_xlim(1e-4, 1.0)
    ax1.set_xlabel(r"Cyclic shear strain $\gamma$ (%)")
    ax1.set_ylabel(r"$G / G_{\max}$ (solid)")
    ax2.set_ylabel(r"Damping ratio $D$ (\%), dashed")
    ax1.set_ylim(0, 1.05)
    ax2.set_ylim(0, 30)
    ax1.set_title(
        f"Modulus reduction + damping — Darendeli (2001), "
        f"$\\sigma'_m={sigma_m_eff}$ kPa, OCR$={OCR}$")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(loc="center left", fontsize=9)
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig


# -----------------------------------------------------------------
# Liquefaction triggering chart -- CSR vs (N1)60,cs
# -----------------------------------------------------------------
def liquefaction_chart(
    data_N: Sequence[float] = None,
    data_CSR: Sequence[float] = None,
    data_FS: Sequence[float] = None,
    Mw: float = 7.5,
    methods=("youd_2001", "idriss_boulanger_2008"),
    ax=None, save_as: str = None,
):
    """Liquefaction-triggering chart: CRR curve(s) plus optional CSR data.

    Parameters
    ----------
    data_N, data_CSR : sequences
        Field data points to overlay (corrected blow counts and CSR).
    data_FS : sequence, optional
        Per-point factor of safety; if given, points are coloured
        red (FS < 1) / green (FS >= 1).
    methods : sequence of str
        CRR curves to draw. Default both Youd-2001 and IB-2008.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    N_arr = np.linspace(1, 30, 200)
    colours = {"youd_2001": "#1f3a93",
               "idriss_boulanger_2008": "#c0392b"}
    labels = {"youd_2001": "NCEER (Youd et al. 2001)",
              "idriss_boulanger_2008": "Idriss \\& Boulanger (2008)"}
    for m in methods:
        crr = []
        for N in N_arr:
            res = liquefaction_crr(N160cs=N, method=m)
            crr.append(res["CRR"])
        ax.plot(N_arr, crr, color=colours.get(m, "k"),
                linewidth=2.0, label=labels.get(m, m))

    # Shaded liquefaction zone
    ax.fill_between(N_arr, 0, [liquefaction_crr(N160cs=n, method="youd_2001")["CRR"]
                                for n in N_arr],
                    color="#c0392b", alpha=0.08, label="Liquefaction zone")

    # Field data overlay
    if data_N is not None and data_CSR is not None:
        if data_FS is not None:
            colors_pts = ["#c0392b" if fs < 1 else "#27ae60" for fs in data_FS]
        else:
            colors_pts = "#000000"
        ax.scatter(data_N, data_CSR, c=colors_pts, s=70,
                    edgecolors="black", linewidths=0.8, zorder=5)

    ax.set_xlim(0, 30)
    ax.set_ylim(0, 0.6)
    ax.set_xlabel(r"Corrected SPT blow count $(N_1)_{60,cs}$")
    ax.set_ylabel(r"CSR or CRR$_{7.5}$")
    ax.set_title(f"Liquefaction triggering chart ($M_w = {Mw}$)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig
