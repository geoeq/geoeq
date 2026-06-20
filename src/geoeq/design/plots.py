"""
Design-module plot helpers.

* ``bearing_capacity_plot`` -- ultimate bearing capacity vs. footing width.
* ``settlement_time_plot`` -- U vs Tv curve (Terzaghi 1-D consolidation).
* ``stress_isobar_plot``  -- pressure bulb beneath a point load.
* ``taylor_chart_plot``   -- stability-number chart for cohesive slopes.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from geoeq.design.bearing import bearing_capacity, bearing_factors
from geoeq.design.boussinesq import stress_bulb
from geoeq.design.settlement import time_factor, consolidation_degree
from geoeq.design.slopes import taylor_stability


# -----------------------------------------------------------------
# Bearing capacity vs. footing width (designer's tool)
# -----------------------------------------------------------------
def bearing_capacity_plot(
    c: float, gamma: float, Df: float, phi: float,
    B_range: Sequence[float] = None, methods: Sequence[str] = None,
    L_over_B: float = 1.0, ax=None, save_as: str = None,
):
    """Plot ultimate bearing capacity vs. footing width B.

    Useful design-chart helper: shows how q_u scales with B and how
    the four classical methods diverge for granular soils.

    Parameters
    ----------
    c, gamma, Df, phi : float
        Soil parameters (kPa, kN/m^3, m, degrees).
    B_range : sequence
        Footing widths to evaluate (m). Default 0.5 .. 5 m.
    methods : sequence of str
        Subset of ``{"terzaghi", "meyerhof", "hansen", "vesic"}``.
        Default: all four.
    L_over_B : float
        Footing aspect ratio (1.0 = square). Use 1e6 for strip.
    """
    import matplotlib.pyplot as plt
    if B_range is None:
        B_range = np.linspace(0.5, 5.0, 30)
    if methods is None:
        methods = ["terzaghi", "meyerhof", "hansen", "vesic"]
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    colours = {"terzaghi": "#888888", "meyerhof": "#1f3a93",
               "hansen": "#27ae60", "vesic": "#c0392b"}
    for m in methods:
        q = []
        for B in B_range:
            L = B * L_over_B if np.isfinite(L_over_B) else None
            res = bearing_capacity(c, gamma, Df, B, phi, L=L, method=m)
            q.append(res["q_u"])
        ax.plot(B_range, q, label=m.title(),
                color=colours.get(m, "k"), linewidth=1.8)
    ax.set_xlabel("Footing width $B$ (m)")
    ax.set_ylabel("Ultimate bearing capacity $q_u$ (kPa)")
    ax.set_title(
        f"$q_u$ vs $B$ — $c={c}$, $\\phi={phi}°$, $D_f={Df}$ m, $\\gamma={gamma}$ kN/m³")
    ax.grid(alpha=0.3)
    ax.legend()
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig


# -----------------------------------------------------------------
# U vs Tv -- Terzaghi 1-D consolidation curve
# -----------------------------------------------------------------
def settlement_time_plot(ax=None, save_as: str = None):
    """The classical Terzaghi consolidation curve: U vs Tv (log scale).

    Plots the approximate piecewise formula used by ``time_factor``:
    Tv = (pi/4) U^2 for U < 0.6, Tv = 1.781 - 0.933 log10(100(1-U))
    for U >= 0.6.
    """
    import matplotlib.pyplot as plt
    U_arr = np.linspace(0.001, 0.99, 200)
    Tv_arr = np.array([time_factor(u) for u in U_arr])
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    ax.plot(Tv_arr, U_arr * 100, color="#1f3a93", linewidth=2.0)
    # Annotated reference points.
    for U in (0.5, 0.9):
        Tv = time_factor(U)
        ax.plot([Tv], [U * 100], "o", color="#c0392b", markersize=6)
        ax.annotate(f"$U={int(U*100)}\\%$\n$T_v={Tv:.3f}$",
                    xy=(Tv, U * 100), xytext=(10, -10),
                    textcoords="offset points", fontsize=9,
                    color="#c0392b")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 3)
    ax.set_ylim(0, 100)
    ax.set_xlabel(r"Time factor $T_v$")
    ax.set_ylabel(r"Average degree of consolidation $U$ (%)")
    ax.set_title("Terzaghi 1-D consolidation: $U$ vs $T_v$")
    ax.grid(True, which="both", alpha=0.3)
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig


# -----------------------------------------------------------------
# Pressure-bulb contour beneath a point load
# -----------------------------------------------------------------
def stress_isobar_plot(P: float = 100.0, z_max: float = 8.0,
                       r_max: float = 4.0, levels=None,
                       ax=None, save_as: str = None):
    """Stress isobars (pressure bulb) beneath a point load.

    Parameters
    ----------
    P : float
        Point load (kN).
    z_max, r_max : float
        Plot extents (m).
    levels : sequence
        Iso-stress contour levels (kPa).
    """
    import matplotlib.pyplot as plt
    res = stress_bulb(P=P, z_range=(0.05, z_max),
                       r_range=(-r_max, r_max), n=120)
    if levels is None:
        # Geometric spacing from max stress in the grid.
        smax = np.nanmax(res["sigma"])
        levels = np.geomspace(max(smax * 0.01, 0.1), smax * 0.9, 8)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))
    else:
        fig = ax.figure
    cs = ax.contourf(res["r"], res["z"], res["sigma"],
                     levels=levels, cmap="RdYlBu_r", alpha=0.8)
    ax.contour(res["r"], res["z"], res["sigma"],
               levels=levels, colors="black", linewidths=0.5)
    cbar = fig.colorbar(cs, ax=ax, label=r"$\Delta\sigma_z$ (kPa)")
    # Mark the load application point.
    ax.plot([0], [0], "v", color="black", markersize=14)
    ax.text(0, -0.2, f"  $P = {P:.0f}$ kN", ha="center", va="bottom",
            fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel(r"Radial distance $r$ (m)")
    ax.set_ylabel(r"Depth $z$ (m)")
    ax.set_title("Boussinesq pressure bulb beneath a point load")
    ax.set_aspect("equal")
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig


# -----------------------------------------------------------------
# Taylor stability number chart
# -----------------------------------------------------------------
def taylor_chart_plot(phi_values=(0, 5, 10, 15, 20, 25),
                     beta_range=None, ax=None, save_as: str = None):
    """Taylor's stability number chart: m vs slope angle beta for
    a series of phi values."""
    import matplotlib.pyplot as plt
    if beta_range is None:
        beta_range = np.linspace(15, 90, 30)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure
    for phi in phi_values:
        m_values = []
        for beta in beta_range:
            # We bypass FS by reading the m value directly.
            res = taylor_stability(phi=phi, c=1, gamma=1, H=1, beta=beta)
            m_values.append(res["m"])
        ax.plot(beta_range, m_values, label=f"$\\phi = {phi}°$",
                linewidth=1.7)
    ax.set_xlabel(r"Slope angle $\beta$ (°)")
    ax.set_ylabel(r"Stability number $m = c / (\gamma H_{cr})$")
    ax.set_title("Taylor's stability number (toe circles)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    if save_as:
        fig.savefig(save_as, dpi=300, bbox_inches="tight")
    return fig
