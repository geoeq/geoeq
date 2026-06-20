r"""
Soil classification systems for GeoEq.

Implements the Unified Soil Classification System (USCS) per ASTM D2487
and the AASHTO classification system per AASHTO M145.  Also provides a
Casagrande plasticity chart plotter.

References
----------
ASTM D2487-17 (2017). *Standard Practice for Classification of Soils for
Engineering Purposes (Unified Soil Classification System).*

AASHTO M145-91 (2012). *Standard Specification for Classification of Soils
and Soil-Aggregate Mixtures for Highway Construction Purposes.*

Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
Ch. 4.
"""

import numpy as np
import matplotlib.pyplot as plt


# -------------------------------------------------------------------
# USCS Classification  (ASTM D2487)
# -------------------------------------------------------------------

# USCS group names keyed by symbol
_USCS_NAMES = {
    "GW": "Well-graded gravel",
    "GP": "Poorly graded gravel",
    "GM": "Silty gravel",
    "GC": "Clayey gravel",
    "GW-GM": "Well-graded gravel with silt",
    "GW-GC": "Well-graded gravel with clay",
    "GP-GM": "Poorly graded gravel with silt",
    "GP-GC": "Poorly graded gravel with clay",
    "GM-GC": "Silty, clayey gravel",
    "SW": "Well-graded sand",
    "SP": "Poorly graded sand",
    "SM": "Silty sand",
    "SC": "Clayey sand",
    "SW-SM": "Well-graded sand with silt",
    "SW-SC": "Well-graded sand with clay",
    "SP-SM": "Poorly graded sand with silt",
    "SP-SC": "Poorly graded sand with clay",
    "SM-SC": "Silty, clayey sand",
    "CL": "Lean clay",
    "ML": "Silt",
    "OL": "Organic clay / organic silt",
    "CH": "Fat clay",
    "MH": "Elastic silt",
    "OH": "Organic clay / organic silt",
    "CL-ML": "Silty clay",
    "Pt": "Peat",
}


def _a_line(LL):
    """Casagrande A-line: PI = 0.73 * (LL - 20)."""
    return 0.73 * (LL - 20.0)


def _classify_coarse(gravel, sand, fines, LL, PL, Cu, Cc):
    """Classify coarse-grained soil (fines < 50%)."""
    PI = (LL - PL) if (LL is not None and PL is not None) else None

    if gravel > sand:
        prefix = "G"
    else:
        prefix = "S"

    if fines < 5:
        if Cu is None or Cc is None:
            raise ValueError(
                "Cu and Cc are required when fines < 5% to distinguish "
                "well-graded from poorly graded."
            )
        if prefix == "G":
            well = (Cu >= 4) and (1 <= Cc <= 3)
        else:
            well = (Cu >= 6) and (1 <= Cc <= 3)
        return (prefix + "W") if well else (prefix + "P")

    elif fines > 12:
        if LL is None or PL is None:
            raise ValueError(
                "LL and PL are required when fines > 12%."
            )
        a_val = _a_line(LL)
        if PI >= 4 and PI >= a_val:
            return prefix + "C"
        elif PI < 4 or PI < a_val:
            return prefix + "M"
        else:
            return prefix + "M"

    else:
        # 5 <= fines <= 12  → dual symbol
        if Cu is None or Cc is None:
            raise ValueError(
                "Cu and Cc required when 5% <= fines <= 12%."
            )
        if prefix == "G":
            well = (Cu >= 4) and (1 <= Cc <= 3)
        else:
            well = (Cu >= 6) and (1 <= Cc <= 3)
        grade = "W" if well else "P"

        if LL is not None and PL is not None:
            a_val = _a_line(LL)
            if PI >= 4 and PI >= a_val:
                fine_sym = "C"
            elif PI < 4 or PI < a_val:
                fine_sym = "M"
            else:
                fine_sym = "M"
        else:
            fine_sym = "M"

        return f"{prefix}{grade}-{prefix}{fine_sym}"


def _classify_fine(LL, PL, organic=False):
    """Classify fine-grained soil (fines >= 50%)."""
    if LL is None or PL is None:
        raise ValueError("LL and PL are required for fine-grained classification.")
    PI = LL - PL
    a_val = _a_line(LL)

    if LL < 50:
        if organic:
            return "OL"
        if PI >= 4 and PI >= a_val:
            return "CL"
        elif PI < 4 and PI < a_val:
            return "ML"
        else:
            return "CL-ML"
    else:
        if organic:
            return "OH"
        if PI >= a_val:
            return "CH"
        else:
            return "MH"


def classify_uscs(
    LL=None,
    PL=None,
    gravel=0.0,
    sand=0.0,
    fines=None,
    Cu=None,
    Cc=None,
    organic=False,
):
    r"""Classify a soil using the Unified Soil Classification System (USCS).

    Implements the flowchart in ASTM D2487-17 and Das (2021), Fig. 4.8.

    Parameters
    ----------
    LL : float, optional
        Liquid limit (%).
    PL : float, optional
        Plastic limit (%).
    gravel : float, optional
        Gravel fraction (% retained on No. 4 / 4.75 mm sieve). Default 0.
    sand : float, optional
        Sand fraction (% passing No. 4 and retained on No. 200 / 0.075 mm).
        Default 0.
    fines : float, optional
        Fines content (% passing No. 200).  If not given, computed as
        ``100 - gravel - sand``.
    Cu : float, optional
        Uniformity coefficient.  Required when fines < 12 %.
    Cc : float, optional
        Coefficient of curvature.  Required when fines < 12 %.
    organic : bool, optional
        ``True`` if oven-dried LL is < 75% of non-dried LL (organic soil
        indicator per ASTM D2487). Default ``False``.

    Returns
    -------
    dict
        ``{'symbol': str, 'name': str}`` — USCS group symbol and name.

    Raises
    ------
    ValueError
        If required inputs for the classification path are missing, or
        grain-size fractions do not sum to ≈ 100 %.

    References
    ----------
    ASTM D2487-17 (2017). *Standard Practice for Classification of Soils
    for Engineering Purposes.*

    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Ch. 4, Fig. 4.8.

    Examples
    --------
    >>> from geoeq.soil.classification import classify_uscs

    Well-graded sand:

    >>> classify_uscs(gravel=5, sand=90, fines=5, Cu=8, Cc=2)
    {'symbol': 'SW', 'name': 'Well-graded sand'}

    Fat clay:

    >>> classify_uscs(LL=70, PL=30, fines=80, gravel=5, sand=15)
    {'symbol': 'CH', 'name': 'Fat clay'}
    """
    gravel = float(gravel)
    sand = float(sand)

    if fines is None:
        fines = 100.0 - gravel - sand
    else:
        fines = float(fines)

    total = gravel + sand + fines
    if abs(total - 100.0) > 1.0:
        raise ValueError(
            f"Grain-size fractions must sum to ~100%, got {total:.1f}% "
            f"(gravel={gravel}, sand={sand}, fines={fines})."
        )

    if LL is not None and PL is not None:
        PI = LL - PL
    else:
        PI = None

    if fines >= 50:
        symbol = _classify_fine(LL, PL, organic=organic)
    else:
        symbol = _classify_coarse(gravel, sand, fines, LL, PL, Cu, Cc)

    name = _USCS_NAMES.get(symbol, symbol)
    return {"symbol": symbol, "name": name}


# -------------------------------------------------------------------
# AASHTO Classification  (AASHTO M145)
# -------------------------------------------------------------------

def _group_index(LL, PL, fines):
    """Compute AASHTO Group Index: GI = (F-35)[0.2+0.005(LL-40)] + 0.01(F-15)(PI-10)."""
    PI = LL - PL
    F = fines

    a = max(F - 35, 0)
    b = max(LL - 40, 0)
    c = max(F - 15, 0)
    d = max(PI - 10, 0)

    GI = a * (0.2 + 0.005 * b) + 0.01 * c * d
    return max(0, round(GI))


def classify_aashto(LL, PL, gravel=0.0, sand=0.0, fines=None):
    r"""Classify a soil using the AASHTO system (M145).

    Follows the left-to-right elimination procedure of Das (2021),
    Table 4.4.

    Parameters
    ----------
    LL : float
        Liquid limit (%).
    PL : float
        Plastic limit (%).
    gravel : float, optional
        % retained on No. 10 (2.0 mm) sieve.  Default 0.
    sand : float, optional
        % passing No. 10 and retained on No. 200.  Default 0.
    fines : float, optional
        % passing No. 200.  If not given, computed as
        ``100 - gravel - sand``.

    Returns
    -------
    dict
        ``{'group': str, 'group_index': int, 'description': str}``

    Raises
    ------
    ValueError
        If LL or PL is missing.

    References
    ----------
    AASHTO M145-91 (2012). *Standard Specification for Classification of
    Soils and Soil-Aggregate Mixtures for Highway Construction Purposes.*

    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Ch. 4, Table 4.4.

    Examples
    --------
    >>> from geoeq.soil.classification import classify_aashto
    >>> classify_aashto(LL=45, PL=22, fines=60)
    {'group': 'A-7-6', 'group_index': 13, 'description': 'Clayey soils'}
    """
    if LL is None or PL is None:
        raise ValueError("LL and PL are required for AASHTO classification.")

    LL = float(LL)
    PL = float(PL)
    gravel = float(gravel)
    sand = float(sand)
    if fines is None:
        fines = 100.0 - gravel - sand
    else:
        fines = float(fines)

    PI = LL - PL
    GI = _group_index(LL, PL, fines)

    # Granular materials (fines <= 35%)
    if fines <= 35:
        if fines <= 15 and LL <= 40 and PI <= 6:
            return {"group": "A-1-a", "group_index": 0,
                    "description": "Stone fragments, gravel, and sand"}
        if fines <= 25 and LL <= 40 and PI <= 6:
            return {"group": "A-1-b", "group_index": 0,
                    "description": "Stone fragments, gravel, and sand"}
        if fines <= 10 and PI == 0:
            return {"group": "A-3", "group_index": 0,
                    "description": "Fine sand"}
        if fines <= 35 and LL <= 40 and PI <= 10:
            return {"group": "A-2-4", "group_index": GI,
                    "description": "Silty or clayey gravel and sand"}
        if fines <= 35 and LL > 40 and PI <= 10:
            return {"group": "A-2-5", "group_index": GI,
                    "description": "Silty or clayey gravel and sand"}
        if fines <= 35 and LL <= 40 and PI > 10:
            return {"group": "A-2-6", "group_index": GI,
                    "description": "Silty or clayey gravel and sand"}
        if fines <= 35 and LL > 40 and PI > 10:
            return {"group": "A-2-7", "group_index": GI,
                    "description": "Silty or clayey gravel and sand"}

    # Silt-clay materials (fines > 35%)
    if LL <= 40 and PI <= 10:
        return {"group": "A-4", "group_index": GI,
                "description": "Silty soils"}
    if LL > 40 and PI <= 10:
        return {"group": "A-5", "group_index": GI,
                "description": "Silty soils"}
    if LL <= 40 and PI > 10:
        return {"group": "A-6", "group_index": GI,
                "description": "Clayey soils"}
    # A-7
    if PI <= (LL - 30):
        return {"group": "A-7-5", "group_index": GI,
                "description": "Clayey soils"}
    else:
        return {"group": "A-7-6", "group_index": GI,
                "description": "Clayey soils"}


# -------------------------------------------------------------------
# Casagrande plasticity chart  (Das Fig. 4.8)
# -------------------------------------------------------------------

def plasticity_chart(
    LL=None,
    PL=None,
    labels=None,
    ax=None,
    save=None,
):
    r"""Plot the Casagrande plasticity chart with A-line and U-line.

    Plots the Casagrande A-line and U-line and (optionally) plots soil
    samples on the chart for visual USCS classification.

    .. math::

        \text{A-line:} \quad PI = 0.73\,(LL - 20)   \quad\text{[Das Fig. 4.8]}

        \text{U-line:} \quad PI = 0.9\,(LL - 8)     \quad\text{[Das Fig. 4.8]}

    Parameters
    ----------
    LL : float or array_like, optional
        Liquid limit(s) (%) for data points to plot.
    PL : float or array_like, optional
        Plastic limit(s) (%) for data points to plot.
    labels : str or list of str, optional
        Labels for each data point.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  If ``None``, a new figure is created.
    save : str, optional
        File path to save the figure (e.g. ``'chart.png'``).

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the chart.

    References
    ----------
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
    Ch. 4, Fig. 4.8.

    Examples
    --------
    >>> from geoeq.soil.classification import plasticity_chart
    >>> ax = plasticity_chart(LL=[35, 60], PL=[18, 25])  # doctest: +SKIP
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    ll_range = np.linspace(0, 120, 300)

    # A-line
    a_pi = 0.73 * (ll_range - 20.0)
    a_pi = np.maximum(a_pi, 0.0)
    ax.plot(ll_range, a_pi, "k-", linewidth=1.5, label="A-line")

    # U-line
    u_pi = 0.9 * (ll_range - 8.0)
    u_pi = np.maximum(u_pi, 0.0)
    ax.plot(ll_range, u_pi, "k--", linewidth=1.2, label="U-line")

    # LL = 50 divider
    ax.axvline(x=50, color="gray", linestyle=":", linewidth=0.8)

    # Zone labels
    ax.text(25, 5, "CL-ML", fontsize=9, ha="center", style="italic")
    ax.text(30, 18, "CL", fontsize=11, ha="center", fontweight="bold")
    ax.text(25, 3, "ML or OL", fontsize=9, ha="center", color="gray")
    ax.text(70, 20, "MH or OH", fontsize=9, ha="center", color="gray")
    ax.text(70, 45, "CH or OH", fontsize=11, ha="center", fontweight="bold")

    # Plot data points
    if LL is not None and PL is not None:
        LL_arr = np.atleast_1d(np.asarray(LL, dtype=float))
        PL_arr = np.atleast_1d(np.asarray(PL, dtype=float))
        PI_arr = LL_arr - PL_arr
        ax.scatter(LL_arr, PI_arr, c="red", s=60, zorder=5, edgecolors="black")

        if labels is not None:
            if isinstance(labels, str):
                labels = [labels]
            for i, lbl in enumerate(labels):
                if i < len(LL_arr):
                    ax.annotate(
                        lbl,
                        (LL_arr[i], PI_arr[i]),
                        textcoords="offset points",
                        xytext=(6, 6),
                        fontsize=8,
                    )

    ax.set_xlabel("Liquid Limit, LL (%)", fontsize=12)
    ax.set_ylabel("Plasticity Index, PI (%)", fontsize=12)
    ax.set_title("Casagrande Plasticity Chart", fontsize=13, fontweight="bold")
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 70)
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_aspect("auto")

    if save:
        ax.figure.savefig(save, dpi=300, bbox_inches="tight")

    return ax
