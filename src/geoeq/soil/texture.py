r"""
USDA soil texture classification for GeoEq.

Implements the USDA textural classification (the "texture triangle")
from sand / silt / clay percentages, and a publication-quality ternary
plot of the twelve texture classes.

The classifier uses the exact USDA boundary inequalities (not a
point-in-polygon approximation), so results are exact at class
boundaries.

References
----------
Soil Science Division Staff (2017). *Soil Survey Manual*, USDA
Handbook 18, Ch. 3 — Examination and Description of Soil Profiles.

Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed.,
Ch. 4 (particle-size classification systems).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as _MplPolygon


_SIN60 = np.sin(np.radians(60.0))


def _to_xy(sand, silt, clay):
    """Ternary (sand, silt, clay) [%] -> cartesian (x, y).

    Convention: bottom axis = % sand (right to left), left axis = % clay,
    right axis = % silt.
    """
    x = silt + clay / 2.0
    y = clay * _SIN60
    return x, y


# -------------------------------------------------------------------
# Classification (exact USDA rules)
# -------------------------------------------------------------------

def classify_usda(sand, silt, clay, atol=0.5):
    r"""Classify soil texture per the USDA textural system.

    Uses the exact boundary inequalities from the USDA Soil Survey
    Manual — e.g. *sand* requires :math:`silt + 1.5\,clay < 15`.

    Parameters
    ----------
    sand : float
        Sand fraction (%), 0.05–2.0 mm.
    silt : float
        Silt fraction (%), 0.002–0.05 mm.
    clay : float
        Clay fraction (%), < 0.002 mm.
    atol : float, optional
        Tolerance on ``sand + silt + clay == 100`` (default 0.5 %).

    Returns
    -------
    dict
        ``texture`` — one of the 12 USDA class names;
        ``sand`` / ``silt`` / ``clay`` — the input percentages.

    Raises
    ------
    ValueError
        If any fraction is negative or the fractions do not sum to 100 %.

    References
    ----------
    Soil Science Division Staff (2017). *Soil Survey Manual*, USDA
    Handbook 18, Ch. 3.

    Examples
    --------
    >>> classify_usda(sand=33, silt=33, clay=34)["texture"]
    'Clay Loam'
    >>> classify_usda(sand=92, silt=5, clay=3)["texture"]
    'Sand'
    """
    for name, v in (("sand", sand), ("silt", silt), ("clay", clay)):
        if v < 0:
            raise ValueError(f"{name} must be non-negative, got {v}")
    total = sand + silt + clay
    if not np.isclose(total, 100.0, atol=atol):
        raise ValueError(
            f"sand + silt + clay must equal 100% (got {total:.1f}%)"
        )

    s, si, c = float(sand), float(silt), float(clay)

    if si + 1.5 * c < 15:
        texture = "Sand"
    elif si + 2.0 * c < 30:
        texture = "Loamy Sand"
    elif (7 <= c <= 20 and s > 52) or (c < 7 and si < 50):
        texture = "Sandy Loam"
    elif 7 <= c <= 27 and 28 <= si < 50 and s <= 52:
        texture = "Loam"
    elif si >= 50 and ((12 <= c < 27) or (si < 80 and c < 12)):
        texture = "Silt Loam"
    elif si >= 80 and c < 12:
        texture = "Silt"
    elif 20 <= c < 35 and si < 28 and s > 45:
        texture = "Sandy Clay Loam"
    elif 27 <= c < 40 and 20 < s <= 45:
        texture = "Clay Loam"
    elif 27 <= c < 40 and s <= 20:
        texture = "Silty Clay Loam"
    elif c >= 35 and s > 45:
        texture = "Sandy Clay"
    elif c >= 40 and si >= 40:
        texture = "Silty Clay"
    elif c >= 40:
        texture = "Clay"
    else:  # pragma: no cover — the rule set tiles the whole triangle
        texture = "Unclassified"

    return {"texture": texture, "sand": sand, "silt": silt, "clay": clay}


# -------------------------------------------------------------------
# Region polygons for plotting (sand, silt, clay vertices)
# -------------------------------------------------------------------
# Derived from the same USDA boundaries — a gap-free tiling of the
# triangle used only for drawing; classification never touches these.

_USDA_REGIONS = {
    "Clay":            [(45, 15, 40), (20, 40, 40), (0, 40, 60), (0, 0, 100), (45, 0, 55)],
    "Silty Clay":      [(0, 60, 40), (20, 40, 40), (0, 40, 60)],
    "Sandy Clay":      [(65, 0, 35), (45, 20, 35), (45, 0, 55)],
    "Clay Loam":       [(45, 28, 27), (45, 15, 40), (20, 40, 40), (20, 53, 27)],
    "Silty Clay Loam": [(20, 53, 27), (20, 40, 40), (0, 60, 40), (0, 73, 27)],
    "Sandy Clay Loam": [(52, 28, 20), (80, 0, 20), (65, 0, 35), (45, 20, 35), (45, 28, 27)],
    "Loam":            [(43, 50, 7), (52, 41, 7), (52, 28, 20), (45, 28, 27), (23, 50, 27)],
    "Silt Loam":       [(23, 50, 27), (0, 73, 27), (0, 88, 12), (8, 80, 12),
                        (20, 80, 0), (50, 50, 0), (43, 50, 7)],
    "Silt":            [(20, 80, 0), (8, 80, 12), (0, 88, 12), (0, 100, 0)],
    "Sandy Loam":      [(70, 30, 0), (50, 50, 0), (43, 50, 7), (52, 41, 7),
                        (52, 28, 20), (80, 0, 20), (85, 0, 15)],
    "Loamy Sand":      [(85, 15, 0), (70, 30, 0), (85, 0, 15), (90, 0, 10)],
    "Sand":            [(100, 0, 0), (85, 15, 0), (90, 0, 10)],
}

_REGION_COLORS = {
    "Clay": "#b06a6f", "Silty Clay": "#a3766c", "Sandy Clay": "#c08a72",
    "Clay Loam": "#b57f62", "Silty Clay Loam": "#9c7a5e",
    "Sandy Clay Loam": "#d19a71", "Loam": "#b5643a", "Silt Loam": "#c9a181",
    "Silt": "#e8dccb", "Sandy Loam": "#d9ab84", "Loamy Sand": "#e0bd93",
    "Sand": "#ecdcc3",
}

# Hand-tuned label anchors (sand, silt, clay) and font sizes.
_LABEL_POS = {
    "Clay": (25, 25, 50), "Silty Clay": (7, 47, 46), "Sandy Clay": (52, 6, 42),
    "Clay Loam": (33, 34, 33), "Silty Clay Loam": (10, 57, 33),
    "Sandy Clay Loam": (60, 13, 27), "Loam": (42, 40, 18),
    "Silt Loam": (20, 65, 15), "Silt": (6, 89, 5),
    "Sandy Loam": (65, 25, 10), "Loamy Sand": (81, 13, 6), "Sand": (93, 4, 3),
}
_LABEL_FS = {
    "Silty Clay": 8, "Sandy Clay": 8, "Silt": 8, "Loamy Sand": 7.5,
    "Sand": 8, "Silty Clay Loam": 8.5,
}


def texture_triangle(sand=None, silt=None, clay=None, labels=None,
                     ax=None, save=None):
    r"""Plot the USDA soil texture triangle, optionally with sample points.

    Draws the twelve USDA texture regions as a ternary diagram with 10 %
    gridlines.  Sample points are classified automatically via
    :func:`classify_usda` and annotated with their texture class.

    Parameters
    ----------
    sand : float or array_like, optional
        Sand fraction(s) (%).
    silt : float or array_like, optional
        Silt fraction(s) (%).
    clay : float or array_like, optional
        Clay fraction(s) (%).
    labels : str or list of str, optional
        Label(s) for the sample point(s).  Defaults to the texture class.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  If ``None``, a new figure is created.
    save : str, optional
        File path to save the figure (e.g. ``'triangle.png'``).

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the chart.

    References
    ----------
    Soil Science Division Staff (2017). *Soil Survey Manual*, USDA
    Handbook 18, Ch. 3, Fig. 3-16.

    Examples
    --------
    >>> from geoeq.soil.texture import texture_triangle
    >>> ax = texture_triangle(sand=40, silt=40, clay=20)  # doctest: +SKIP
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9.5, 8.8))

    # region fills + names
    for name, verts in _USDA_REGIONS.items():
        xy = [_to_xy(*v) for v in verts]
        ax.add_patch(_MplPolygon(
            xy, closed=True, facecolor=_REGION_COLORS[name],
            edgecolor="#2b2b2b", linewidth=1.3, alpha=0.93, zorder=1,
        ))
        lx, ly = _to_xy(*_LABEL_POS[name])
        ax.text(lx, ly, name.upper(), ha="center", va="center",
                fontsize=_LABEL_FS.get(name, 9.5), fontweight="bold",
                color="#1a1a1a", zorder=3)

    # 10 % gridlines
    for v in range(10, 100, 10):
        ax.plot(*zip(_to_xy(100 - v, 0, v), _to_xy(0, 100 - v, v)),
                color="#444", lw=0.3, alpha=0.4, zorder=2)
        ax.plot(*zip(_to_xy(100 - v, v, 0), _to_xy(0, v, 100 - v)),
                color="#444", lw=0.3, alpha=0.4, zorder=2)
        ax.plot(*zip(_to_xy(v, 100 - v, 0), _to_xy(v, 0, 100 - v)),
                color="#444", lw=0.3, alpha=0.4, zorder=2)

    # outer frame
    tri = [_to_xy(100, 0, 0), _to_xy(0, 100, 0), _to_xy(0, 0, 100),
           _to_xy(100, 0, 0)]
    ax.plot([p[0] for p in tri], [p[1] for p in tri],
            color="black", lw=2.2, zorder=4)

    # tick labels
    for v in range(0, 101, 10):
        x, y = _to_xy(100 - v, 0, v)              # clay — left edge
        ax.text(x - 2.2, y + 1.2, f"{v}", ha="right", va="center", fontsize=8)
        x, y = _to_xy(0, 100 - v, v)              # silt — right edge
        ax.text(x + 2.2, y + 1.2, f"{100 - v}", ha="left", va="center",
                fontsize=8)
        x, y = _to_xy(v, 100 - v, 0)              # sand — bottom edge
        ax.text(x, y - 3.5, f"{v}", ha="center", va="top", fontsize=8)

    # axis titles + direction arrows
    ax.text(13.5, 45, "PERCENT CLAY", rotation=60, fontsize=11,
            fontweight="bold", ha="center", va="center")
    ax.annotate("", xy=_to_xy(52, 0, 48), xytext=_to_xy(72, 0, 28),
                arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.text(86.5, 45, "PERCENT SILT", rotation=-60, fontsize=11,
            fontweight="bold", ha="center", va="center")
    ax.annotate("", xy=_to_xy(0, 72, 28), xytext=_to_xy(0, 52, 48),
                arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.text(50, -11.5, "PERCENT SAND", fontsize=11, fontweight="bold",
            ha="center", va="top")
    ax.annotate("", xy=(38, -9.5), xytext=(62, -9.5),
                arrowprops=dict(arrowstyle="->", lw=1.4))

    # sample points
    if sand is not None:
        sand_arr = np.atleast_1d(np.asarray(sand, dtype=float))
        silt_arr = np.atleast_1d(np.asarray(silt, dtype=float))
        clay_arr = np.atleast_1d(np.asarray(clay, dtype=float))
        if isinstance(labels, str):
            labels = [labels]

        for i, (sa, si, cl) in enumerate(zip(sand_arr, silt_arr, clay_arr)):
            res = classify_usda(sa, si, cl)
            x, y = _to_xy(sa, si, cl)
            ax.plot(x, y, "o", ms=11, mfc="#ffe600", mec="black",
                    mew=1.7, zorder=6)
            name = labels[i] if labels is not None and i < len(labels) \
                else res["texture"]
            # short radial offset away from the triangle centroid
            cx, cy = 50.0, 100 * _SIN60 / 3
            vx, vy = x - cx, y - cy
            norm = max(np.hypot(vx, vy), 1e-6)
            dx, dy = 16 * vx / norm, 16 * vy / norm
            ha = "left" if dx >= 0 else "right"
            ax.annotate(
                f"{name}: {res['texture']}\n"
                f"(sa {sa:.0f} / si {si:.0f} / cl {cl:.0f})",
                (x, y), xytext=(x + dx, y + dy), fontsize=8.5,
                fontweight="bold", ha=ha, zorder=7,
                arrowprops=dict(arrowstyle="-", lw=0.9, shrinkA=0, shrinkB=6),
                bbox=dict(boxstyle="round,pad=0.35", fc="#fffef0",
                          ec="black", lw=0.8, alpha=0.95),
            )

    ax.set_xlim(-16, 116)
    ax.set_ylim(-17, 94)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("USDA Soil Texture Classification", fontsize=14,
                 fontweight="bold", pad=14)

    if save:
        ax.figure.savefig(save, dpi=300, bbox_inches="tight")

    return ax
