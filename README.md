# GeoEq

**A clean, validated Python library for geotechnical engineering calculations.**

The installable package name is **`geoeq`** (like **NumPy** / `numpy`). **GeoEq** is the project name used in documentation and publications.

GeoEq brings common geotechnical formulas from laboratory practice and textbooks into one consistent, validated, and well-documented Python package. It is built for engineers, students, and researchers who want to replace repetitive spreadsheet work with clean, reusable code.

```bash
pip install geoeq
```

Requires Python 3.9+. Dependencies: `numpy`, `matplotlib`, `scipy`.

---

## 🧪 What GeoEq Can Do

### Ch. 2 — Particle Size Analysis
| Function | Description |
| :--- | :--- |
| `sp.sieve_ana()` | Sieve analysis with ASTM, BS, and IS standard mapping |
| `sp.hydro_ana()` | Full ASTM D7928 hydrometer analysis (Stokes' Law) |
| `sp.grain_size_plot()` | Publication-ready semi-log grain size distribution plots |
| `sp.grain_d10/d30/d60()` | Log-linear interpolated Dₓ values in mm |
| `sp.grain_Cu()` | Uniformity Coefficient (D₆₀ / D₁₀) |
| `sp.grain_Cc()` | Coefficient of Curvature (D₃₀² / D₁₀·D₆₀) |

### Ch. 3 — Soil Properties
| Function | Description |
| :--- | :--- |
| `sp.void_ratio()` | Void ratio from porosity, or volumes |
| `sp.porosity()` | Porosity from void ratio, or volumes |
| `sp.density()` | Dry, saturated, bulk, submerged unit weight — any unit system |
| `sp.saturation()` | Degree of saturation from phase params or volumes |
| `sp.water_content()` | Water content from masses or phase parameters |
| `sp.relative_density()` | Dr from void ratio limits or dry density limits |
| `sp.atterberg()` | Plasticity Index, Liquidity Index, Consistency Index |

---

## 📖 Quick Start

```python
import geoeq as sp

# --- Ch. 2: Particle Size Analysis ---

# 1. Sieve Analysis (ASTM, BS, or IS designations)
res = sp.sieve_ana(
    opening=["3/4\"", "#4", "#10", "#40", "#200"],
    mass_retained=[0, 40, 60, 105, 56],
    standard="ASTM",
    total_mass=300.0
)
print(res["percent_finer"])  # array of % passing values

# 2. Hydrometer Analysis (ASTM D7928)
h_res = sp.hydro_ana(
    reading=[52, 48, 44, 36, 25],
    time=[1, 2, 5, 15, 60],    # minutes
    T=22.0,                    # temperature °C
    Gs=2.65,
    Ws=50.0                    # dry mass of soil (g)
)
print(h_res["diameter"])      # equivalent particle diameters in mm

# 3. Publication-Ready Grain Size Plot
sp.grain_size_plot(
    {"Sieve": res, "Hydrometer": h_res},
    smooth=True,           # 1000-point PCHIP interpolation
    annotation=True,       # label each data source in legend
    D_para=True,           # draw D10, D30, D60 projection lines
    Cu_para=True,
    Cc_para=True,
    param_pos="top right", # or "top left", "bottom right", (x, y)
    color="black",
    marker="o",
    save_as="grain_size.pdf"  # PNG, SVG, PDF, EPS all supported
)

# 4. Extract Parameters
d10 = sp.grain_d10(res)
cu  = sp.grain_Cu(res)
cc  = sp.grain_Cc(res)
print(f"D10={d10:.3f} mm  Cu={cu:.2f}  Cc={cc:.2f}")
```

```python
# --- Ch. 3: Soil Properties ---

# Void ratio & porosity
e = sp.void_ratio(n=0.42)        # e = n / (1 - n)
n = sp.porosity(e=0.72)          # n = e / (1 + e)

# Unit weights — all 4 kinds, any unit system
gamma_d   = sp.density(Gs=2.65, e=0.72, kind="dry",       unit="kN/m3")
gamma_sat = sp.density(Gs=2.65, e=0.72, kind="saturated", unit="kN/m3")
gamma_sub = sp.density(Gs=2.65, e=0.72, kind="submerged", unit="kN/m3")
gamma_b   = sp.density(Gs=2.65, e=0.72, S=0.8, kind="bulk", unit="pcf")

# Phase relations
S = sp.saturation(w=0.18, Gs=2.65, e=0.72)          # S = w·Gs / e
w = sp.water_content(S=0.80, Gs=2.65, e=0.72)        # w = S·e / Gs

# Relative density (Dr)
Dr = sp.relative_density(e=0.60, e_max=0.85, e_min=0.45, kind="void")
print(f"Dr = {Dr:.1%}")   # → Dr = 62.5%  (Medium Dense)

# Atterberg limits
PI = sp.atterberg(LL=48, PL=22, kind="PI")
all_idx = sp.atterberg(LL=48, PL=22, w=35, kind="all")
print(all_idx)  # → {"PI": 26, "LI": 0.50, "CI": 0.50}
```

---

## 📊 Grain Size Plot — Full Customization

`sp.grain_size_plot()` wraps Matplotlib and passes all standard **Line2D keyword arguments** directly to the renderer. Everything you can do with `ax.plot()`, you can do here.

```python
sp.grain_size_plot(
    data,
    # --- GeoEq Parameters ---
    smooth=True,           # bool: high-res PCHIP curve (1000 points)
    annotation=True,       # bool: show Sieve / Hydro legend
    D_para=True,           # bool: D10, D30, D60 red dotted projection lines
    Cu_para=True,          # bool: show Cu in parameter box
    Cc_para=True,          # bool: show Cc in parameter box
    param_pos="top right", # str or (x, y): parameter box position
    save_as="figure.pdf",  # str: any Matplotlib format (.png .svg .pdf .eps)
    ax=None,               # matplotlib Axes: embed into an existing figure
    # --- All Matplotlib kwargs ---
    color="#1e4d8f",
    linewidth=1.8,
    linestyle="-",
    alpha=0.9,
    marker="o",
    markersize=6,
    markerfacecolor="white",
    markeredgecolor="#1e4d8f",
    markeredgewidth=1.2,
)
```

### Multi-Dataset (Combined Sieve + Hydrometer)
Pass a **dictionary** of named datasets. GeoEq stitches them into one gapless smooth curve and applies automatic marker cycling per source — squares for the first dataset, stars for the second, etc.

```python
sp.grain_size_plot(
    {"Sieve Analysis": s_res, "Hydrometer Analysis": h_res},
    smooth=True,
    annotation=True,    # shows both names in the legend
)
```

### Embedding in a Multi-Panel Report Figure
The function accepts an existing `ax` and returns the `fig` object for further customization:

```python
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

sp.grain_size_plot(s_res, ax=ax1, smooth=True, color="navy")

ax2.plot(depth, settlement)
ax2.set_title("Settlement Curve")

plt.savefig("combined_report.pdf", bbox_inches="tight", dpi=300)
```

---

## 🗂️ Return Value Reference

All analysis functions return a **plain Python dictionary** so results work seamlessly with Pandas, NumPy, and Matplotlib.

### `sp.sieve_ana()` returns:
| Key | Description |
| :--- | :--- |
| `"diameter"` | Nominal sieve opening in mm (same as `"opening"`) |
| `"opening"` | Alias for `"diameter"` |
| `"mass_retained"` | Mass on each sieve (g) |
| `"percent_retained"` | % of total mass on each sieve |
| `"cumulative_retained"` | Cumulative % retained |
| `"percent_finer"` | **% passing — your primary Y-axis data** |
| `"total_mass"` | Total specimen mass used |

### `sp.hydro_ana()` returns:
| Key | Description |
| :--- | :--- |
| `"diameter"` | Equivalent particle diameter D (mm) via Stokes' Law |
| `"percent_finer"` | Corrected percent finer P (%) |
| `"reading"` | Corrected hydrometer readings (Rc) |
| `"depth"` | Effective depth L (cm) for each reading |
| `"time"` | Elapsed time (minutes) |

---

## 🏗️ Design Principles

- **Flat API** — everything is importable as `import geoeq as sp`; no long module chains.
- **Validated inputs** — every function checks that values are physically meaningful (e.g., porosity cannot exceed 1, saturation cannot exceed 1).
- **Traceable formulas** — docstrings cite the exact source (textbook section, ASTM standard).
- **No magic** — returns plain dicts and Matplotlib figures; nothing hidden.
- **Publication quality** — 300 DPI export, professional shaded zones, red-dotted Dx lines, and a clean parameter box out of the box.

---

## 🗺️ Roadmap

| Version | Scope |
|---------|-------|
| **v0.0.4** ✅ **current** | Particle size analysis — sieve, hydrometer, professional plotting |
| v0.0.3 | Atterberg limits, unified density function |
| v0.0.1–v0.0.2 | Soil phase relationships (16 functions) |
| v0.1.x | Effective stress, pore pressure, seepage |
| v0.2.x | Shallow foundations, bearing capacity, settlement |
| v0.3.x | Soil classification (USCS/AASHTO), SPT/CPT correlations |
| v0.4.x | Earth pressure, retaining walls, consolidation |
| v0.5.x | Pile capacity, slope stability, liquefaction |
| v1.0.0 | Stable API, full documentation, public open-source |

---

## Running the Tests

```bash
pip install -e ".[dev]"
pytest
```

---

## Contributing

The repository is **private** until approximately **v0.5**. Until then, feedback and feature requests are welcome.

---

## License

MIT License. See [LICENSE](LICENSE) for details.  
Copyright © 2026 Ripon Chandra Malo
