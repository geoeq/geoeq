# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

- Effective stress and pore pressure calculations
- Seepage and hydraulic gradient functions
- USCS / AASHTO soil classification strings

---

## [0.0.4] - 2026-04-04

### Added — Particle Size Distribution Suite (Chapter 2)

#### Analysis Functions
- **`sp.sieve_ana(opening, mass_retained, standard, total_mass)`**
  - Full sieve analysis with automatic designation-to-mm mapping
  - Supports: ASTM E11 (`#4`, `#10`, `#200`, `3/4"`, etc.), BS 410, IS 460
  - Returns: `diameter`, `opening`, `mass_retained`, `percent_retained`, `cumulative_retained`, `percent_finer`, `total_mass`

- **`sp.hydro_ana(reading, time, T, Gs, Ws, Fm, Fz)`**
  - Full ASTM D7928 hydrometer sedimentation analysis
  - Implements Stokes' Law for equivalent particle diameter D
  - Applies meniscus correction (Fm), zero correction (Fz), and temperature corrections
  - Calibrated for the standard **152H hydrometer** model
  - Returns: `diameter`, `percent_finer`, `reading`, `depth`, `time`

#### Visualization — `sp.grain_size_plot()`
- **Smooth curves**: High-resolution PCHIP interpolation (1000 points) for a "perfectly smooth, single-line" academic curve
- **Separated line + markers**: Smooth line is plotted independently from original data point markers, preventing the chunky "step" appearance of raw data
- **Multi-dataset mode**: Accept `{"Sieve": s_res, "Hydro": h_res}` — datasets are stitched into one gapless curve; automatic distinct marker cycling per source (squares, stars, triangles…)
- **Annotation** (`annotation=True`): Names each data source in the legend
- **D-parameter projection lines** (`D_para=True`): Red dotted horizontal + vertical lines at D₁₀, D₃₀, D₆₀ with intersection markers
- **Professional shading**: Gravel (grey), Sand (beige), Fines (light blue) zone shading with boundary lines
- **Customizable parameter box** (`param_pos`): accepts `"top right"`, `"top left"`, `"bottom right"`, `"bottom left"`, or a precise `(x, y)` tuple in axes-fraction coordinates
- **Full Matplotlib pass-through**: All `Line2D` kwargs work — `color`, `linewidth`, `linestyle`, `alpha`, `marker`, `markersize`, `markerfacecolor`, `markeredgecolor`, `markeredgewidth`
- **Multi-panel embedding**: Accepts an existing `ax` argument and returns `fig` for downstream customization
- **Publication-quality export** (`save_as`): 300 DPI, tight bounding box; supports `.png`, `.svg`, `.pdf`, `.eps`

#### Parameter Extraction
- **`sp.grain_d10(data)`** — D₁₀ effective size via log-linear interpolation
- **`sp.grain_d30(data)`** — D₃₀ via log-linear interpolation
- **`sp.grain_d60(data)`** — D₆₀ controlling size via log-linear interpolation
- **`sp.grain_Cu(data)`** — Uniformity Coefficient Cᵤ = D₆₀ / D₁₀
- **`sp.grain_Cc(data)`** — Coefficient of Curvature Cᶜ = D₃₀² / (D₁₀ × D₆₀)

#### Dependencies Added
- `numpy` — array operations for grain size calculations
- `matplotlib` — professional plot rendering
- `scipy` — PCHIP monotonic interpolation (`PchipInterpolator`)

#### Flat API Exports
All new functions are available directly under `import geoeq as sp` (no submodule path required).

### Changed
- `sieve_ana()` return dict now includes a `"diameter"` key (alias for `"opening"`) for unified compatibility with the plotting function
- Duplicate `return fig` / `save_as` blocks in `viz/grain_size.py` cleaned up
- `pyproject.toml` version bumped to `0.0.4`

### Documentation
- **`README.md`** completely rewritten with full API table, return value reference, customization guide, and multi-panel embedding example
- **`user-guide.html`** restructured to match Braja Das chapter order (Ch.2 before Ch.3) with deeply detailed sections for every function argument and output key
- **`index.html`** homepage updated with new feature card, refreshed Quick Start code block, and roadmap card updated to v0.0.4

---

## [0.0.3] - 2026-03-30

### Added

- `sp.density()`: Unified density calculation function supporting variable state modes (`kind`), multiple output formats (`unit`), and basic (`mass`, `volume`) calculations.
- `sp.atterberg()`: Unified function for plasticity index (`PI`), liquidity index (`LI`), and consistency index (`CI`).

---

## [0.0.2] - 2026-03-27

### Changed

- README and package docstrings use the **GeoEq** project name; the PyPI and
  import name remains `geoeq` (similar to NumPy / `numpy`).

---

## [0.0.1] - 2026-03-27

### Added

- `geoeq.soil.properties` module with 16 soil phase-relationship functions:
  - void ratio / porosity conversions
  - degree of saturation
  - water content
  - dry, saturated, bulk, and submerged unit weights
  - dry and bulk density
  - void ratio from water content and from dry unit weight
  - relative density
  - plasticity index, liquidity index, consistency index
- `geoeq.core.validation` module with input validation helpers
- `geoeq.core.constants` module with physical constants (gravity, gamma_w)
- `geoeq.core.types` module with `Soil` dataclass
- Test suite with 43 tests covering all public functions
- MIT license
- README with installation, usage, API reference, and roadmap
