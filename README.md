# GeoEq

A Python library for geotechnical engineering calculations.

The installable package name is **`geoeq`** (like **NumPy** / `numpy`). **GeoEq**
is the project name used in documentation.

**GeoEq** brings common geotechnical formulas into one consistent, validated,
and well-documented Python package. It is built for engineers, students, and
researchers who want to replace repetitive spreadsheet work with clean,
reusable code.

---

## Status

**v0.0.3** is the current release. It covers **soil property and phase-relationship formulas** --
the building blocks for everything else in geotechnical engineering.

Future releases will add effective stress, bearing capacity, settlement,
consolidation, earth pressure, slope stability, pile capacity, liquefaction
assessment, and more.

---

## Installation

```bash
pip install geoeq
```

Requires Python 3.9 or later. No external dependencies.

---

## Quick start

```python
from geoeq.soil import properties as sp

# Unified density function (with unit support and mass/volume options)
gamma_d = sp.density(Gs=2.65, e=0.72, kind="dry", unit="kN/m3")
print(f"Dry unit weight: {gamma_d:.2f} kN/m3")

# Choose output unit directly (kN/m3, kg/m3, pcf, g/cm3)
rho_bulk = sp.density(Gs=2.65, e=0.72, S=0.8, kind="bulk", unit="kg/m3")
print(f"Bulk density: {rho_bulk:.1f} kg/m3")

# Or calculate directly from mass and volume
sample_density = sp.density(mass=5.50, volume=0.003, unit="kN/m3")

# Phase relations
n = sp.porosity_from_void_ratio(0.72)
S = sp.degree_of_saturation(w=0.18, Gs=2.65, e=0.72)
w = sp.water_content(S=0.80, Gs=2.65, e=0.72)
e = sp.void_ratio_from_porosity(0.42)

# Unified Atterberg limits
PI = sp.atterberg(LL=48, PL=22, kind="PI")
all_indices = sp.atterberg(LL=48, PL=22, w=35, kind="all")
```

---

## What is included in v0.0.3

### Soil property and phase-relationship formulas

All functions validate inputs, document units in docstrings, and reference
their source methods.

**Void ratio and porosity**
- `void_ratio_from_porosity(n)` -- convert porosity to void ratio
- `porosity_from_void_ratio(e)` -- convert void ratio to porosity

**Degree of saturation and water content**
- `degree_of_saturation(w, Gs, e)` -- compute S from water content, Gs, and e
- `water_content(S, Gs, e)` -- compute w from saturation, Gs, and e

**Unit weights**
- `dry_unit_weight(Gs, e)` -- dry unit weight (kN/m3)
- `saturated_unit_weight(Gs, e)` -- saturated unit weight (kN/m3)
- `bulk_unit_weight(Gs, e, S)` -- moist/total unit weight (kN/m3)
- `submerged_unit_weight(Gs, e)` -- buoyant unit weight (kN/m3)

**Density**
- `dry_density(Gs, e)` -- dry density (kg/m3)
- `bulk_density(Gs, e, S)` -- bulk density (kg/m3)

**Unified API functions**
- `density(Gs, e, S, mass, volume, kind, unit)` -- all density/unit weight formulas combined, supports custom units and fallback combinations.

**Void ratio from other quantities**
- `void_ratio_from_water_content(w, Gs, S)` -- compute e from w, Gs, S
- `void_ratio_from_dry_unit_weight(gamma_d, Gs)` -- compute e from gamma_d

**Relative density**
- `relative_density(e, e_max, e_min)` -- density index for granular soils

**Atterberg limits and consistency**
- `plasticity_index(LL, PL)` -- PI = LL - PL
- `liquidity_index(w, PL, PI)` -- LI
- `consistency_index(w, LL, PI)` -- CI
- `atterberg(LL, PL, w, kind)` -- unified limits function

### Core utilities

- **Validation** -- input checking (positive, non-negative, range, fraction, angle)
- **Constants** -- gravity, unit weight of water, atmospheric pressure
- **Types** -- `Soil` dataclass for passing soil properties between functions

---

## Design principles

- **Validated inputs.** Every function checks that values are physically
  meaningful. Negative void ratios, impossible saturation, and bad angles
  raise clear errors.
- **Documented units.** Every docstring states what units are expected and
  returned.
- **Traceable methods.** Docstrings reference the source formula or textbook.
- **No hidden dependencies.** GeoEq uses only the Python standard library.
- **Tested.** Every formula is tested against hand-calculated values and
  textbook examples.

---

## Roadmap

| Version | Scope |
|---------|-------|
| v0.0.1  | Soil properties and phase relations |
| v0.1.x  | Effective stress, pore pressure, seepage |
| v0.2.x  | Shallow foundations, bearing capacity, settlement |
| v0.3.x  | Soil classification (USCS), SPT/CPT correlations |
| v0.4.x  | Earth pressure, consolidation, slope stability |
| v0.5.x  | Pile capacity, liquefaction, advanced methods |
| v1.0.0  | Stable API, full documentation, examples |

---

## Running the tests

From a clone of this repository:

```bash
pip install -e ".[dev]"
pytest
```

If you only install pytest, tests still run from the project root (see `pythonpath` in `pyproject.toml`):

```bash
pip install pytest
pytest
```

The wheel on PyPI does not include the test files; run tests from a source checkout.

---

## Contributing

The repository is **private** until around **v0.5**; it will then be opened for
public contributions. Until then, feedback and feature requests are appreciated.

Every contribution should include the code, a test, a documentation update,
and a reference for the method.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

Copyright (c) 2026 Ripon Chandra Malo
