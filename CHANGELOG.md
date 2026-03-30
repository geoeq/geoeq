# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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
