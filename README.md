<!-- GeoEq — Open-source Python library for geotechnical engineering calculations, soil mechanics, foundation design, SPT, CPT, bearing capacity, settlement analysis, slope stability, and liquefaction assessment. Apache 2.0 licensed. -->

<div align="center">

<a href="https://geoeq.org">
<img src="https://raw.githubusercontent.com/geoeq/geoeq/main/docs/assets/hero-animated.svg" alt="GeoEq — Geotechnical Engineering Python Library: soil mechanics, bearing capacity, SPT, CPT, settlement, liquefaction" width="900" />
</a>

<br/>

<p>
  <a href="https://pypi.org/project/geoeq/"><img src="https://img.shields.io/pypi/v/geoeq.svg?style=flat&color=1e4d8f&label=PyPI" alt="PyPI version"></a>
  <a href="https://pypi.org/project/geoeq/"><img src="https://img.shields.io/pypi/dm/geoeq.svg?style=flat&color=2563a8&label=Downloads/month" alt="Monthly PyPI downloads"></a>
  <a href="https://pepy.tech/project/geoeq"><img src="https://static.pepy.tech/badge/geoeq?style=flat" alt="Total PyPI downloads"></a>
  <a href="https://pypi.org/project/geoeq/"><img src="https://img.shields.io/pypi/pyversions/geoeq.svg?style=flat&color=4a90d9" alt="Python 3.9 | 3.10 | 3.11 | 3.12 | 3.13"></a>
  <a href="https://www.apache.org/licenses/LICENSE-2.0"><img src="https://img.shields.io/badge/License-Apache_2.0-1e4d8f.svg?style=flat" alt="Apache 2.0 License"></a>
</p>
<p>
  <a href="https://geoeq.org"><img src="https://img.shields.io/badge/Docs-geoeq.org-27ae60?style=flat&logo=readthedocs&logoColor=white" alt="Documentation — geoeq.org"></a>
  <a href="https://github.com/geoeq/geoeq/stargazers"><img src="https://img.shields.io/github/stars/geoeq/geoeq?style=flat&color=e8a825&label=Stars" alt="GitHub stars"></a>
  <img src="https://img.shields.io/badge/Tests-563_passing-27ae60?style=flat" alt="563 tests passing">
  <img src="https://img.shields.io/badge/Functions-170+-1e4d8f?style=flat" alt="170+ geotechnical functions">
</p>

<h3>Geotechnical engineering, solved in Python.</h3>

<p><i>A clean, validated, Apache 2.0-licensed Python library for the complete onshore geotechnical workflow —<br/>
laboratory testing, site investigation, foundation design, soil dynamics, and data exchange —<br/>
under a single flat namespace.</i></p>

<p>
  <a href="#installation"><b>Install</b></a> ·
  <a href="#quick-start"><b>Quick Start</b></a> ·
  <a href="#capabilities"><b>Capabilities</b></a> ·
  <a href="https://geoeq.org"><b>Documentation</b></a> ·
  <a href="https://geoeq.org/user-guide.html"><b>User Guide</b></a> ·
  <a href="https://geoeq.org/tutorials.html"><b>Tutorials</b></a>
</p>

</div>

---

## Why GeoEq?

Geotechnical analysis still runs on spreadsheets — one tab per soil layer, copy-pasted Boussinesq equations, magic-number unit conversions, and `=IF()` ladders nobody else can read six months later. Reports are non-reproducible. Calculations get re-derived from scratch every project.

**GeoEq replaces that entire workflow with one Python import.**

```python
import geoeq as ge

# Soil classification, bearing capacity, liquefaction — all in one namespace
ge.density(Gs=2.65, e=0.72, kind="saturated", unit="kN/m3")
ge.bearing_capacity(c=10, gamma=18, Df=1, B=2, phi=30, method="meyerhof")
ge.liquefaction_fos(CSR=0.20, CRR=0.18, Mw=7.0)
```

170+ validated functions. Plain dicts in, plain dicts out. Every formula cites its textbook source. Every function is tested against published values. No class hierarchies, no deep imports, no magic. Apache 2.0 licensed.

---

## Installation

```bash
pip install geoeq
```

Python 3.9+ on any platform. Dependencies: `numpy`, `matplotlib`, `scipy` — nothing else.

```python
import geoeq as ge
# Every function is a top-level attribute — no submodule imports needed
```

---

## Quick Start

A complete site screening in fifteen lines — soil profiling, bearing capacity, and liquefaction check:

```python
import geoeq as ge

# 1. Define a layered soil profile with water table
p = ge.SoilProfile([
    (0, 2,  ge.Soil("Fill",       gamma=18)),
    (2, 8,  ge.Soil("Soft Clay",  gamma=17, gamma_sat=18.5,
                                  phi=0, c=25, Cc=0.27, e=0.92)),
    (8, 20, ge.Soil("Dense Sand", gamma=19, gamma_sat=20.5, phi=35)),
], water_table=2.0)

print(p.stress_at(10))
# {'sigma': 188.0, 'u': 78.48, 'sigma_eff': 109.52}

# 2. Bearing capacity — Meyerhof method for a 3 m square footing
bc = ge.bearing_capacity(c=25, gamma=18.5 - 9.81, Df=2, B=3, L=3,
                         phi=0, method="meyerhof")
print(f"q_u = {bc['q_u']:.0f} kPa")    # q_u = 192 kPa

# 3. Liquefaction triggering (NCEER simplified procedure)
csr = ge.liquefaction_csr(amax=0.25, sigma_v=120, sigma_v_eff=70, z=6, Mw=7.0)
crr = ge.liquefaction_crr(N160cs=12, method="youd_2001")
fs  = ge.liquefaction_fos(csr["CSR"], crr["CRR"], Mw=7.0)
print(fs)   # {'FS': 0.58, 'liquefies': True, ...}
```

---

## Capabilities

Full per-function reference at [geoeq.org](https://geoeq.org). Below is what ships in v0.1.3.

### Soil Properties & Classification

| Function | Description |
| :-- | :-- |
| `ge.void_ratio()` · `ge.porosity()` · `ge.specific_gravity()` | Phase volumetrics |
| `ge.density()` | Dry / saturated / bulk / submerged — kN/m³ or pcf |
| `ge.saturation()` · `ge.water_content()` | Phase relations |
| `ge.relative_density()` | D<sub>r</sub> from void ratio or dry density limits |
| `ge.atterberg()` · `ge.activity()` · `ge.sensitivity()` · `ge.liquidity_index()` | Atterberg limits & index properties |
| `ge.classify_uscs()` · `ge.classify_aashto()` · `ge.plasticity_chart()` | ASTM D2487 USCS · AASHTO M145 classification |

### Laboratory Testing

| Suite | Functions |
| :-- | :-- |
| **Grain size analysis** | `sieve_ana` · `hydro_ana` · `grain_size_plot` · `grain_d10/d30/d60` · `grain_Cu` · `grain_Cc` |
| **Shear strength** | `direct_shear` · `triaxial` (UU/CU/CD) · `unconfined` · `mohr_circle` · `direct_shear_plot` |
| **Consolidation** | `oedometer` · `preconsolidation` (Casagrande) · `compression_index` (6 methods) · `cv` (log-time & root-time) |
| **Compaction** | `proctor` · `zav_line` · `saturation_line` · `relative_compaction` · `proctor_plot` |
| **Permeability** | `constant_head` · `falling_head` |
| **Atterberg testing** | `liquid_limit_test` · `flow_curve_plot` |
| **CBR** | `cbr_test` · `cbr_plot` |

### Site Investigation

| Suite | Functions |
| :-- | :-- |
| **Standard Penetration Test (SPT)** | `spt_n60` · `spt_n160` (3 methods) · `spt_n160cs` · `spt_friction_angle` (3) · `spt_su` (2) · `spt_dr` (3) · `spt_modulus` (6 soils) |
| **Cone Penetration Test (CPT)** | `cpt_normalize` · `cpt_ic` · `cpt_sbt` (Robertson) · `cpt_friction_angle` · `cpt_su` · `cpt_dr` · `cpt_modulus` · `cpt_sbt_plot` |
| **Field vane shear** | `vane_su` · `vane_correction` (Bjerrum) · `vane_remolded` |
| **Pressuremeter (PMT)** | `pmt_parameters` · `pmt_modulus` · `pmt_su` · `pmt_bearing` · `pmt_settlement` · `pmt_ko` |
| **Plate load test (PLT)** | `plt_bearing` (3 criteria) · `plt_subgrade_modulus` · `plt_settlement_correction` · `plt_elastic_modulus` |
| **Pile load test** | `davisson` · `chin` · `de_beer` · `hansen_80` · `fhwa_5_percent` · `case_method` · `hiley` · `danish_formula` · `enr` |
| **Field permeability** | `slug_test` (Hvorslev) · `pumping_test_confined` · `pumping_test_unconfined` (Thiem) · `lefranc_test` |
| **Field CBR / DCP** | `dcp_cbr` (Webster · TRL · Kleyn) · `field_cbr_test` |

### Foundation & Geotechnical Design

| Suite | Functions |
| :-- | :-- |
| **Effective stress analysis** | `total_stress` · `pore_pressure` · `effective_stress` · `capillary_rise` · `stress_plot` |
| **Seepage & flow** | `darcy_flow` · `hydraulic_gradient` · `critical_gradient` · `equivalent_k` · `flow_net` |
| **Stress distribution** | `boussinesq_point/line/strip/circular/rect` · `westergaard_point` · `newmark_influence` · `stress_2to1` · `stress_isobar_plot` |
| **Bearing capacity** | Terzaghi · Meyerhof · Hansen · Vesic — shape, depth, inclination factors · `bearing_capacity_plot` |
| **Settlement analysis** | `settlement_immediate` (Janbu) · `settlement_primary` (NC/OC) · `settlement_secondary` · `settlement_schmertmann` · `time_factor` · `consolidation_degree` |
| **Earth pressure** | `K0` (Jaky · Mayne-Kulhawy) · `Ka` / `Kp` (Rankine · Coulomb) · `earth_pressure` · `tension_crack_depth` |
| **Retaining walls** | `wall_overturning` · `wall_sliding` · `wall_bearing` (kern check) · `sheet_pile` (Blum) |
| **Pile foundation design** | `pile_end_bearing` (Meyerhof · Vesic · Skempton) · `pile_skin_friction` (α · β · λ) · `pile_capacity` · `pile_group_efficiency` · `pile_settlement` |
| **Slope stability** | `infinite_slope` · `culmann` · `taylor_stability` · `bishop` (iterative) · `taylor_chart_plot` |

### Soil Dynamics & Earthquake Engineering

| Function | Description |
| :-- | :-- |
| `gmax()` · `gmax_hardin()` | Small-strain shear modulus G<sub>max</sub> |
| `modulus_reduction()` | G/G<sub>max</sub> — Darendeli (2001), Vucetic-Dobry (1991) |
| `damping_ratio()` | Equivalent-linear damping (Darendeli, Hardin-Drnevich) |
| `depth_reduction()` | r<sub>d</sub> — Idriss (1999), Liao-Whitman (1986), Cetin (2004) |
| `liquefaction_csr()` | Cyclic stress ratio (Seed-Idriss 1971) |
| `liquefaction_crr()` | CRR — Youd (2001), Idriss-Boulanger (2008), Andrus-Stokoe (2000) |
| `magnitude_scaling_factor()` | MSF — Idriss · NCEER · Boulanger-Idriss (2014) |
| `liquefaction_fos()` | Factor of safety against liquefaction |

### Layered Ground Model

```python
clay = ge.Soil("Soft Clay", gamma=17, gamma_sat=18.5, phi=0, c=25, e=0.9)

p = ge.SoilProfile([
    (0, 2,  ge.Soil("Fill",       gamma=18)),
    (2, 8,  clay),
    (8, 20, ge.Soil("Dense Sand", gamma=19, gamma_sat=20.5, phi=35)),
], water_table=1.5)

p.effective_stress(10)   # σ' at z = 10 m
p.plot()                 # publication-quality stress profile
```

### Data I/O

| Function | Description |
| :-- | :-- |
| `read_csv()` | Geotech CSV with auto-header detection |
| `read_ags()` | AGS4 format (UK / AU / NZ ground investigation) |
| `read_gef()` | GEF-CPT format (Dutch standard) |
| `CPT()` · `.from_gef()` · `.from_ags()` | CPT container class |

---

## Design Principles

| Principle | What it means |
| :-- | :-- |
| **Flat API** | `import geoeq as ge` — every function at top level. No deep import chains. |
| **Validated inputs** | Physical meaningfulness checks (porosity ≤ 1, saturation ≤ 1) with clear error messages. |
| **Traceable formulas** | Every docstring cites the textbook or paper. No mystery constants. |
| **Plain returns** | `dict` results, Matplotlib figures. Inspect, slice, feed into Pandas. |
| **Test-backed** | 563 tests across 170+ functions. Spot-checked against Das, Bowles, NCEER. |
| **Publication quality** | 300 DPI figures, semi-log axes, ASTM particle zones, out of the box. |
| **Apache 2.0 forever** | Consulting, research, teaching, enterprise. Patent protection included. |

---

## Tutorials

27 Jupyter notebooks covering all 163+ functions with synthetic data, step-by-step explanations, and real-world geotechnical problems:

| Chapter | Topic | Functions |
| :-- | :-- | :-- |
| 01 | Soil Properties & Classification | `density`, `atterberg`, `classify_uscs`, `classify_aashto` |
| 02 | Phase Relations | `void_ratio`, `porosity`, `saturation`, `water_content` |
| 03 | Grain Size Analysis | `sieve_ana`, `hydro_ana`, `grain_size_plot` |
| 04 | Shear Strength | `direct_shear`, `triaxial`, `mohr_circle` |
| 05–06 | Consolidation & Compaction | `oedometer`, `proctor`, `cv` |
| 07–10 | Permeability, Atterberg, CBR | `constant_head`, `falling_head`, `cbr_test` |
| 11–17 | Site Investigation (SPT, CPT, PMT, PLT) | Full SPT/CPT/PMT/PLT/vane suites |
| 18–20 | Stress Analysis & Seepage | `boussinesq_*`, `darcy_flow`, `flow_net` |
| 21–27 | Design (Bearing, Settlement, Piles, Slopes) | `bearing_capacity`, `settlement_*`, `bishop` |

Browse notebooks: [`test_and_tutorial/`](test_and_tutorial/)

---

## Running the Tests

```bash
git clone https://github.com/geoeq/geoeq.git
cd geoeq
pip install -e ".[dev]"
pytest                       # → 563 passed
```

---

## Contributing

Pull requests are welcome once the project reaches v1.0; until then, please open an issue first.

- **Formula bug or unit error?** [Open an issue](https://github.com/geoeq/geoeq/issues) with a textbook citation.
- **Feature request?** Open an issue tagged `enhancement` with the use case.
- **Using GeoEq in research or teaching?** Open an issue tagged `usage` — we'd love to hear about it.

---

## Citation

```bibtex
@software{geoeq2026,
  author       = {Malo, Ripon Chandra},
  title        = {GeoEq: A Python Library for Geotechnical Engineering},
  year         = {2026},
  version      = {0.1.3},
  url          = {https://github.com/geoeq/geoeq},
  license      = {Apache-2.0}
}
```

---

## License

Apache 2.0 — free for personal, commercial, consulting, and enterprise use forever, with explicit patent protection. See [LICENSE](LICENSE).

Copyright © 2026 Ripon Chandra Malo · University of Utah

---

## Built With

This package was built using [**projectmem**](https://github.com/riponcm/projectmem) — a local-first memory and judgment layer for AI coding agents that tracks decisions, prevents repeated mistakes, and maintains project context across sessions.

---

<div align="center">

**[geoeq.org](https://geoeq.org)** · **[Documentation](https://geoeq.org/user-guide.html)** · **[GitHub](https://github.com/geoeq/geoeq)** · **[PyPI](https://pypi.org/project/geoeq/)**

<sub><i>GeoEq — geotechnical engineering, solved in Python.</i></sub>
<br/>
<sub>Soil mechanics · Foundation engineering · SPT · CPT · Bearing capacity · Settlement · Liquefaction · Slope stability</sub>
<br/><br/>
<sub>If GeoEq saved you a spreadsheet, please <a href="https://github.com/geoeq/geoeq">⭐ star the repository</a> — it's the cheapest way to help.</sub>

</div>
