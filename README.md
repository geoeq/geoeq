<div align="center">

<img src="https://raw.githubusercontent.com/geoeq/geoeq/main/docs/assets/logo-animated.svg" alt="GeoEq" width="160" />

# GeoEq

**Geotechnical engineering, solved in Python.**

A clean, validated, MIT-licensed library for the complete onshore geotechnical workflow — laboratory characterisation, site investigation, engineering design, soil dynamics, and data exchange — under a single flat namespace.

<p>
  <a href="https://pypi.org/project/geoeq/"><img src="https://img.shields.io/pypi/v/geoeq.svg?color=1e4d8f&label=pypi" alt="PyPI version"></a>
  <a href="https://pypi.org/project/geoeq/"><img src="https://img.shields.io/pypi/pyversions/geoeq.svg?color=2563a8" alt="Python versions"></a>
  <a href="https://pypi.org/project/geoeq/"><img src="https://img.shields.io/pypi/dm/geoeq.svg?color=3b82c4&label=downloads" alt="PyPI downloads"></a>
  <a href="https://github.com/geoeq/geoeq/stargazers"><img src="https://img.shields.io/github/stars/geoeq/geoeq?style=flat&color=e8a825&label=stars" alt="GitHub stars"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/license-MIT-1e4d8f.svg" alt="MIT License"></a>
  <img src="https://img.shields.io/badge/tests-563_passing-27ae60" alt="563 tests passing">
  <img src="https://img.shields.io/badge/functions-170+-1e4d8f" alt="170+ functions">
</p>

<p>
  <a href="#installation"><b>Install</b></a> ·
  <a href="#quick-start"><b>Quick start</b></a> ·
  <a href="#capabilities"><b>Capabilities</b></a> ·
  <a href="https://geoeq.github.io"><b>Website</b></a> ·
  <a href="https://geoeq.github.io/user-guide.html"><b>Guide</b></a>
</p>

<br/>

<img src="https://raw.githubusercontent.com/geoeq/geoeq/main/docs/assets/hero.png" alt="GeoEq triptych — stress profile, bearing capacity vs footing width, liquefaction triggering chart" width="900" />

</div>

---

## The problem

Geotechnical analysis still runs on spreadsheets. One tab per soil layer, copy-pasted Boussinesq equations, magic-number unit conversions, `=IF()` ladders nobody else can read six months later. Reports are non-reproducible. Calculations get re-done from scratch every project. New engineers re-derive Terzaghi every time they touch a footing.

## The solution

`geoeq` is one flat, validated Python library covering the full onshore geotechnical workflow — lab tests through site investigation through design through dynamics — in 170+ functions under a single import. Plain dictionaries in, plain dictionaries out. Every formula cites its textbook source. Every function is validated against published values.

```python
import geoeq as ge

ge.density(Gs=2.65, e=0.72, kind="saturated", unit="kN/m3")           # 19.6 kN/m^3
ge.spt_friction_angle(N=15, sigma_v=80, method="hatanaka")            # 33.2 deg
ge.bearing_capacity(c=10, gamma=18, Df=1, B=2, phi=30, method="meyerhof")
ge.liquefaction_fos(CSR=0.20, CRR=0.18, Mw=7.0)                       # {'FS': 0.95, ...}
```

---

## Installation

```bash
pip install geoeq
```

Python 3.9+ on any platform. Dependencies: `numpy`, `matplotlib`, `scipy` — nothing else.

```python
import geoeq as ge
```

After import, every public function is a top-level attribute of `ge`. You never need to remember which submodule a function lives in.

---

## Quick start

A complete end-to-end screen of a multi-layer site, in fifteen lines.

```python
import geoeq as ge

# 1. Define a layered profile with a water table
p = ge.SoilProfile([
    (0, 2,  ge.Soil("Fill",       gamma=18)),
    (2, 8,  ge.Soil("Soft Clay",  gamma=17, gamma_sat=18.5,
                                  phi=0, c=25, Cc=0.27, e=0.92)),
    (8, 20, ge.Soil("Dense Sand", gamma=19, gamma_sat=20.5, phi=35)),
], water_table=2.0)

print(p.stress_at(10))
# {'sigma': 188.0, 'u': 78.48, 'sigma_eff': 109.52}

# 2. Bearing capacity of a 3 m square footing
bc = ge.bearing_capacity(c=25, gamma=18.5 - 9.81, Df=2, B=3, L=3,
                         phi=0, method="meyerhof")
print(f"q_u = {bc['q_u']:.0f} kPa")    # q_u = 192 kPa

# 3. Liquefaction triggering check (NCEER simplified procedure)
csr = ge.liquefaction_csr(amax=0.25, sigma_v=120, sigma_v_eff=70, z=6, Mw=7.0)
crr = ge.liquefaction_crr(N160cs=12, method="youd_2001")
fs  = ge.liquefaction_fos(csr["CSR"], crr["CRR"], Mw=7.0)
print(fs)
# {'FS': 0.58, 'liquefies': True, ...}
```

---

## Capabilities

A working summary of what ships in v0.1.2. The complete per-function reference lives in the [handbook](docs/GeoEq-Handbook.tex) and the [online guide](https://docs.geoeq.org).

### Soil properties and classification

| Function | Description |
| :-- | :-- |
| `ge.void_ratio()` · `ge.porosity()` · `ge.specific_gravity()` | Phase volumetrics |
| `ge.density()` | One function — dry / saturated / bulk / submerged, in kN/m³ or pcf |
| `ge.saturation()` · `ge.water_content()` | Phase relations |
| `ge.relative_density()` | Dr from void ratio or dry density limits |
| `ge.atterberg()` · `ge.activity()` · `ge.sensitivity()` · `ge.liquidity_index()` | Index correlations |
| `ge.classify_uscs()` · `ge.classify_aashto()` · `ge.plasticity_chart()` | ASTM D2487 · AASHTO M145 |

### Laboratory testing

| Suite | Functions |
| :-- | :-- |
| **Particle size** | `sieve_ana` · `hydro_ana` · `grain_size_plot` · `grain_d10/d30/d60` · `grain_Cu` · `grain_Cc` |
| **Shear strength** | `direct_shear` · `triaxial` · `unconfined` · `mohr_circle` |
| **Consolidation** | `oedometer` · `preconsolidation` · `compression_index` · `cv` (root + log time) |
| **Compaction** | `proctor` · `zav_line` · `saturation_line` · `relative_compaction` |
| **Permeability** | `constant_head` · `falling_head` |
| **Atterberg test** | `liquid_limit_test` · `flow_curve_plot` |
| **CBR** | `cbr_test` · `cbr_plot` |

### Site investigation

| Suite | Functions |
| :-- | :-- |
| **SPT** | `spt_n60` · `spt_n160` (3 methods) · `spt_n160cs` · `spt_friction_angle` (3) · `spt_su` (2) · `spt_dr` (3) · `spt_modulus` (6 soils) |
| **CPT** | `cpt_normalize` · `cpt_ic` · `cpt_sbt` (Robertson) · `cpt_friction_angle` · `cpt_su` · `cpt_dr` · `cpt_modulus` · `cpt_sbt_plot` |
| **Field vane** | `vane_su` · `vane_correction` (Bjerrum) · `vane_remolded` |
| **Pressuremeter** | `pmt_parameters` · `pmt_modulus` · `pmt_su` · `pmt_bearing` · `pmt_settlement` · `pmt_ko` |
| **Plate load** | `plt_bearing` (3 criteria) · `plt_subgrade_modulus` · `plt_settlement_correction` · `plt_elastic_modulus` |
| **Pile load tests** | `davisson` · `chin` · `de_beer` · `hansen_80` · `fhwa_5_percent` · `case_method` · `hiley` · `danish_formula` · `enr` |
| **Field permeability** | `slug_test` (Hvorslev) · `pumping_test_confined` · `pumping_test_unconfined` (Thiem) · `lefranc_test` |
| **Field CBR / DCP** | `dcp_cbr` (Webster · TRL · Kleyn) · `field_cbr_test` |

### Engineering design

| Suite | Functions |
| :-- | :-- |
| **Effective stress** | `total_stress` · `pore_pressure` · `effective_stress` · `capillary_rise` · `stress_plot` |
| **Seepage** | `darcy_flow` · `hydraulic_gradient` · `critical_gradient` · `equivalent_k` · `flow_net` |
| **Stress distribution** | `boussinesq_point/line/strip/circular/rect` · `westergaard_point` · `newmark_influence` · `stress_2to1` · `stress_isobar_plot` |
| **Bearing capacity** | Terzaghi · Meyerhof · Hansen · Vesic — shape, depth, inclination factors · `bearing_capacity_plot` |
| **Settlement** | `settlement_immediate` (Janbu) · `settlement_primary` (NC/OC/crossing) · `settlement_secondary` · `settlement_schmertmann` (1978) · `time_factor` ↔ `consolidation_degree` · `settlement_time_plot` |
| **Earth pressure** | `K0` (Jaky · Mayne-Kulhawy) · `Ka` / `Kp` (Rankine · Coulomb) · `earth_pressure` · `tension_crack_depth` · `earth_pressure_plot` |
| **Retaining walls** | `wall_overturning` · `wall_sliding` · `wall_bearing` (kern check) · `sheet_pile` (Blum) |
| **Pile design** | `pile_end_bearing` (Meyerhof · Vesic · Skempton) · `pile_skin_friction` (alpha · beta · lambda) · `pile_capacity` · `pile_group_efficiency` · `pile_settlement` (Vesic) |
| **Slope stability** | `infinite_slope` (with seepage) · `culmann` · `taylor_stability` · `bishop` (iterative) · `taylor_chart_plot` |

### Soil dynamics

| Function | Description |
| :-- | :-- |
| `gmax()` · `gmax_hardin()` | G<sub>max</sub> from V<sub>s</sub> or Hardin & Black (1968) |
| `modulus_reduction()` | G/G<sub>max</sub> — Darendeli (2001), Vucetic-Dobry (1991) |
| `damping_ratio()` | Equivalent-linear damping (Darendeli + Hardin-Drnevich) |
| `depth_reduction()` | r<sub>d</sub> from Idriss (1999), Liao-Whitman (1986), Cetin et al. (2004) |
| `liquefaction_csr()` | Seed-Idriss (1971) cyclic stress ratio |
| `liquefaction_crr()` | CRR — Youd et al. (2001), Idriss-Boulanger (2008) SPT/CPT, Andrus-Stokoe (2000) V<sub>s</sub> |
| `magnitude_scaling_factor()` | Idriss · NCEER · Boulanger-Idriss (2014) |
| `liquefaction_fos()` | FS<sub>L</sub> = CRR · MSF · K<sub>σ</sub> · K<sub>α</sub> / CSR |
| `gmax_curves_plot()` · `liquefaction_chart()` | Publication-quality charts |

### Layered ground model

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
| `read_csv()` | Geotech CSV with auto-header and units-row detection |
| `read_ags()` | AGS4 format (UK / Australia / NZ ground-investigation standard) |
| `read_gef()` | GEF-CPT format (Dutch standard) |
| `CPT()` · `.from_gef()` · `.from_ags()` | CPT container class |

---

## Design principles

- **Flat API.** `import geoeq as ge`, then call functions. No deep import chains, no class hierarchies.
- **Validated inputs.** Every function checks physical meaningfulness (porosity ≤ 1, saturation ≤ 1, etc.) with engineer-readable errors.
- **Traceable formulas.** Every docstring cites the textbook section or original paper. No mystery constants.
- **No magic.** Plain `dict` returns, Matplotlib figures. Inspect, slice, feed into Pandas, do whatever you would do with NumPy output.
- **Test-backed.** 563 tests across 170+ functions, all passing. Textbook values are spot-checked: Fadum I<sub>5</sub>(1,1) = 0.1752; Meyerhof N-factors at φ = 30°; Boussinesq circular at z = R; Rankine K<sub>a</sub> / K<sub>p</sub>; Hardin & Drnevich plasticity exponent k.
- **Publication quality.** 300 DPI figures, semi-log axes, shaded ASTM particle zones, red-dashed D<sub>x</sub> projection lines out of the box.
- **MIT forever.** Commercial consulting, in-house tools, research, teaching. No exceptions, no dual licensing.

---

## Running the tests

```bash
git clone https://github.com/geoeq/geoeq.git
cd geoeq
pip install -e ".[dev]"
pytest                       # -> 563 passed
```

---

## Contributing

The repository is open for issue reports and feedback. Pull requests are welcome once the project reaches v1.0; until then, please open an issue first to discuss any proposed change.

- Found a formula bug or a unit mistake? [Open an issue](https://github.com/geoeq/geoeq/issues) with a textbook citation.
- Want a function prioritised? Open an issue tagged `enhancement` describing the use case.
- Using GeoEq in coursework or research? Open an issue tagged `usage` — we would like to hear about real-world uses.

---

## Citation

If you use GeoEq in academic work, please cite:

```bibtex
@software{geoeq2026,
  author       = {Malo, Ripon Chandra},
  title        = {GeoEq: A Python Library for Geotechnical Engineering},
  year         = {2026},
  version      = {0.1.2},
  url          = {https://github.com/geoeq/geoeq},
  license      = {MIT}
}
```

---

## License

MIT — free for personal, commercial, consulting, and enterprise use forever. See [LICENSE](LICENSE).

Copyright © 2026 Ripon Chandra Malo · University of Utah

---

<div align="center">
<sub><i>GeoEq — geotechnical engineering, solved in Python.</i></sub>
<br/>
<sub>If GeoEq saved you a spreadsheet, please <a href="https://github.com/geoeq/geoeq">star the repository</a> — it is the cheapest way to help.</sub>
</div>
