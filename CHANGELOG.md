# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

- Improved Schmertmann settlement chart-based Iz reading
- Bishop's rigorous (non-simplified) method
- Janbu's simplified slope-stability method
- Cetin et al. (2018) updated CRR curves
- Newmark sliding-block displacement

---

## [0.1.3] - 2026-07-18

### Changed
- **License changed from MIT to Apache 2.0** — same permissive freedoms,
  now with an explicit patent grant and contributor patent protection.
  Added `NOTICE` file per Apache 2.0 convention.
- Maintainer contact email updated to `support@matily.org`.
- Security policy now directs reports to GitHub Issues instead of email.
- README redesigned: animated SVG hero banner, download badges,
  SEO-optimised structure, and a tutorials index.

### Note
No functional code changes. All 563 tests still pass.

---

## [0.1.2] - 2026-05-20

### Fixed
- Added explicit `src/geoeq/viz/__init__.py` so the `geoeq.viz`
  submodule is recognised by static documentation generators
  (mkdocstrings / griffe). Without this file, `geoeq.viz` was a
  namespace package — fine at runtime but invisible to static
  parsers, which caused the documentation site's API reference
  generation to fail on `ge.grain_size_plot` and related symbols.

### Note
No functional code changes from v0.1.1. This is a packaging
hygiene release to support the public documentation build at
`docs.geoeq.org`. All 563 tests still pass.

---

## [0.1.1] - 2026-05-20

### Removed
- Internal planning documents removed from the public repository
  (`Implementation.md`, `PUBLISHING.md`). These remain available
  locally for maintainers.
- README "Contributing" section no longer references the internal
  roadmap document.

### Changed
- `pyproject.toml` version bumped to `0.1.1`.
- `MANIFEST.in` no longer ships `Implementation.md` in the sdist.
- `.gitignore` extended to keep maintainer-only documents out of
  the public repository.

### Note
This release contains no functional code changes from v0.1.0. It is
a documentation-hygiene release to remove internal planning
material from the public distribution.

---

## [0.1.0] - 2026-05-20

### Added — Complete Engineering Design Module (`ge.design`)

The full Phase 1–Phase 5 design scope shipped in one release. **~75 new
functions across 3 new modules.** Tests increased from 404 to 563.

#### New Module — `ge.profile` (Layered Soil Profiles)
- **`ge.Soil(...)`** — dataclass for one soil layer (gamma, phi, c, e, k, Es, mu, Cc, Cr, OCR, Su)
- **`ge.SoilProfile(layers, water_table)`** — layered profile with stress computations
  - `total_stress(z)`, `pore_pressure(z)`, `effective_stress(z)`, `stress_at(z)`
  - `layer_at(z)`, `add_layer()`, `to_dataframe()`, `plot()`, iteration
- **`ge.mesh(profile, dz)`** — calculation grid
- **`ge.log_plot(boreholes)`** — multi-borehole log plot
- *Reference:* Das (2010) Ch. 5; Terzaghi (1943).

#### New Module — `ge.design.stress` (Effective Stress)
- **`ge.total_stress(gamma, depth)`** — sum(gamma_i * H_i)
- **`ge.pore_pressure(z, z_w)`** — hydrostatic u
- **`ge.effective_stress(sigma, u)`** — Terzaghi principle
- **`ge.capillary_rise(D10, e)`** — Hazen empirical
- *Reference:* Terzaghi (1943); Das (2010) Ch. 5.

#### New Module — `ge.design.seepage` (Seepage & Flow)
- **`ge.darcy_flow(k, i, A)`** — Q = kiA
- **`ge.hydraulic_gradient(dh, L)`**
- **`ge.critical_gradient(Gs, e)`** — quicksand/heave check
- **`ge.equivalent_k(k_layers, H_layers, direction)`** — horizontal & vertical
- **`ge.flow_net(Nf, Nd, k, dh)`**
- *Reference:* Darcy (1856); Das (2010) Ch. 5.

#### New Module — `ge.design.boussinesq` (Stress Distribution)
- **`ge.boussinesq_point/line/strip/circular/rect()`** — full Boussinesq family
- **`ge.westergaard_point()`** — stratified medium
- **`ge.newmark_influence(m, n)`** — Fadum/Newmark influence factor
- **`ge.stress_2to1(q, B, L, z)`** — engineering approximation
- **`ge.stress_bulb()`** — visualization grid
- *Reference:* Boussinesq (1885), Westergaard (1938), Newmark (1935), Fadum (1948).

#### New Module — `ge.design.bearing` (Bearing Capacity, 4 methods)
- **`ge.bearing_factors(phi, method)`** — Nc, Nq, N_gamma for Terzaghi / Meyerhof / Hansen / Vesic
- **`ge.bearing_capacity()`** — unified solver with shape/depth/inclination factors
- **`ge.bearing_shape_factors`**, **`ge.bearing_depth_factors`**, **`ge.bearing_inclination_factors`**
- **`ge.bearing_allowable(q_u, FS)`**
- *Reference:* Terzaghi (1943); Meyerhof (1963); Hansen (1970); Vesic (1973); Das (2014) Ch. 4.

#### New Module — `ge.design.settlement` (Settlement, 5 methods)
- **`ge.settlement_immediate()`** — Janbu / Steinbrenner elastic settlement
- **`ge.settlement_primary()`** — 1-D Terzaghi (NC / OC-below / OC-crossing)
- **`ge.settlement_secondary()`** — Mesri & Godlewski Cα log(t)
- **`ge.settlement_schmertmann()`** — strain influence diagram for sands
- **`ge.time_factor(U)`**, **`ge.consolidation_degree(Tv)`**, **`ge.consolidation_time()`**
- *Reference:* Terzaghi (1925); Janbu et al. (1956); Schmertmann (1970, 1978); Mesri & Godlewski (1977).

#### New Module — `ge.design.earth_pressure` (Lateral Earth Pressure)
- **`ge.K0()`** — Jaky (1944) and Mayne–Kulhawy (1982)
- **`ge.Ka()`**, **`ge.Kp()`** — Rankine and Coulomb formulations
- **`ge.earth_pressure()`** — distribution with surcharge + water + cohesion
- **`ge.tension_crack_depth()`** — for cohesive soils
- **`ge.earth_pressure_plot()`**
- *Reference:* Rankine (1857); Coulomb (1776); Jaky (1944); Mayne & Kulhawy (1982).

#### New Module — `ge.design.walls` (Retaining Walls)
- **`ge.wall_overturning()`** — FS = M_R / M_O
- **`ge.wall_sliding()`** — FS with friction + cohesion + passive resistance
- **`ge.wall_bearing()`** — eccentric loading, kern check
- **`ge.sheet_pile()`** — cantilever sheet pile (Blum's method)
- *Reference:* Das (2014) Ch. 8; Blum (1931).

#### New Module — `ge.design.piles` (Pile Design)
- **`ge.pile_end_bearing()`** — Meyerhof, Vesic, Skempton clay
- **`ge.pile_skin_friction()`** — α (Tomlinson) / β (Burland) / λ (Vijayvergiya–Focht)
- **`ge.pile_capacity()`** — Q_ult = Q_p + Q_s + FS
- **`ge.pile_group_efficiency()`** — Converse–Labarre and Feld
- **`ge.pile_settlement()`** — Vesic 3-component method
- *Reference:* Meyerhof (1976); Tomlinson (1957); Burland (1973); Vijayvergiya & Focht (1972); Vesic (1977); Converse & Labarre (1947).

#### New Module — `ge.design.slopes` (Slope Stability, 4 methods)
- **`ge.infinite_slope()`** — with and without seepage parallel to slope
- **`ge.culmann()`** — planar failure critical height
- **`ge.taylor_stability()`** — Taylor's stability number chart fit
- **`ge.bishop()`** — Bishop's simplified method (iterative)
- *Reference:* Culmann (1875); Taylor (1937); Bishop (1955); Das (2010) Ch. 13.

#### New Module — `ge.dynamics` (Soil Dynamics & Seismic)
- **`ge.gmax(Vs, gamma)`** — Gmax = ρVs²
- **`ge.gmax_hardin(e, sigma_m, OCR, PI, soil_type)`** — Hardin & Black (1968) / Hardin & Drnevich (1972)
- **`ge.modulus_reduction()`** — G/Gmax curves (Darendeli 2001, Vucetic–Dobry 1991)
- **`ge.damping_ratio()`** — Darendeli (2001) with Hardin–Drnevich hyperbolic damping
- **`ge.depth_reduction(z, Mw)`** — rd from Idriss (1999), Liao–Whitman (1986), Cetin et al. (2004)
- **`ge.liquefaction_csr(amax, sigma_v, sigma_v_eff, Mw, z)`** — Seed & Idriss (1971)
- **`ge.liquefaction_crr()`** — CRR from SPT (Youd et al. 2001, Idriss–Boulanger 2008), CPT (IB2008), or Vs (Andrus & Stokoe 2000)
- **`ge.magnitude_scaling_factor(Mw)`** — Idriss (1999), NCEER, Boulanger–Idriss (2014)
- **`ge.liquefaction_fos(CSR, CRR, Mw)`** — full FS_L = CRR·MSF·K_σ·K_α / CSR
- **`ge.response_spectrum(T, Sa)`** — plot helper
- *Reference:* Hardin & Black (1968); Darendeli (2001); Seed & Idriss (1971); Youd et al. (2001); Idriss & Boulanger (2008); Andrus & Stokoe (2000).

#### New Module — `ge.io` (Data Readers)
- **`ge.read_csv(path, units_row)`** — geotech CSV with auto-header detection
- **`ge.read_ags(path)`** — AGS4 (UK ground-investigation standard)
- **`ge.read_gef(path)`** — GEF-CPT (Dutch format)
- **`ge.CPT(depth, qc, fs, u2)`** + `CPT.from_gef()` / `CPT.from_ags()` — container class
- *Reference:* AGS (2017) v4.0.4 spec; NEN GEF-CPT-Report spec.

### Removed
- Orphan `src/geoeq/insitu/` directory (superseded by `src/geoeq/site/` in v0.0.7;
  was never re-exported and contained no tests).

### Changed
- `pyproject.toml` version bumped to `0.1.0`.
- `geoeq/__init__.py` now exports **170+ top-level functions / classes** across 9 submodules
  (`soil`, `lab`, `site`, `viz`, `core`, `profile`, `design`, `dynamics`, `io`).

### Testing
- **563 tests passing** (up from 404 in v0.0.7).
- New test files: `test_profile.py`, `test_design_stress.py`, `test_design_boussinesq.py`,
  `test_design_bearing.py`, `test_design_settlement.py`, `test_design_earth_walls.py`,
  `test_design_piles_slopes.py`, `test_dynamics.py`, `test_io.py`.
- Textbook-value cross-checks against Das (2010, 2014), e.g.:
  - Fadum's table 6.5 (m=1, n=1 -> I=0.1752) ✓
  - Boussinesq circular at z=R (sigma/q = 0.6465) ✓
  - Meyerhof Nc, Nq, N_gamma at phi=30 (30.14, 18.40, 15.67) ✓
  - Hansen N_gamma at phi=30 (14.39) ✓
  - Westergaard at r=0, mu=0 (2P/πz²) ✓
  - Rankine Ka/Kp at phi=30 (1/3, 3) ✓

---

## [0.0.7] - 2026-05-03

### Added — Complete Site Investigation Module (`ge.site`)

#### New Module — `ge.site.spt` (Standard Penetration Test)
- **`ge.spt_n60(N, ER, Cb, Cs, Cr)`** — Energy-corrected blow count N₆₀ [Skempton, 1986]
- **`ge.spt_n160(N60, sigma_v, method)`** — Overburden-corrected (N₁)₆₀; 3 methods: Liao & Whitman, Skempton, Peck
- **`ge.spt_n160cs(N160, FC)`** — Fines-content correction for liquefaction [Youd et al., 2001]
- **`ge.spt_friction_angle(N, sigma_v, method)`** — φ' from SPT; 3 methods: Hatanaka, Kulhawy, Peck
- **`ge.spt_su(N60, method, PI)`** — Undrained shear strength; 2 methods: Stroud, Hara
- **`ge.spt_dr(N160, method)`** — Relative density; 3 methods: Meyerhof, Skempton, Kulhawy
- **`ge.spt_modulus(N, soil_type)`** — Elastic modulus from Bowles (1996); 6 soil types

#### New Module — `ge.site.cpt` (Cone Penetration Test)
- **`ge.cpt_normalize(qt, fs, sigma_v, sigma_v_eff, u2, u0)`** — Normalized Qt, Fr, Bq
- **`ge.cpt_ic(Qt, Fr)`** — Soil Behaviour Type Index Ic [Robertson, 2009]
- **`ge.cpt_sbt(Qt, Fr)`** — SBT zone classification (zones 2–7)
- **`ge.cpt_friction_angle(qt, sigma_v_eff)`** — φ' from CPT [Robertson & Campanella, 1983]
- **`ge.cpt_su(qt, sigma_v, Nkt)`** — Undrained shear strength [Lunne et al., 1997]
- **`ge.cpt_dr(qc, sigma_v_eff, C0, C1, C2)`** — Relative density [Baldi et al., 1986]
- **`ge.cpt_modulus(qt, sigma_v, Ic, alpha)`** — Constrained modulus [Robertson, 2009]
- **`ge.cpt_sbt_plot(Qt, Fr)`** — Robertson SBT chart visualization

#### New Module — `ge.site.vane` (Field Vane Shear Test)
- **`ge.vane_su(T, D, H)`** — Su from torque [ASTM D2573]
- **`ge.vane_correction(Su_fvt, PI)`** — Bjerrum (1972) correction factor
- **`ge.vane_remolded(T_peak, T_remolded, D, H)`** — Sensitivity and classification

#### New Module — `ge.site.pmt` (Pressuremeter Test)
- **`ge.pmt_parameters(pressure, volume)`** — Extract Em and pL from p–V curve
- **`ge.pmt_modulus(Em, alpha)`** — Young's modulus from PMT [Ménard, 1975]
- **`ge.pmt_su(pL, p0, sigma_h0)`** — Su from limit pressure [Baguelin et al., 1978]
- **`ge.pmt_bearing(pL, p0, sigma_v, shape)`** — Bearing capacity (Ménard method)
- **`ge.pmt_settlement(q, B, Em, alpha, shape)`** — Settlement estimate
- **`ge.pmt_ko(p0, sigma_v, u0)`** — K₀ from in-situ horizontal stress

#### New Module — `ge.site.plt_test` (Plate Load Test)
- **`ge.plt_bearing(pressure, settlement, failure_criterion)`** — Bearing capacity; 3 methods: log-log, tangent, settlement
- **`ge.plt_subgrade_modulus(pressure, settlement)`** — Modulus of subgrade reaction ks
- **`ge.plt_settlement_correction(s_plate, B_plate, B_foundation, soil_type)`** — Terzaghi & Peck size correction
- **`ge.plt_elastic_modulus(pressure, settlement, plate_diameter)`** — E from PLT
- **`ge.plt_plot(pressure, settlement)`** — Pressure–settlement curve

#### New Module — `ge.site.pile_load` (Pile Load Tests)
- **`ge.davisson(load, settlement, diameter, length, area, E)`** — Davisson's Offset Limit (1972)
- **`ge.chin(load, settlement)`** — Chin's Hyperbolic Method (1970)
- **`ge.de_beer(load, settlement)`** — De Beer's Double-Log Method (1967)
- **`ge.hansen_80(load, settlement)`** — Hansen's 80% Method (1963)
- **`ge.fhwa_5_percent(load, settlement, diameter)`** — FHWA 5% diameter criterion
- **`ge.case_method(F1, F2, v1, v2, impedance, Jc)`** — Case Method dynamic capacity
- **`ge.hiley(W, h, eta, s, c)`** — Hiley Driving Formula
- **`ge.danish_formula(W, h, eta, s, L, A, E)`** — Danish Driving Formula
- **`ge.enr(W, h, s)`** — Engineering News Record Formula
- **`ge.pile_impedance(E, A, c)`** — Pile impedance Z = EA/c
- **`ge.beta_integrity(measured, known)`** — Pile integrity beta ratio

#### New Module — `ge.site.field_perm` (Field Permeability)
- **`ge.slug_test(r, R, Le, T0)`** — Hvorslev (1951) method
- **`ge.pumping_test_confined(Q, h1, h2, r1, r2, H)`** — Thiem equation (confined)
- **`ge.pumping_test_unconfined(Q, h1, h2, r1, r2)`** — Thiem equation (unconfined)
- **`ge.lefranc_test(Q, H, D)`** — Lefranc borehole permeability test

#### New Module — `ge.site.field_cbr` (Field CBR / DCP)
- **`ge.dcp_cbr(dcpi, method)`** — CBR from DCP; 3 methods: Webster, TRL, Kleyn
- **`ge.field_cbr_test(penetration, load)`** — In-situ CBR from load–penetration

### Changed
- `pyproject.toml` version bumped to `0.0.7`
- `pyproject.toml` keywords expanded with SPT, CPT, pressuremeter, pile load test, DCP, site investigation
- `geoeq/__init__.py` exports all 46 new site functions at top level (92 total)

### Testing
- **404 tests** (up from 240) covering:
  - SPT: energy correction, 3 overburden methods, fines correction, friction angle (3 methods), Su (2 methods), Dr (3 methods), modulus (6 soil types)
  - CPT: normalization, Ic calculation, SBT zone classification, friction angle, Su, Dr, modulus, SBT plot
  - Vane: torque-to-Su conversion, Bjerrum correction, sensitivity classification
  - PMT: parameter extraction, modulus, Su, bearing capacity, settlement, K₀
  - PLT: 3 failure criteria, subgrade modulus, Terzaghi & Peck size correction, elastic modulus
  - Pile load: Davisson offset, Chin hyperbolic, De Beer log-log, Hansen 80%, FHWA 5%, Case method, Hiley, Danish, ENR, impedance, integrity
  - Field permeability: Hvorslev slug test, Thiem (confined/unconfined), Lefranc
  - Field CBR: 3 DCP–CBR correlations, in-situ CBR at 2.54/5.08 mm

---

## [0.0.6] - 2026-05-01

### Added — Complete Laboratory Testing Module (`ge.lab`)

#### New Module — `ge.lab.shear` (Shear Strength Tests)
- **`ge.direct_shear(normal_stress, shear_stress)`** — Mohr–Coulomb fitting from direct shear data → c, φ, R² (ASTM D3080)
- **`ge.triaxial(sigma3, delta_sigma, kind)`** — UU / CU / CD triaxial processing → c, φ (ASTM D2850/D4767)
- **`ge.unconfined(qu)`** — Unconfined compression → Su = qu/2 with consistency classification (ASTM D2166)
- **`ge.mohr_circle(sigma1, sigma3)`** — Draw Mohr circles with automatic failure envelope fitting
- **`ge.direct_shear_plot(normal_stress, shear_stress)`** — σ–τ failure envelope visualization

#### New Module — `ge.lab.consolidation` (Consolidation Tests)
- **`ge.oedometer(stress, void_ratio)`** — Process 1-D consolidation data → Cc, Cr (ASTM D2435)
- **`ge.preconsolidation(stress, void_ratio, method)`** — Casagrande method for preconsolidation pressure pc
- **`ge.compression_index(method, **params)`** — 6 empirical correlations: Terzaghi, Skempton, Rendon, Nishida, Nagaraj, Hough
- **`ge.cv(time, deformation, method, H_dr)`** — Coefficient of consolidation via log-time or root-time method
- **`ge.oedometer_plot(stress, void_ratio)`** — e–log(σ') curve with Cc, Cr, and pc annotation
- **`ge.cv_plot(time, deformation, method)`** — Time–deformation plot with cv extraction

#### New Module — `ge.lab.compaction` (Compaction Tests)
- **`ge.proctor(water_content, dry_density)`** — Standard/Modified Proctor → w_opt, γ_d_max (ASTM D698/D1557)
- **`ge.zav_line(Gs, w_range)`** — Zero Air Voids line at S = 100% (Das Eq. 5.12)
- **`ge.saturation_line(Gs, S, w_range)`** — Constant-saturation contour on compaction plot
- **`ge.relative_compaction(gamma_d, gamma_d_max)`** — RC = γ_d / γ_d_max × 100%
- **`ge.proctor_plot(water_content, dry_density, Gs)`** — Compaction curve with ZAV line and S = 60%/80% contours

#### New Module — `ge.lab.permeability` (Permeability Tests)
- **`ge.constant_head(Q, L, A, h, t)`** — Constant-head permeability k (ASTM D2434, Das Eq. 5.11)
- **`ge.falling_head(a, L, A, h1, h2, t)`** — Falling-head permeability k (ASTM D5084, Das Eq. 5.13)
- **`ge.permeability_plot()`** — Q vs i (constant-head) or ln(h) vs t (falling-head)

#### New Module — `ge.lab.atterberg_test` (Atterberg Limits Processing)
- **`ge.liquid_limit_test(blow_count, water_content)`** — Casagrande cup flow curve → LL at N = 25 (ASTM D4318)
- **`ge.flow_curve_plot(blow_count, water_content)`** — Semi-log flow curve with LL marker

#### New Module — `ge.lab.cbr` (California Bearing Ratio)
- **`ge.cbr_test(penetration, load)`** — CBR at 2.54 mm and 5.08 mm penetration (ASTM D1883)
- **`ge.cbr_plot(penetration, load)`** — Load–penetration curve with standard reference markers

### Changed
- `pyproject.toml` version bumped to `0.0.6`
- `pyproject.toml` keywords expanded with lab test terms
- `geoeq/__init__.py` exports all 21 new lab functions at top level (46 total)
- Import convention standardized to `import geoeq as ge` / `ge.` prefix across all documentation

### Testing
- **240 tests** (up from 144) covering:
  - Direct shear fitting with hand-calculated Mohr–Coulomb examples
  - Triaxial UU (φ = 0, Su) and CD/CU (p-q space fitting)
  - Unconfined compression consistency classification
  - Mohr circle plotting (smoke tests)
  - Oedometer Cc/Cr extraction from Das-style load step data
  - Casagrande preconsolidation pressure determination
  - All 6 compression index empirical correlations
  - cv from log-time and root-time methods
  - Proctor compaction peak fitting (w_opt, γ_d_max)
  - ZAV line and saturation line calculations
  - Relative compaction (scalar and array)
  - Constant-head and falling-head permeability (Das Eq. 5.11, 5.13)
  - Liquid limit flow curve fitting
  - CBR at 2.54 mm and 5.08 mm with standard reference loads
  - All visualization functions (smoke tests)

---

## [0.0.5] - 2026-05-01

### Added — Complete Soil Module (Chapters 3 & 4)

#### New Functions — Index Properties
- **`ge.specific_gravity(Ms, Vs)`** — Compute Gs from mass and volume of solids
- **`ge.activity(PI, clay_fraction)`** — Activity of clay (Skempton, 1953); classifies Inactive / Normal / Active
- **`ge.sensitivity(Su_undisturbed, Su_remolded)`** — Sensitivity of clay (St); classifies Insensitive through Quick clay
- **`ge.liquidity_index(w, PL, PI)`** — Standalone liquidity index convenience function

#### New Functions — Soil Classification
- **`ge.classify_uscs(LL, PL, gravel, sand, fines, Cu, Cc, organic)`**
  - Full USCS classification per ASTM D2487-17 flowchart
  - Supports all group symbols: GW, GP, GM, GC, SW, SP, SM, SC, CL, ML, CH, MH, OL, OH, CL-ML
  - Handles dual symbols for 5–12% fines zone (e.g. SW-SM, GP-GC)
  - Auto-computes fines from gravel + sand if not provided
  - Returns `{'symbol': str, 'name': str}`

- **`ge.classify_aashto(LL, PL, gravel, sand, fines)`**
  - Full AASHTO M145 classification with left-to-right elimination
  - Groups: A-1-a, A-1-b, A-3, A-2-4 through A-2-7, A-4 through A-7-6
  - Computes Group Index: GI = (F-35)[0.2+0.005(LL-40)] + 0.01(F-15)(PI-10)
  - Returns `{'group': str, 'group_index': int, 'description': str}`

- **`ge.plasticity_chart(LL, PL, labels, ax, save)`**
  - Casagrande plasticity chart with A-line and U-line
  - Plots soil samples on the chart for visual USCS classification
  - Supports custom axes, labels, and file export

#### Improvements — Existing Functions
- **All 11 soil functions** now have complete NumPy-style docstrings with:
  - LaTeX formulas with Das (2021) equation numbers
  - Full parameter descriptions with types, units, and valid ranges
  - Reference citations (Das, Holtz, Skempton, ASTM)
  - 2+ working usage examples per function
- **NumPy array support** — All functions accept both scalars and numpy arrays via `np.asarray()`; return `float` for scalar inputs and `ndarray` for array inputs
- **Validation upgraded** — `core/validation.py` now handles numpy arrays correctly using `np.any()`

### Changed
- `pyproject.toml` version bumped to `0.0.5`
- `pyproject.toml` keywords expanded with "soil classification", "USCS", "AASHTO", "Atterberg limits"
- `geoeq/__init__.py` exports all 7 new functions at top level

### Testing
- **144 tests** (up from ~70) covering:
  - Textbook hand-calculated examples (Das 2021)
  - Edge cases (e=0 solid rock, S=0 dry, S=1 saturated)
  - Invalid input validation (ValueError for physically impossible values)
  - NumPy array inputs (verify output shapes)
  - Roundtrip consistency checks (void_ratio ↔ porosity, saturation ↔ water_content)
  - All USCS soil types (GW through CL-ML, dual symbols)
  - AASHTO groups (A-1-a through A-7-6 with Group Index verification)
  - Plasticity chart smoke tests

---

## [0.0.4] - 2026-04-04

### Added — Particle Size Distribution Suite (Chapter 2)

#### Analysis Functions
- **`ge.sieve_ana(opening, mass_retained, standard, total_mass)`**
  - Full sieve analysis with automatic designation-to-mm mapping
  - Supports: ASTM E11 (`#4`, `#10`, `#200`, `3/4"`, etc.), BS 410, IS 460
  - Returns: `diameter`, `opening`, `mass_retained`, `percent_retained`, `cumulative_retained`, `percent_finer`, `total_mass`

- **`ge.hydro_ana(reading, time, T, Gs, Ws, Fm, Fz)`**
  - Full ASTM D7928 hydrometer sedimentation analysis
  - Implements Stokes' Law for equivalent particle diameter D
  - Applies meniscus correction (Fm), zero correction (Fz), and temperature corrections
  - Calibrated for the standard **152H hydrometer** model
  - Returns: `diameter`, `percent_finer`, `reading`, `depth`, `time`

#### Visualization — `ge.grain_size_plot()`
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
- **`ge.grain_d10(data)`** — D₁₀ effective size via log-linear interpolation
- **`ge.grain_d30(data)`** — D₃₀ via log-linear interpolation
- **`ge.grain_d60(data)`** — D₆₀ controlling size via log-linear interpolation
- **`ge.grain_Cu(data)`** — Uniformity Coefficient Cᵤ = D₆₀ / D₁₀
- **`ge.grain_Cc(data)`** — Coefficient of Curvature Cᶜ = D₃₀² / (D₁₀ × D₆₀)

#### Dependencies Added
- `numpy` — array operations for grain size calculations
- `matplotlib` — professional plot rendering
- `scipy` — PCHIP monotonic interpolation (`PchipInterpolator`)

#### Flat API Exports
All new functions are available directly under `import geoeq as ge` (no submodule path required).

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

- `ge.density()`: Unified density calculation function supporting variable state modes (`kind`), multiple output formats (`unit`), and basic (`mass`, `volume`) calculations.
- `ge.atterberg()`: Unified function for plasticity index (`PI`), liquidity index (`LI`), and consistency index (`CI`).

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
