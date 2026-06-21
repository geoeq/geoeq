"""
GeoEq -- A Python library for geotechnical engineering calculations.

A clean, validated, MIT-licensed library covering the complete
onshore geotechnical workflow under one flat namespace:
soil properties, lab testing, site investigation, layered ground
modelling, engineering design, soil dynamics, and data I/O.

    import geoeq as ge
    ge.bearing_capacity(c=10, gamma=18, Df=1, B=2, phi=30)

Homepage:      https://geoeq.github.io
Repository:    https://github.com/geoeq/geoeq
Documentation: https://geoeq.github.io/user-guide.html
License:       MIT
"""

__version__ = "0.1.2"
__author__ = "Ripon Chandra Malo"
__license__ = "MIT"

from geoeq.soil.properties import (
    void_ratio,
    porosity,
    saturation,
    water_content,
    specific_gravity,
    density,
    relative_density,
    atterberg,
    activity,
    sensitivity,
    liquidity_index,
)
from geoeq.soil.classification import (
    classify_uscs,
    classify_aashto,
    plasticity_chart,
)
from geoeq.soil.sieve import sieve_ana
from geoeq.soil.hydrometer import hydro_ana
from geoeq.soil.grain_size import (
    grain_d10,
    grain_d30,
    grain_d60,
    grain_Cu,
    grain_Cc,
    process_combined_data,
)
from geoeq.viz.grain_size import grain_size_plot

# Lab — Shear strength
from geoeq.lab.shear import (
    direct_shear,
    triaxial,
    unconfined,
    mohr_circle,
    direct_shear_plot,
)

# Lab — Consolidation
from geoeq.lab.consolidation import (
    oedometer,
    preconsolidation,
    compression_index,
    cv,
    oedometer_plot,
    cv_plot,
)

# Lab — Compaction
from geoeq.lab.compaction import (
    proctor,
    zav_line,
    saturation_line,
    relative_compaction,
    proctor_plot,
)

# Lab — Permeability
from geoeq.lab.permeability import (
    constant_head,
    falling_head,
    permeability_plot,
)

# Lab — Atterberg limits test
from geoeq.lab.atterberg_test import (
    liquid_limit_test,
    flow_curve_plot,
)

# Lab — CBR
from geoeq.lab.cbr import (
    cbr_test,
    cbr_plot,
)

# Site — SPT
from geoeq.site.spt import (
    spt_n60,
    spt_n160,
    spt_n160cs,
    spt_friction_angle,
    spt_su,
    spt_dr,
    spt_modulus,
)

# Site — CPT
from geoeq.site.cpt import (
    cpt_normalize,
    cpt_ic,
    cpt_sbt,
    cpt_friction_angle,
    cpt_su,
    cpt_dr,
    cpt_modulus,
    cpt_sbt_plot,
)

# Site — Vane shear
from geoeq.site.vane import (
    vane_su,
    vane_correction,
    vane_remolded,
)

# Site — Pressuremeter
from geoeq.site.pmt import (
    pmt_parameters,
    pmt_modulus,
    pmt_su,
    pmt_bearing,
    pmt_settlement,
    pmt_ko,
)

# Site — Plate load test
from geoeq.site.plt_test import (
    plt_bearing,
    plt_subgrade_modulus,
    plt_settlement_correction,
    plt_elastic_modulus,
    plt_plot,
)

# Site — Pile load testing
from geoeq.site.pile_load import (
    davisson,
    chin,
    de_beer,
    hansen_80,
    fhwa_5_percent,
    case_method,
    hiley,
    danish_formula,
    enr,
    pile_impedance,
    beta_integrity,
)

# Site — Field permeability
from geoeq.site.field_perm import (
    slug_test,
    pumping_test_confined,
    pumping_test_unconfined,
    lefranc_test,
)

# Site — Field CBR / DCP
from geoeq.site.field_cbr import (
    dcp_cbr,
    field_cbr_test,
)

# --------------------------------------------------------------------
# Profile (v0.1.0) -- layered soil profiles
# --------------------------------------------------------------------
from geoeq.profile import Soil, SoilProfile, mesh, log_plot

# --------------------------------------------------------------------
# Design (v0.1.0) -- effective stress, Boussinesq, bearing, settlement,
# earth pressure, walls, piles, slopes, seepage
# --------------------------------------------------------------------
from geoeq.design import (
    # Stress
    total_stress, pore_pressure, effective_stress,
    capillary_rise, stress_plot,
    # Seepage
    darcy_flow, hydraulic_gradient, critical_gradient,
    equivalent_k, flow_net,
    # Stress distribution
    boussinesq_point, boussinesq_line, boussinesq_strip,
    boussinesq_circular, boussinesq_rect, westergaard_point,
    newmark_influence, stress_2to1, stress_bulb,
    # Bearing
    bearing_factors, bearing_capacity, bearing_allowable,
    bearing_shape_factors, bearing_depth_factors,
    bearing_inclination_factors,
    # Settlement
    settlement_immediate, settlement_primary, settlement_secondary,
    settlement_schmertmann, time_factor, consolidation_degree,
    consolidation_time,
    # Earth pressure
    K0, Ka, Kp, earth_pressure,
    tension_crack_depth, earth_pressure_plot,
    # Walls
    wall_overturning, wall_sliding, wall_bearing, sheet_pile,
    # Piles (design side)
    pile_end_bearing, pile_skin_friction, pile_capacity,
    pile_settlement, pile_group_efficiency,
    # Slopes
    infinite_slope, culmann, taylor_stability, bishop,
    # Design plots
    bearing_capacity_plot, settlement_time_plot,
    stress_isobar_plot, taylor_chart_plot,
)

# --------------------------------------------------------------------
# Dynamics (v0.1.0) -- Gmax, modulus reduction, damping, liquefaction
# --------------------------------------------------------------------
from geoeq.dynamics import (
    gmax, gmax_hardin,
    modulus_reduction, damping_ratio,
    depth_reduction,
    liquefaction_csr, liquefaction_crr, liquefaction_fos,
    magnitude_scaling_factor, response_spectrum,
    gmax_curves_plot, liquefaction_chart,
)

# --------------------------------------------------------------------
# I/O (v0.1.0) -- CSV, AGS4, GEF, CPT container
# --------------------------------------------------------------------
from geoeq.io import read_csv, read_ags, read_gef, CPT

__all__ = [
    # Phase relations
    "void_ratio",
    "porosity",
    "saturation",
    "water_content",
    "specific_gravity",
    # Density
    "density",
    "relative_density",
    # Atterberg & index properties
    "atterberg",
    "activity",
    "sensitivity",
    "liquidity_index",
    # Classification
    "classify_uscs",
    "classify_aashto",
    "plasticity_chart",
    # Grain size
    "sieve_ana",
    "hydro_ana",
    "grain_d10",
    "grain_d30",
    "grain_d60",
    "grain_Cu",
    "grain_Cc",
    "process_combined_data",
    "grain_size_plot",
    # Lab — Shear strength
    "direct_shear",
    "triaxial",
    "unconfined",
    "mohr_circle",
    "direct_shear_plot",
    # Lab — Consolidation
    "oedometer",
    "preconsolidation",
    "compression_index",
    "cv",
    "oedometer_plot",
    "cv_plot",
    # Lab — Compaction
    "proctor",
    "zav_line",
    "saturation_line",
    "relative_compaction",
    "proctor_plot",
    # Lab — Permeability
    "constant_head",
    "falling_head",
    "permeability_plot",
    # Lab — Atterberg test
    "liquid_limit_test",
    "flow_curve_plot",
    # Lab — CBR
    "cbr_test",
    "cbr_plot",
    # Site — SPT
    "spt_n60",
    "spt_n160",
    "spt_n160cs",
    "spt_friction_angle",
    "spt_su",
    "spt_dr",
    "spt_modulus",
    # Site — CPT
    "cpt_normalize",
    "cpt_ic",
    "cpt_sbt",
    "cpt_friction_angle",
    "cpt_su",
    "cpt_dr",
    "cpt_modulus",
    "cpt_sbt_plot",
    # Site — Vane shear
    "vane_su",
    "vane_correction",
    "vane_remolded",
    # Site — Pressuremeter
    "pmt_parameters",
    "pmt_modulus",
    "pmt_su",
    "pmt_bearing",
    "pmt_settlement",
    "pmt_ko",
    # Site — Plate load test
    "plt_bearing",
    "plt_subgrade_modulus",
    "plt_settlement_correction",
    "plt_elastic_modulus",
    "plt_plot",
    # Site — Pile load testing
    "davisson",
    "chin",
    "de_beer",
    "hansen_80",
    "fhwa_5_percent",
    "case_method",
    "hiley",
    "danish_formula",
    "enr",
    "pile_impedance",
    "beta_integrity",
    # Site — Field permeability
    "slug_test",
    "pumping_test_confined",
    "pumping_test_unconfined",
    "lefranc_test",
    # Site — Field CBR / DCP
    "dcp_cbr",
    "field_cbr_test",
    # ---- v0.1.0 ----
    # Profile
    "Soil", "SoilProfile", "mesh", "log_plot",
    # Design — Stress
    "total_stress", "pore_pressure", "effective_stress",
    "capillary_rise", "stress_plot",
    # Design — Seepage
    "darcy_flow", "hydraulic_gradient", "critical_gradient",
    "equivalent_k", "flow_net",
    # Design — Boussinesq / stress distribution
    "boussinesq_point", "boussinesq_line", "boussinesq_strip",
    "boussinesq_circular", "boussinesq_rect", "westergaard_point",
    "newmark_influence", "stress_2to1", "stress_bulb",
    # Design — Bearing capacity
    "bearing_factors", "bearing_capacity", "bearing_allowable",
    "bearing_shape_factors", "bearing_depth_factors",
    "bearing_inclination_factors",
    # Design — Settlement
    "settlement_immediate", "settlement_primary", "settlement_secondary",
    "settlement_schmertmann", "time_factor", "consolidation_degree",
    "consolidation_time",
    # Design — Earth pressure
    "K0", "Ka", "Kp", "earth_pressure",
    "tension_crack_depth", "earth_pressure_plot",
    # Design — Walls
    "wall_overturning", "wall_sliding", "wall_bearing", "sheet_pile",
    # Design — Piles
    "pile_end_bearing", "pile_skin_friction", "pile_capacity",
    "pile_settlement", "pile_group_efficiency",
    # Design — Slopes
    "infinite_slope", "culmann", "taylor_stability", "bishop",
    # Dynamics
    "gmax", "gmax_hardin",
    "modulus_reduction", "damping_ratio",
    "depth_reduction",
    "liquefaction_csr", "liquefaction_crr", "liquefaction_fos",
    "magnitude_scaling_factor", "response_spectrum",
    "gmax_curves_plot", "liquefaction_chart",
    # I/O
    "read_csv", "read_ags", "read_gef", "CPT",
]
