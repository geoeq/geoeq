"""
ge.design -- Geotechnical engineering analysis & design.

The big module: stress, bearing, settlement, earth pressure, piles,
slopes, retaining walls, seepage.

All formulas cite their textbook source or original paper in the
function docstring. SI units throughout (kN, kPa, m, degrees) unless
explicitly noted.
"""

from geoeq.design.stress import (
    total_stress,
    pore_pressure,
    effective_stress,
    capillary_rise,
    stress_plot,
)
from geoeq.design.seepage import (
    darcy_flow,
    hydraulic_gradient,
    critical_gradient,
    equivalent_k,
    flow_net,
)
from geoeq.design.boussinesq import (
    boussinesq_point,
    boussinesq_line,
    boussinesq_strip,
    boussinesq_circular,
    boussinesq_rect,
    westergaard_point,
    newmark_influence,
    stress_2to1,
    stress_bulb,
)
from geoeq.design.bearing import (
    bearing_factors,
    bearing_capacity,
    bearing_allowable,
    bearing_shape_factors,
    bearing_depth_factors,
    bearing_inclination_factors,
)
from geoeq.design.settlement import (
    settlement_immediate,
    settlement_primary,
    settlement_secondary,
    settlement_schmertmann,
    time_factor,
    consolidation_degree,
    consolidation_time,
)
from geoeq.design.earth_pressure import (
    K0,
    Ka,
    Kp,
    earth_pressure,
    tension_crack_depth,
    earth_pressure_plot,
)
from geoeq.design.walls import (
    wall_overturning,
    wall_sliding,
    wall_bearing,
    sheet_pile,
)
from geoeq.design.piles import (
    pile_end_bearing,
    pile_skin_friction,
    pile_capacity,
    pile_settlement,
    pile_group_efficiency,
)
from geoeq.design.slopes import (
    infinite_slope,
    culmann,
    taylor_stability,
    bishop,
)
from geoeq.design.plots import (
    bearing_capacity_plot,
    settlement_time_plot,
    stress_isobar_plot,
    taylor_chart_plot,
)

__all__ = [
    # Stress
    "total_stress", "pore_pressure", "effective_stress",
    "capillary_rise", "stress_plot",
    # Seepage
    "darcy_flow", "hydraulic_gradient", "critical_gradient",
    "equivalent_k", "flow_net",
    # Stress distribution
    "boussinesq_point", "boussinesq_line", "boussinesq_strip",
    "boussinesq_circular", "boussinesq_rect", "westergaard_point",
    "newmark_influence", "stress_2to1", "stress_bulb",
    # Bearing
    "bearing_factors", "bearing_capacity", "bearing_allowable",
    "bearing_shape_factors", "bearing_depth_factors",
    "bearing_inclination_factors",
    # Settlement
    "settlement_immediate", "settlement_primary", "settlement_secondary",
    "settlement_schmertmann", "time_factor", "consolidation_degree",
    "consolidation_time",
    # Earth pressure
    "K0", "Ka", "Kp", "earth_pressure",
    "tension_crack_depth", "earth_pressure_plot",
    # Walls
    "wall_overturning", "wall_sliding", "wall_bearing", "sheet_pile",
    # Piles
    "pile_end_bearing", "pile_skin_friction", "pile_capacity",
    "pile_settlement", "pile_group_efficiency",
    # Slopes
    "infinite_slope", "culmann", "taylor_stability", "bishop",
    # Plots
    "bearing_capacity_plot", "settlement_time_plot",
    "stress_isobar_plot", "taylor_chart_plot",
]
