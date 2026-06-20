"""
ge.dynamics -- Soil dynamics and seismic engineering.

Small-strain shear modulus, modulus reduction and damping curves
(Darendeli 2001), and the simplified-procedure liquefaction triggering
analysis (Seed-Idriss 1971 / Youd et al. 2001 / Idriss-Boulanger 2008).

References
----------
* Hardin, B. O., Black, W. L. (1968). "Vibration modulus of normally
  consolidated clays." *J. Soil Mech. Found. Eng.*, ASCE, 94(SM2), 353-369.
* Hardin, B. O., Drnevich, V. P. (1972). "Shear modulus and damping in
  soils: design equations and curves." *JSMFE*, 98(SM7), 667-692.
* Darendeli, M. B. (2001). *Development of a new family of normalized
  modulus reduction and material damping curves*. PhD thesis, UT Austin.
* Seed, H. B., Idriss, I. M. (1971). "Simplified procedure for evaluating
  soil liquefaction potential." *J. Soil Mech. Found. Eng.*, 97(SM9).
* Youd, T. L., Idriss, I. M., et al. (2001). "Liquefaction resistance of
  soils: summary report from 1996 NCEER and 1998 NCEER/NSF workshops on
  evaluation of liquefaction resistance of soils." *JGGE*, 127(10), 817-833.
* Idriss, I. M., Boulanger, R. W. (2008). *Soil liquefaction during
  earthquakes*. EERI Monograph MNO-12.
"""

from geoeq.dynamics.modulus import (
    gmax,
    gmax_hardin,
    modulus_reduction,
    damping_ratio,
)
from geoeq.dynamics.liquefaction import (
    depth_reduction,
    liquefaction_csr,
    liquefaction_crr,
    liquefaction_fos,
    magnitude_scaling_factor,
)
from geoeq.dynamics.response import response_spectrum
from geoeq.dynamics.plots import gmax_curves_plot, liquefaction_chart

__all__ = [
    "gmax", "gmax_hardin",
    "modulus_reduction", "damping_ratio",
    "depth_reduction",
    "liquefaction_csr", "liquefaction_crr", "liquefaction_fos",
    "magnitude_scaling_factor",
    "response_spectrum",
    "gmax_curves_plot", "liquefaction_chart",
]
