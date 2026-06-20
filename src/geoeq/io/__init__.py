"""
ge.io -- Data readers for geotechnical file formats.

* ``read_csv``  -- CSV with auto-header detection.
* ``read_ags``  -- AGS4 (Association of Geotechnical & Geoenvironmental
                   Specialists) -- the UK / Australian / NZ ground-investigation
                   data standard.
* ``read_gef``  -- GEF-CPT (Dutch geotechnical data exchange format).
* ``CPT``       -- light-weight container class for a single CPT sounding.

References
----------
* AGS (2017). *Electronic transfer of geotechnical and geoenvironmental
  data*. Version 4.0.4. https://www.ags.org.uk/data-format/
* GEF (Geef Eenvoudige Formaat) specification, NEN.
"""

from geoeq.io.csv_reader import read_csv
from geoeq.io.ags_reader import read_ags
from geoeq.io.gef_reader import read_gef
from geoeq.io.cpt_container import CPT

__all__ = ["read_csv", "read_ags", "read_gef", "CPT"]
