"""
GEF-CPT reader (Dutch geotechnical exchange format).

A GEF file has a header section of ``#KEY= value`` lines, ending with
``#EOH=``, followed by columnar data. For CPT files the relevant header
keys are:

* ``#COLUMN`` -- number of data columns
* ``#COLUMNINFO=  index, unit, name, quantity`` -- one per column
* ``#COLUMNSEPARATOR`` -- field separator
* ``#COLUMNVOID``  -- placeholder values to treat as NaN
* ``#RECORDSEPARATOR`` -- line separator (often newline)

Reference: NEN 5104 / GEF-CPT-Report -- https://publicwiki.deltares.nl/
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np


def read_gef(path, encoding: str = "latin-1") -> Dict:
    """Read a Dutch GEF-CPT file.

    Parameters
    ----------
    path : str or Path
        Path to the .gef file.
    encoding : str
        File encoding (most GEF files are 'latin-1'). Default.

    Returns
    -------
    dict
        ``{'header': {...}, 'columns': [...], 'units': [...],
        'data': {colname: ndarray, ...}}``.

    Reference
    ---------
    NEN -- GEF file specification (CPT-Report).
    """
    path = Path(path)
    text = path.read_text(encoding=encoding, errors="replace")

    header: Dict[str, str] = {}
    column_info: List[dict] = []
    sep = " "
    void_values: List[float] = []
    in_header = True
    data_lines: List[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if in_header:
            if line.startswith("#EOH"):
                in_header = False
                continue
            if line.startswith("#"):
                # parse #KEY= value
                if "=" in line:
                    key, _, value = line[1:].partition("=")
                    key = key.strip().upper()
                    value = value.strip()
                    if key == "COLUMNSEPARATOR":
                        sep = value if value else " "
                    elif key == "COLUMNINFO":
                        toks = [t.strip() for t in value.split(",")]
                        if len(toks) >= 3:
                            try:
                                idx = int(toks[0])
                            except ValueError:
                                idx = len(column_info) + 1
                            column_info.append({
                                "index": idx,
                                "unit": toks[1] if len(toks) > 1 else "",
                                "name": toks[2] if len(toks) > 2 else f"col{idx}",
                                "quantity": toks[3] if len(toks) > 3 else "",
                            })
                    elif key == "COLUMNVOID":
                        toks = [t.strip() for t in value.split(",")]
                        try:
                            void_values.append(float(toks[-1]))
                        except (ValueError, IndexError):
                            pass
                    else:
                        header[key] = value
        else:
            if line:
                data_lines.append(line)

    # Parse data
    rows = []
    for line in data_lines:
        parts = [p for p in line.replace("\t", " ").split(sep) if p] \
            if sep != " " else line.split()
        try:
            vals = [float(p) for p in parts]
            rows.append(vals)
        except ValueError:
            continue

    arr = np.array(rows, dtype=float) if rows else \
        np.empty((0, len(column_info)))
    # Replace void values with NaN.
    for v in void_values:
        arr[arr == v] = np.nan

    column_info.sort(key=lambda c: c["index"])
    names = [c["name"] for c in column_info]
    units = [c["unit"] for c in column_info]
    data = {n: arr[:, i] for i, n in enumerate(names) if i < arr.shape[1]}

    return {
        "header": header, "columns": names, "units": units,
        "data": data, "column_info": column_info,
    }
