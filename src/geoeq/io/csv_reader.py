"""
Generic CSV reader with auto-detection of headers and units rows.

Geotech CSV files often have:
  * 1 line of metadata (project name etc.)
  * 1 header line (depth, qc, fs, u2 ...)
  * 1 units line (m, MPa, kPa, kPa ...)
  * data

This reader scans the first ~20 lines, picks the row with the most
non-numeric tokens as the header, and returns a dict of columns plus
the raw DataFrame (if pandas is available).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np


def read_csv(path, sep: str = ",", skiprows: Optional[int] = None,
              encoding: str = "utf-8", units_row: bool = False) -> dict:
    """Read a geotech CSV file.

    Parameters
    ----------
    path : str or Path
        File path.
    sep : str
        Field separator. Default ','.
    skiprows : int, optional
        Number of header rows to skip. If None, auto-detected.
    encoding : str
        File encoding. Default 'utf-8'.
    units_row : bool
        If True, treat the row after the header as a units row (returned
        in the 'units' dict) and skip it for numeric parsing.

    Returns
    -------
    dict
        ``{'data': dict[col -> ndarray], 'df': pd.DataFrame|None,
        'units': dict|None, 'header_row': int}``
    """
    path = Path(path)
    lines = path.read_text(encoding=encoding).splitlines()

    # Auto-detect header row: first line whose every field, when stripped
    # of brackets, fails to parse as float.
    if skiprows is None:
        skiprows = 0
        for i, line in enumerate(lines[:20]):
            toks = [t.strip() for t in line.split(sep)]
            if not toks or all(t == "" for t in toks):
                continue
            nonnum = sum(1 for t in toks if not _is_number(t))
            if nonnum == len(toks) and len(toks) >= 2:
                skiprows = i
                break

    header = [t.strip() for t in lines[skiprows].split(sep)]
    data_start = skiprows + 1
    units = None
    if units_row and data_start < len(lines):
        units_toks = [t.strip() for t in lines[data_start].split(sep)]
        if not all(_is_number(t) for t in units_toks):
            units = dict(zip(header, units_toks))
            data_start += 1

    rows = []
    for line in lines[data_start:]:
        if not line.strip():
            continue
        toks = [t.strip() for t in line.split(sep)]
        if len(toks) != len(header):
            continue
        try:
            rows.append([float(t) if _is_number(t) else np.nan for t in toks])
        except ValueError:
            continue
    arr = np.array(rows, dtype=float) if rows else np.empty((0, len(header)))
    data = {col: arr[:, i] for i, col in enumerate(header)}

    df = None
    try:
        import pandas as pd
        df = pd.DataFrame(data)
    except ImportError:
        pass

    return {"data": data, "df": df, "units": units,
            "header_row": skiprows, "n_rows": len(rows)}


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False
