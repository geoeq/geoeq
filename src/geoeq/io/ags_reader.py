"""
AGS4 reader -- minimal but functional.

AGS4 is the UK ground-investigation data standard. A file is a UTF-8
text file with groups of quoted comma-separated lines:

    "GROUP","XXXX"
    "HEADING","LOCA_ID","LOCA_NATE","LOCA_NATN",...
    "UNIT","","m","m",...
    "TYPE","ID","2DP","2DP",...
    "DATA","BH1","123.45","678.90",...

This implementation extracts the data rows of each group into a dict
keyed by group name and returns column dicts.

Spec: AGS4 -- https://www.ags.org.uk/data-format/
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List


def read_ags(path, encoding: str = "utf-8") -> Dict[str, dict]:
    """Parse an AGS4 file and return all groups as dicts.

    Parameters
    ----------
    path : str or Path
        AGS file path.
    encoding : str
        File encoding (AGS4 specifies UTF-8). Default 'utf-8'.

    Returns
    -------
    dict
        ``{group_name: {'headings': [...], 'units': [...], 'data': [{...}, ...]}}``

    Reference
    ---------
    AGS (2017) -- AGS4 data-transfer specification.
    """
    path = Path(path)
    text = path.read_text(encoding=encoding, errors="replace")

    groups: Dict[str, dict] = {}
    current_group = None
    headings: List[str] = []
    units: List[str] = []
    rows: List[dict] = []

    for line_tokens in csv.reader(text.splitlines(), quotechar='"'):
        if not line_tokens:
            continue
        tag = line_tokens[0]
        if tag == "GROUP":
            # Commit previous group.
            if current_group is not None:
                groups[current_group] = {
                    "headings": headings, "units": units, "data": rows}
            current_group = line_tokens[1] if len(line_tokens) > 1 else None
            headings, units, rows = [], [], []
        elif tag == "HEADING":
            headings = list(line_tokens[1:])
        elif tag == "UNIT":
            units = list(line_tokens[1:])
        elif tag == "TYPE":
            # Field type row -- ignored for read.
            pass
        elif tag == "DATA":
            values = list(line_tokens[1:])
            row = dict(zip(headings, values))
            rows.append(row)

    # Commit last.
    if current_group is not None:
        groups[current_group] = {
            "headings": headings, "units": units, "data": rows}

    return groups
