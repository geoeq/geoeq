"""
CPT data container -- holds depth, qc, fs, u2 arrays plus methods to
apply normalization and classification.

This is a thin wrapper that composes the functional API from
``ge.site.cpt`` so users can write::

    cpt = CPT.from_gef("file.gef")
    cpt.normalize(sigma_v_profile)
    cpt.classify()
    cpt.plot_sbt()
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np


class CPT:
    """Container for a single CPT sounding.

    Attributes
    ----------
    depth : ndarray
        Depth (m).
    qc : ndarray
        Cone tip resistance (MPa or kPa -- callers must be consistent).
    fs : ndarray
        Sleeve friction (kPa).
    u2 : ndarray
        Pore pressure behind cone (kPa). Default zeros.
    title : str
        Sounding name.
    """

    def __init__(
        self,
        depth: Sequence[float],
        qc: Sequence[float],
        fs: Sequence[float] = None,
        u2: Sequence[float] = None,
        title: str = "CPT",
    ):
        self.depth = np.asarray(depth, dtype=float)
        self.qc = np.asarray(qc, dtype=float)
        self.fs = (np.asarray(fs, dtype=float)
                   if fs is not None else np.zeros_like(self.depth))
        self.u2 = (np.asarray(u2, dtype=float)
                   if u2 is not None else np.zeros_like(self.depth))
        self.title = title
        self._normalized = None

    @classmethod
    def from_gef(cls, path: Union[str, Path]) -> "CPT":
        """Build a CPT from a GEF file."""
        from geoeq.io.gef_reader import read_gef
        gef = read_gef(path)
        d = gef["data"]
        # Common GEF column names: depth, qc, fs, u2 or "Sondeerlengte"...
        depth = _pick(d, ("depth", "z", "sondeerlengte"))
        qc = _pick(d, ("qc", "qc1", "konusw"))
        fs = _pick(d, ("fs", "fs1", "wrijving"), default=None)
        u2 = _pick(d, ("u2", "u", "waterdruk"), default=None)
        return cls(depth=depth, qc=qc, fs=fs, u2=u2,
                    title=Path(path).stem)

    @classmethod
    def from_ags(cls, path: Union[str, Path]) -> "CPT":
        """Build a CPT from an AGS4 file (uses STCN / SCPT / CPTU groups)."""
        from geoeq.io.ags_reader import read_ags
        ags = read_ags(path)
        group = None
        for g in ("SCPT", "STCN", "CPTU", "CPTC"):
            if g in ags:
                group = ags[g]
                break
        if group is None:
            raise ValueError("No CPT-like group (SCPT/STCN/CPTU) in AGS file.")
        depth_key = next((h for h in group["headings"]
                          if "DPTH" in h.upper() or "DEPTH" in h.upper()),
                         group["headings"][0])
        qc_key = next((h for h in group["headings"]
                       if "QC" in h.upper() or "RES" in h.upper()), None)
        fs_key = next((h for h in group["headings"]
                       if "FS" in h.upper()), None)
        u_key = next((h for h in group["headings"]
                      if h.upper() in ("SCPT_U2", "CPTU_U2", "U2")), None)
        rows = group["data"]
        def col(key):
            if key is None:
                return None
            vals = []
            for r in rows:
                try:
                    vals.append(float(r.get(key, "nan")))
                except (TypeError, ValueError):
                    vals.append(np.nan)
            return np.asarray(vals)
        return cls(depth=col(depth_key), qc=col(qc_key),
                    fs=col(fs_key), u2=col(u_key),
                    title=Path(path).stem)

    def __len__(self) -> int:
        return len(self.depth)

    def __repr__(self) -> str:
        return (f"CPT({self.title!r}, n={len(self)}, "
                f"z={self.depth.min():.1f}..{self.depth.max():.1f} m)")


def _pick(d: dict, candidates: Sequence[str], default=None) -> Optional[np.ndarray]:
    keys_lower = {k.lower(): k for k in d.keys()}
    for c in candidates:
        if c.lower() in keys_lower:
            return d[keys_lower[c.lower()]]
    return default
