"""
Seepage and groundwater flow.

References
----------
* Darcy, H. (1856). *Les fontaines publiques de la ville de Dijon*.
* Das (2010), Ch. 5.
* Budhu (2011), Ch. 9.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from geoeq.core.constants import GAMMA_WATER
from geoeq.core.validation import check_positive


def darcy_flow(k: float, i: float, A: float) -> float:
    """Darcy's law: Q = k * i * A (m^3/s).

    Parameters
    ----------
    k : float
        Hydraulic conductivity (m/s).
    i : float
        Hydraulic gradient (-, dimensionless).
    A : float
        Cross-sectional area (m^2).

    Reference
    ---------
    Darcy (1856); Das (2010) Eq. 5.11.
    """
    check_positive(k, "k")
    check_positive(A, "A")
    return k * i * A


def hydraulic_gradient(dh: float, L: float) -> float:
    """Hydraulic gradient i = delta h / L.

    Parameters
    ----------
    dh : float
        Head loss (m).
    L : float
        Flow-path length (m).

    Reference
    ---------
    Das (2010) Eq. 5.10.
    """
    check_positive(L, "L")
    return dh / L


def critical_gradient(Gs: float, e: float) -> float:
    """Critical (boiling) hydraulic gradient i_cr = (Gs - 1)/(1 + e).

    When i >= i_cr in upward flow, effective stress vanishes and quicksand
    (heave) develops.

    Reference
    ---------
    Das (2010) Eq. 5.18; Terzaghi & Peck (1948).
    """
    check_positive(Gs, "Gs")
    check_positive(e, "e")
    if Gs <= 1:
        raise ValueError("Gs must be > 1.")
    return (Gs - 1) / (1 + e)


def equivalent_k(
    k_layers: Sequence[float],
    H_layers: Sequence[float],
    direction: str = "horizontal",
) -> float:
    """Equivalent permeability of a layered system.

    horizontal: k_eq = sum(k_i * H_i) / sum(H_i)        (parallel flow)
    vertical:   k_eq = sum(H_i) / sum(H_i / k_i)        (series flow)

    Reference
    ---------
    Das (2010) Eq. 5.14 (horizontal), 5.15 (vertical).
    """
    k = np.asarray(k_layers, dtype=float)
    H = np.asarray(H_layers, dtype=float)
    if k.shape != H.shape:
        raise ValueError("k_layers and H_layers must have the same length.")
    check_positive(k, "k_layers")
    check_positive(H, "H_layers")
    direction = direction.lower()
    if direction in ("horizontal", "h", "parallel"):
        return float(np.sum(k * H) / np.sum(H))
    elif direction in ("vertical", "v", "series"):
        return float(np.sum(H) / np.sum(H / k))
    raise ValueError("direction must be 'horizontal' or 'vertical'.")


def flow_net(Nf: float, Nd: float, k: float, dh: float, L: float = 1.0) -> float:
    """Seepage discharge from a flow net.

    Q = k * dh * (Nf / Nd) * L   (m^3/s per metre length normal to plane)

    Parameters
    ----------
    Nf : float
        Number of flow channels.
    Nd : float
        Number of equipotential drops.
    k : float
        Hydraulic conductivity (m/s).
    dh : float
        Total head loss across the system (m).
    L : float
        Length perpendicular to the flow plane (m). Default 1 m.

    Reference
    ---------
    Das (2010) Eq. 5.20.
    """
    check_positive(Nf, "Nf")
    check_positive(Nd, "Nd")
    check_positive(k, "k")
    return k * dh * (Nf / Nd) * L
