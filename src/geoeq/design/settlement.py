"""
Foundation settlement.

Three components, each a separate function:

* Immediate (elastic) settlement -- Janbu / Steinbrenner.
* Primary consolidation -- 1-D Terzaghi for NC or OC soils.
* Secondary compression -- Cα log(t2/t1).

Plus Schmertmann's strain-influence method for sandy soils and the
time-dependent functions (Tv <-> U).

References
----------
* Terzaghi, K. (1925). *Erdbaumechanik auf bodenphysikalischer Grundlage*.
* Janbu, N., Bjerrum, L., Kjaernsli, B. (1956). "Veiledning ved losning
  av fundamenteringsoppgaver." NGI Publication 16.
* Schmertmann, J. H. (1970). "Static cone to compute static settlement
  over sand." *J. Soil Mech. Found. Eng.*, ASCE, 96(SM3), 1011-1043.
* Schmertmann, J. H., Hartman, J. P., Brown, P. R. (1978). "Improved
  strain influence factor diagrams." *J. Geotech. Eng.*, 104(GT8), 1131-1135.
* Mesri, G., Godlewski, P. M. (1977). "Time- and stress-compressibility
  interrelationship." *J. Geotech. Eng. Div.*, ASCE, 103(GT5), 417-430.
* Das (2014), Ch. 5.
"""

from __future__ import annotations

from typing import Sequence, Union

import numpy as np

from geoeq.core.validation import check_positive, check_non_negative


# -----------------------------------------------------------------
# 1. Immediate / elastic settlement
# -----------------------------------------------------------------
def settlement_immediate(
    q: float, B: float, Es: float, mu: float = 0.3,
    Ip: float = None, shape: str = "rigid_square",
) -> float:
    """Elastic (immediate) settlement under a uniformly loaded footing.

        S_i = q * B * (1 - mu^2) / Es * Ip

    Parameters
    ----------
    q : float
        Surface pressure (kPa).
    B : float
        Footing width (m) -- shorter dimension for rectangles.
    Es : float
        Soil elastic modulus (kPa).
    mu : float
        Poisson's ratio (-). Default 0.3.
    Ip : float, optional
        Influence factor. If None, picked from ``shape``.
    shape : str
        'rigid_square'  -> Ip = 0.88
        'rigid_circle'  -> Ip = 0.79
        'flexible_centre' -> Ip = 1.12
        'flexible_corner' -> Ip = 0.56
        'rigid_strip'   -> Ip = 2.0 (approx, depth-dependent)

    Reference
    ---------
    Janbu et al. (1956); Das (2014) Eq. 5.20, Table 5.1.
    """
    check_positive(B, "B")
    check_positive(Es, "Es")
    if Ip is None:
        table = {
            "rigid_square": 0.88, "rigid_circle": 0.79,
            "flexible_centre": 1.12, "flexible_center": 1.12,
            "flexible_corner": 0.56, "rigid_strip": 2.0,
        }
        if shape not in table:
            raise ValueError(
                f"shape must be one of {list(table.keys())}, got {shape!r}.")
        Ip = table[shape]
    return float(q * B * (1 - mu ** 2) / Es * Ip)


# -----------------------------------------------------------------
# 2. Primary consolidation (1-D Terzaghi)
# -----------------------------------------------------------------
def settlement_primary(
    Cc: float, e0: float, H: float, sigma0: float, delta_sigma: float,
    Cr: float = None, sigma_pc: float = None, kind: str = "auto",
) -> float:
    """1-D primary consolidation settlement.

    Three regimes:

    1. **NC** (normally consolidated):
       S = (Cc * H / (1 + e0)) * log10((sigma0 + delta_sigma) / sigma0)

    2. **OC**, final stress below preconsolidation:
       S = (Cr * H / (1 + e0)) * log10((sigma0 + delta_sigma) / sigma0)

    3. **OC**, final stress crosses preconsolidation:
       S = (Cr*H/(1+e0)) * log10(sigma_pc / sigma0)
         + (Cc*H/(1+e0)) * log10((sigma0+delta_sigma)/sigma_pc)

    Parameters
    ----------
    Cc : float
        Compression index (-).
    e0 : float
        Initial void ratio (-).
    H : float
        Drainage path length (m) -- full layer thickness for one-way drainage
        or half-thickness for two-way drainage; pass actual layer thickness.
    sigma0 : float
        Initial effective vertical stress (kPa) at mid-depth.
    delta_sigma : float
        Stress increase at mid-depth (kPa).
    Cr : float
        Recompression (swell) index (-). Required for OC cases.
    sigma_pc : float
        Preconsolidation pressure (kPa). Required for OC cases.
    kind : str
        'NC', 'OC', or 'auto' (decide from sigma_pc).

    Returns
    -------
    S : float
        Primary settlement (m).

    Reference
    ---------
    Terzaghi (1925); Das (2014) Eq. 5.31-5.34.
    """
    check_positive(H, "H")
    check_positive(sigma0, "sigma0")
    check_positive(delta_sigma, "delta_sigma")
    if e0 <= 0:
        raise ValueError("e0 must be positive.")

    sigma_f = sigma0 + delta_sigma

    if kind.lower() == "auto":
        if sigma_pc is None or sigma_pc <= sigma0:
            kind = "NC"
        elif sigma_f <= sigma_pc:
            kind = "OC_below"
        else:
            kind = "OC_cross"
    else:
        kind = kind.upper()

    if kind == "NC":
        return float(Cc * H / (1 + e0) * np.log10(sigma_f / sigma0))

    if Cr is None:
        raise ValueError("Cr is required for an OC analysis.")
    if sigma_pc is None:
        raise ValueError("sigma_pc is required for an OC analysis.")

    if kind in ("OC_below", "OC"):
        return float(Cr * H / (1 + e0) * np.log10(sigma_f / sigma0))

    if kind == "OC_cross":
        S1 = Cr * H / (1 + e0) * np.log10(sigma_pc / sigma0)
        S2 = Cc * H / (1 + e0) * np.log10(sigma_f / sigma_pc)
        return float(S1 + S2)

    raise ValueError(
        f"kind must be 'NC', 'OC', 'OC_below', 'OC_cross', or 'auto', "
        f"got {kind!r}.")


# -----------------------------------------------------------------
# 3. Secondary compression (creep)
# -----------------------------------------------------------------
def settlement_secondary(
    C_alpha: float, H: float, t1: float, t2: float, e_p: float = None,
) -> float:
    """Secondary compression S_s = C_alpha_eps * H * log10(t2/t1).

    Parameters
    ----------
    C_alpha : float
        Secondary compression index. If ``e_p`` is given this is the
        e-vs-log-t form C_alpha (slope of e vs log t after EOP);
        converted to strain form via C_alpha_eps = C_alpha / (1 + e_p).
        Otherwise treated directly as the strain form.
    H : float
        Layer thickness (m).
    t1, t2 : float
        Reference (end of primary) and target times (any consistent unit).
    e_p : float, optional
        Void ratio at end-of-primary, for converting C_alpha -> C_alpha_eps.

    Reference
    ---------
    Mesri & Godlewski (1977); Das (2014) Eq. 5.51.
    """
    check_positive(H, "H")
    check_positive(t1, "t1")
    check_positive(t2, "t2")
    if t2 <= t1:
        raise ValueError("t2 must be > t1.")
    if e_p is not None:
        C_alpha_eps = C_alpha / (1 + e_p)
    else:
        C_alpha_eps = C_alpha
    return float(C_alpha_eps * H * np.log10(t2 / t1))


# -----------------------------------------------------------------
# 4. Schmertmann (1970/1978) strain-influence settlement (granular soils)
# -----------------------------------------------------------------
def settlement_schmertmann(
    q_net: float, B: float, Es: Union[float, Sequence[float]],
    layers: Sequence[tuple] = None, Df: float = 0.0,
    sigma_v_at_base: float = None, gamma: float = 18.0,
    t_years: float = 1.0, shape: str = "square",
) -> dict:
    """Schmertmann's strain-influence settlement for sandy soils.

        S = C1 * C2 * q_net * sum_i ( Iz_i / Es_i * dz_i )

    where:
        C1 = 1 - 0.5 * (sigma_v' / q_net)            -- embedment correction
        C2 = 1 + 0.2 * log10(t_years / 0.1)          -- creep correction

    The strain-influence diagram Iz is taken from Schmertmann et al. (1978):

    * Square / circular (axisymmetric):  Iz peaks at 0.5*B (Iz_peak),
      = 0.1 at z=0, 0 at z=2B.
    * Strip (plane strain):             Iz peaks at 1.0*B,
      = 0.2 at z=0, 0 at z=4B.

    Iz_peak = 0.5 + 0.1 * sqrt(q_net / sigma_v_at_peak).

    Parameters
    ----------
    q_net : float
        Net surface pressure (kPa) above sigma_v at footing base.
    B : float
        Footing width (m).
    Es : float or sequence
        Elastic modulus of subsoil (kPa). If a single number, applied
        to one layer of depth 0 .. (2B or 4B). If sequence, must match ``layers``.
    layers : sequence of (top, bot) in metres below the footing.
        Optional. If None and ``Es`` is scalar, the function builds a
        single-layer model spanning the full influence depth.
    Df : float
        Footing depth (m).
    sigma_v_at_base : float, optional
        Effective vertical stress at the base of the footing (kPa).
        If None, estimated as ``gamma * Df``.
    gamma : float
        Soil unit weight, only used if sigma_v_at_base is None.
    t_years : float
        Time after construction for the creep correction (yrs). Default 1.
    shape : str
        'square' / 'circle' (axisymmetric) or 'strip' (plane strain).

    Returns
    -------
    dict
        {'S': settlement_metres, 'C1', 'C2', 'Iz_peak', 'layers': [...]}

    Reference
    ---------
    Schmertmann (1970); Schmertmann et al. (1978); Das (2014) Ch. 5.10.
    """
    check_positive(q_net, "q_net")
    check_positive(B, "B")
    shape = shape.lower()
    if shape in ("square", "circle", "circular", "axisymmetric"):
        z_peak = 0.5 * B
        z_max = 2.0 * B
        Iz_z0 = 0.1
    elif shape == "strip":
        z_peak = 1.0 * B
        z_max = 4.0 * B
        Iz_z0 = 0.2
    else:
        raise ValueError("shape must be 'square', 'circle', or 'strip'.")

    if sigma_v_at_base is None:
        sigma_v_at_base = gamma * Df

    # sigma_v at z_peak (initial effective stress).
    sigma_v_peak = sigma_v_at_base + gamma * z_peak
    Iz_peak = 0.5 + 0.1 * np.sqrt(q_net / max(sigma_v_peak, 1e-6))

    # Build the Iz(z) function: linear from (0, Iz_z0) to (z_peak, Iz_peak)
    # then linear to (z_max, 0).
    def Iz(z):
        if z < 0 or z > z_max:
            return 0.0
        if z <= z_peak:
            return Iz_z0 + (Iz_peak - Iz_z0) * (z / z_peak)
        return Iz_peak * (1 - (z - z_peak) / (z_max - z_peak))

    # Build layers if not provided.
    if layers is None:
        if np.isscalar(Es):
            layers = [(0.0, z_max)]
            Es = [float(Es)]
        else:
            raise ValueError(
                "If Es is a sequence, layers must be provided.")
    Es_arr = np.atleast_1d(np.asarray(Es, dtype=float))
    if len(Es_arr) != len(layers):
        raise ValueError("Es and layers must have the same length.")

    # Integrate sum(Iz / Es * dz) using midpoint of each layer.
    integ = 0.0
    layer_info = []
    for (top, bot), Esi in zip(layers, Es_arr):
        dz = bot - top
        z_mid = 0.5 * (top + bot)
        contribution = Iz(z_mid) / Esi * dz
        integ += contribution
        layer_info.append({"top": top, "bot": bot, "Es": Esi,
                           "Iz_mid": Iz(z_mid), "dS": contribution})

    C1 = max(0.5, 1 - 0.5 * sigma_v_at_base / q_net)
    C2 = max(1.0, 1 + 0.2 * np.log10(max(t_years, 0.1) / 0.1))
    S = C1 * C2 * q_net * integ

    return {
        "S": float(S), "C1": float(C1), "C2": float(C2),
        "Iz_peak": float(Iz_peak), "z_peak": float(z_peak),
        "z_max": float(z_max), "layers": layer_info,
    }


# -----------------------------------------------------------------
# 5. Time-rate (Terzaghi 1-D consolidation theory)
# -----------------------------------------------------------------
def time_factor(U: float) -> float:
    """Time factor Tv from average degree of consolidation U (0..1).

    Approximate Terzaghi solution (Das Eq. 5.46):
        Tv = (pi/4) * U^2                   for U < 0.6
        Tv = 1.781 - 0.933 * log10(100*(1-U))   for U >= 0.6

    Reference
    ---------
    Terzaghi (1925); Das (2014) Eq. 5.46.
    """
    if not 0 <= U < 1:
        raise ValueError("U must be in [0, 1).")
    if U < 0.6:
        return (np.pi / 4) * U ** 2
    return 1.781 - 0.933 * np.log10(100 * (1 - U))


def consolidation_degree(Tv: float) -> float:
    """Average degree of consolidation U (0..1) from time factor Tv.

    Inverse of ``time_factor`` (Das Eq. 5.46).
    """
    if Tv < 0:
        raise ValueError("Tv must be non-negative.")
    if Tv <= 0.286:  # boundary at U = 0.6
        return float(np.sqrt(4 * Tv / np.pi))
    return float(1 - 10 ** ((1.781 - Tv) / 0.933) / 100)


def consolidation_time(U: float, Hdr: float, cv: float) -> float:
    """Time to reach degree of consolidation U.

        t = Tv * Hdr^2 / cv

    Parameters
    ----------
    U : float
        Target degree of consolidation (0..1).
    Hdr : float
        Drainage path length (m). Half the layer thickness for double drainage.
    cv : float
        Coefficient of consolidation (m^2/s, or any consistent unit).

    Returns
    -------
    t : float
        Time (same time units as cv).

    Reference
    ---------
    Das (2014) Eq. 5.43.
    """
    Tv = time_factor(U)
    check_positive(Hdr, "Hdr")
    check_positive(cv, "cv")
    return float(Tv * Hdr ** 2 / cv)
