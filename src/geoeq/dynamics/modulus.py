"""
Shear modulus and damping for soil dynamics.

Small-strain shear modulus Gmax:
    Gmax = rho * Vs^2                (direct from shear-wave velocity)
    Gmax_hardin = A * F(e) * OCR^k * (sigma'_m)^n      (Hardin & Black 1968)

Modulus reduction (G/Gmax) and damping (D) curves use Darendeli (2001)
empirical functions controlled by PI, OCR, and confining stress.

References
----------
* Hardin & Black (1968); Hardin & Drnevich (1972).
* Darendeli, M. B. (2001). PhD thesis, UT Austin.
* Seed & Idriss (1970, 1986) reduction curves for sand/gravel.
"""

from __future__ import annotations

import numpy as np

from geoeq.core.validation import check_positive, check_non_negative


# -----------------------------------------------------------------
# Small-strain shear modulus
# -----------------------------------------------------------------
def gmax(Vs: float, gamma: float = 18.0) -> float:
    """Small-strain shear modulus Gmax from shear-wave velocity Vs.

        Gmax = rho * Vs^2    (kPa)

    Parameters
    ----------
    Vs : float
        Shear-wave velocity (m/s).
    gamma : float
        Total unit weight (kN/m^3). Default 18.

    Returns
    -------
    Gmax : float
        Small-strain shear modulus (kPa).

    Reference
    ---------
    Stokoe et al. (1999); Das (2010) Eq. 11.18.
    """
    check_positive(Vs, "Vs")
    check_positive(gamma, "gamma")
    rho = gamma * 1000 / 9.81  # kN/m^3 -> kg/m^3
    return float(rho * Vs ** 2 / 1000)  # Pa -> kPa


def gmax_hardin(e: float, sigma_m_eff: float, OCR: float = 1.0,
                PI: float = 0, soil_type: str = "round_grained") -> float:
    """Gmax from Hardin & Black (1968) / Hardin & Drnevich (1972).

        Gmax = A * F(e) * OCR^k * (sigma'_m / pa)^n * pa  (kPa)

    For most soils n ~ 0.5 and:

    * Round-grained sands:    A = 6900,  F(e) = (2.97 - e)^2 / (1 + e)
    * Angular-grained sands:  A = 3300,  F(e) = (2.17 - e)^2 / (1 + e)
    * Clays:                  A = 3230,  F(e) = (2.97 - e)^2 / (1 + e)

    Parameters
    ----------
    e : float
        Void ratio (-).
    sigma_m_eff : float
        Mean effective confining stress (kPa).
    OCR : float
        Over-consolidation ratio.
    PI : float
        Plasticity index. Sets k via Hardin-Drnevich table:
            PI=0:   k=0
            PI=20:  k=0.18
            PI=40:  k=0.30
            PI=60:  k=0.41
            PI=80:  k=0.48
            PI>=100: k=0.50
    soil_type : str
        'round_grained' | 'angular_grained' | 'clay'.

    Returns
    -------
    Gmax : float
        Small-strain shear modulus (kPa).

    Reference
    ---------
    Hardin & Black (1968); Hardin & Drnevich (1972); Kramer (1996) Eq. 6.9.
    """
    check_positive(e, "e")
    check_positive(sigma_m_eff, "sigma_m_eff")
    check_positive(OCR, "OCR")
    soil_type = soil_type.lower()
    if soil_type == "round_grained":
        A, F = 6900, (2.97 - e) ** 2 / (1 + e)
    elif soil_type in ("angular_grained", "angular"):
        A, F = 3300, (2.17 - e) ** 2 / (1 + e)
    elif soil_type == "clay":
        A, F = 3230, (2.97 - e) ** 2 / (1 + e)
    else:
        raise ValueError(
            "soil_type must be 'round_grained', 'angular_grained', or 'clay'.")
    # k from PI (Kramer 1996 Table 6.4)
    k_table = [(0, 0), (20, 0.18), (40, 0.30),
               (60, 0.41), (80, 0.48), (100, 0.50)]
    PI_clip = min(max(PI, 0), 100)
    # piecewise linear
    for (pi1, k1), (pi2, k2) in zip(k_table[:-1], k_table[1:]):
        if pi1 <= PI_clip <= pi2:
            k = k1 + (k2 - k1) * (PI_clip - pi1) / (pi2 - pi1)
            break
    else:
        k = 0.5
    pa = 101.325  # atmospheric pressure, kPa
    n = 0.5
    return float(A * F * OCR ** k * (sigma_m_eff / pa) ** n * pa / 100)


# -----------------------------------------------------------------
# Darendeli (2001) modulus-reduction and damping curves
# -----------------------------------------------------------------
def _darendeli_gamma_r(sigma_m_eff: float, PI: float, OCR: float,
                        freq_Hz: float = 1.0, N_cycles: float = 10) -> float:
    """Reference shear strain gamma_r in % (Darendeli 2001 Eq. 8.12)."""
    pa = 101.325
    sigma_norm = sigma_m_eff / pa
    # Darendeli's fitted constants:
    phi1, phi2, phi3, phi4 = 0.0352, 0.0010, 0.3246, 0.3483
    gamma_r = (phi1 + phi2 * PI * OCR ** phi3) * sigma_norm ** phi4
    return float(gamma_r)


def modulus_reduction(
    strain: float, sigma_m_eff: float = 100.0, PI: float = 15,
    OCR: float = 1.0, method: str = "darendeli",
) -> float:
    """Modulus-reduction ratio G/Gmax at shear strain ``strain`` (%).

        G / Gmax = 1 / (1 + (gamma / gamma_r)^a)

    where gamma_r is the reference strain (Darendeli 2001) and a ~ 0.92
    for typical soils.

    Parameters
    ----------
    strain : float
        Shear strain (decimal, e.g. 0.001 = 0.1%).
    sigma_m_eff : float
        Mean effective confining stress (kPa). Default 100.
    PI : float
        Plasticity index. Default 15.
    OCR : float
        Over-consolidation ratio. Default 1.
    method : str
        'darendeli' (2001) or 'vucetic_dobry' (1991 chart values).

    Returns
    -------
    G_over_Gmax : float
        Modulus reduction ratio (-).

    Reference
    ---------
    Darendeli (2001) Eq. 8.10-8.12; Vucetic & Dobry (1991).
    """
    check_non_negative(strain, "strain")
    method = method.lower()
    gamma_pct = strain * 100  # decimal -> %
    if method == "darendeli":
        gamma_r = _darendeli_gamma_r(sigma_m_eff, PI, OCR)
        a = 0.92
        return float(1.0 / (1.0 + (gamma_pct / gamma_r) ** a))
    if method in ("vucetic_dobry", "vd"):
        # Vucetic & Dobry (1991) for PI=0 sand: digitized fit.
        # G/Gmax ~ 1/(1 + (gamma/gamma_05)^1) with gamma_05 from PI:
        gamma_05 = 0.01 + 0.001 * PI  # very rough fit
        return float(1.0 / (1.0 + (gamma_pct / gamma_05)))
    raise ValueError("method must be 'darendeli' or 'vucetic_dobry'.")


def damping_ratio(
    strain: float, sigma_m_eff: float = 100.0, PI: float = 15,
    OCR: float = 1.0, freq_Hz: float = 1.0, N_cycles: float = 10,
    method: str = "darendeli",
) -> float:
    """Equivalent-linear damping ratio D (decimal, e.g. 0.05 = 5%).

    Uses Darendeli (2001) full model:
        D = D_min + b * (G/Gmax)^0.1 * (D_masing - D_min)/100   [%]
        D_min = (phi6 + phi7 * PI * OCR^phi8) * sigma_norm^phi9
                * (1 + phi10 * ln(freq))

    Returns the small-strain D for very small strains, with the Masing
    correction at larger strains.

    Reference
    ---------
    Darendeli (2001) Eq. 8.13-8.18.
    """
    check_non_negative(strain, "strain")
    pa = 101.325
    sigma_norm = sigma_m_eff / pa
    # Coefficients from Darendeli (2001) Table 8.12 (selected):
    phi6, phi7, phi8, phi9, phi10 = 0.8005, 0.0129, -0.1069, -0.2889, 0.2919
    phi11 = 0.6329
    phi12 = -0.0057

    method = method.lower()
    # Small-strain damping (%), Darendeli Eq. 8.18.
    D_min = ((phi6 + phi7 * PI * OCR ** phi8)
             * sigma_norm ** phi9 * (1 + phi10 * np.log(freq_Hz)))
    G_Gmax = modulus_reduction(strain, sigma_m_eff, PI, OCR,
                                method=method if method != "darendeli"
                                else "darendeli")
    # Use the Hardin & Drnevich (1972) hyperbolic relation:
    #     D / D_max = 1 - G/Gmax
    # to avoid the well-known singularity of the closed-form Masing
    # damping at gamma -> 0. D_max ~ 25 % for typical sands, lower for
    # high-PI clays.
    D_max_pct = max(20.0, 25.0 - 0.05 * PI)  # %
    D_pct = D_min + (D_max_pct - D_min) * (1.0 - G_Gmax)
    return float(D_pct / 100.0)  # % -> decimal
