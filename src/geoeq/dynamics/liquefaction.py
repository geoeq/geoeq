"""
Simplified-procedure liquefaction triggering analysis.

CSR (Cyclic Stress Ratio) from earthquake loading,
CRR (Cyclic Resistance Ratio) from SPT or CPT correlations,
Factor of safety FS_L = CRR / CSR adjusted by magnitude scaling factor.

Two CRR families are supported:

* **NCEER (Youd et al. 2001)** -- SPT-based; the most widely cited
  simplified procedure.
* **Idriss & Boulanger (2008)** -- updated SPT and CPT relations.

References
----------
* Seed, H. B., Idriss, I. M. (1971). "Simplified procedure for evaluating
  soil liquefaction potential." *J. Soil Mech. Found. Eng.*, 97(SM9).
* Youd, T. L., Idriss, I. M., et al. (2001). "Liquefaction resistance of
  soils." *JGGE*, 127(10), 817-833.
* Idriss, I. M., Boulanger, R. W. (2008). *Soil liquefaction during
  earthquakes*. EERI Monograph 12.
* Idriss, I. M. (1999). "An update to the Seed-Idriss simplified procedure
  for evaluating liquefaction potential."
"""

from __future__ import annotations

import numpy as np

from geoeq.core.validation import check_positive, check_non_negative


# -----------------------------------------------------------------
# 1. Depth-reduction factor rd
# -----------------------------------------------------------------
def depth_reduction(z: float, method: str = "idriss_1999",
                     Mw: float = 7.5) -> float:
    """Stress-reduction factor rd that accounts for the flexibility of the
    soil column under earthquake shear stress.

    Parameters
    ----------
    z : float
        Depth (m).
    method : str
        'liao_whitman_1986' -- simple bilinear (deprecated by NCEER).
        'idriss_1999'       -- magnitude-dependent (NCEER-recommended).
        'cetin_2004'        -- empirical, magnitude-dependent.
    Mw : float
        Moment magnitude (only used by idriss_1999/cetin).

    Returns
    -------
    rd : float

    Reference
    ---------
    Idriss (1999); Youd et al. (2001) Eq. 4.
    """
    check_non_negative(z, "z")
    method = method.lower()
    if method == "liao_whitman_1986":
        if z <= 9.15:
            return float(1 - 0.00765 * z)
        if z <= 23:
            return float(1.174 - 0.0267 * z)
        return float(0.744 - 0.008 * z)
    if method == "idriss_1999":
        # Idriss (1999) magnitude-dependent rd (Youd et al. 2001 Eq. 4)
        alpha = -1.012 - 1.126 * np.sin(z / 11.73 + 5.133)
        beta = 0.106 + 0.118 * np.sin(z / 11.28 + 5.142)
        rd = np.exp(alpha + beta * Mw)
        return float(rd)
    if method == "cetin_2004":
        # Cetin et al. (2004) -- depth- and magnitude-dependent.
        num = (1 + (-23.013 - 2.949 * z + 0.999 * Mw + 0.0525 * z * Mw)
               / (16.258 + 0.201 * np.exp(0.341 * (-z + 0.0785 * Mw + 7.586))))
        den = (1 + (-23.013 + 0.999 * Mw)
               / (16.258 + 0.201 * np.exp(0.341 * (0.0785 * Mw + 7.586))))
        return float(num / den)
    raise ValueError("method must be one of 'liao_whitman_1986', "
                     "'idriss_1999', or 'cetin_2004'.")


# -----------------------------------------------------------------
# 2. Cyclic Stress Ratio (CSR)
# -----------------------------------------------------------------
def liquefaction_csr(
    amax: float, sigma_v: float, sigma_v_eff: float,
    Mw: float = 7.5, z: float = None, rd: float = None,
    g: float = 9.81,
) -> dict:
    """Cyclic stress ratio induced by the earthquake.

        CSR = 0.65 * (amax / g) * (sigma_v / sigma_v_eff) * rd

    Parameters
    ----------
    amax : float
        Peak horizontal ground acceleration as a fraction of g
        (e.g. 0.25 for 0.25g).
    sigma_v : float
        Total vertical stress at depth (kPa).
    sigma_v_eff : float
        Effective vertical stress at depth (kPa).
    Mw : float
        Moment magnitude (used only for rd if rd is None).
    z : float
        Depth (m). Required if rd is None.
    rd : float
        Pre-computed depth-reduction factor. If None, computed from z.

    Returns
    -------
    dict
        ``{'CSR': ..., 'rd': ...}``.

    Reference
    ---------
    Seed & Idriss (1971); Youd et al. (2001) Eq. 1.
    """
    check_positive(sigma_v, "sigma_v")
    check_positive(sigma_v_eff, "sigma_v_eff")
    check_positive(amax, "amax")
    amax_over_g = float(amax)
    if rd is None:
        if z is None:
            raise ValueError("Provide rd or z (and Mw).")
        rd = depth_reduction(z, method="idriss_1999", Mw=Mw)
    CSR = 0.65 * amax_over_g * (sigma_v / sigma_v_eff) * rd
    return {"CSR": float(CSR), "rd": float(rd), "amax_g": float(amax_over_g)}


# -----------------------------------------------------------------
# 3. Cyclic Resistance Ratio (CRR)
# -----------------------------------------------------------------
def liquefaction_crr(
    N160cs: float = None, qc1Ncs: float = None,
    Vs1: float = None, FC: float = 0,
    method: str = "youd_2001", Mw: float = 7.5,
) -> dict:
    """Cyclic resistance ratio at Mw=7.5 from SPT, CPT, or Vs.

    Methods:
    * **youd_2001** -- NCEER SPT (Youd et al. 2001 Eq. 5):
        CRR7.5 = 1/(34 - N160cs) + N160cs/135 + 50/(10*N160cs+45)^2 - 1/200
    * **idriss_boulanger_2008** -- updated SPT (Eq. 4):
        CRR7.5 = exp[ N160cs/14.1 + (N160cs/126)^2 - (N160cs/23.6)^3
                    + (N160cs/25.4)^4 - 2.8 ]
    * **idriss_boulanger_2008_cpt** -- CPT-based (Eq. 6):
        CRR7.5 = exp[ qc1Ncs/540 + (qc1Ncs/67)^2 - (qc1Ncs/80)^3
                    + (qc1Ncs/114)^4 - 3 ]
    * **andrus_stokoe_2000** -- Vs-based (Eq. 9):
        CRR7.5 = 0.022(Vs1/100)^2 + 2.8*(1/(215-Vs1) - 1/215)

    Returns
    -------
    dict with 'CRR', 'method'.

    Reference
    ---------
    Youd et al. (2001); Idriss & Boulanger (2008); Andrus & Stokoe (2000).
    """
    method = method.lower()
    if method in ("youd_2001", "nceer", "youd"):
        if N160cs is None:
            raise ValueError("N160cs required for Youd 2001.")
        if N160cs >= 30:
            CRR = 2.0  # capped (non-liquefiable)
        else:
            CRR = (1 / (34 - N160cs)
                   + N160cs / 135
                   + 50 / (10 * N160cs + 45) ** 2
                   - 1 / 200)
        out = {"CRR": float(CRR), "method": "youd_2001"}
    elif method in ("idriss_boulanger_2008", "ib2008", "ib"):
        if N160cs is None:
            raise ValueError("N160cs required for IB2008 SPT.")
        x = N160cs
        if x > 37.5:  # cap at upper bound
            x = 37.5
        CRR = np.exp(x / 14.1
                     + (x / 126) ** 2
                     - (x / 23.6) ** 3
                     + (x / 25.4) ** 4 - 2.8)
        out = {"CRR": float(CRR), "method": "idriss_boulanger_2008"}
    elif method in ("idriss_boulanger_2008_cpt", "ib2008_cpt", "cpt"):
        if qc1Ncs is None:
            raise ValueError("qc1Ncs required for IB2008 CPT.")
        x = qc1Ncs
        CRR = np.exp(x / 540
                     + (x / 67) ** 2
                     - (x / 80) ** 3
                     + (x / 114) ** 4 - 3.0)
        out = {"CRR": float(CRR), "method": "idriss_boulanger_2008_cpt"}
    elif method in ("andrus_stokoe_2000", "as2000", "vs"):
        if Vs1 is None:
            raise ValueError("Vs1 required for Andrus & Stokoe 2000.")
        if FC <= 5:
            Vs1_star = 215
        elif FC < 35:
            Vs1_star = 215 - 0.5 * (FC - 5)
        else:
            Vs1_star = 200
        if Vs1 >= Vs1_star:
            CRR = 2.0
        else:
            CRR = 0.022 * (Vs1 / 100) ** 2 + 2.8 * (1 / (Vs1_star - Vs1) - 1 / Vs1_star)
        out = {"CRR": float(CRR), "method": "andrus_stokoe_2000",
               "Vs1_star": float(Vs1_star)}
    else:
        raise ValueError(
            "method must be 'youd_2001', 'idriss_boulanger_2008', "
            "'idriss_boulanger_2008_cpt', or 'andrus_stokoe_2000'.")
    return out


# -----------------------------------------------------------------
# 4. Magnitude scaling factor
# -----------------------------------------------------------------
def magnitude_scaling_factor(Mw: float,
                              method: str = "idriss_1999") -> float:
    """Magnitude scaling factor MSF to adjust CRR_7.5 to other magnitudes.

    Methods:
    * **idriss_1999**:  MSF = 6.9 * exp(-Mw/4) - 0.058  (capped at 1.8)
    * **nceer**:        MSF = (10^2.24)/Mw^2.56 (Youd et al. 2001 Eq. 24)
    * **boulanger_idriss_2014**:
        MSF = 1 + (MSFmax - 1)*(8.64 * exp(-Mw/4) - 1.325)

    Reference
    ---------
    Idriss (1999); Youd et al. (2001); Boulanger & Idriss (2014).
    """
    check_positive(Mw, "Mw")
    method = method.lower()
    if method in ("idriss_1999", "idriss"):
        return float(min(1.8, 6.9 * np.exp(-Mw / 4) - 0.058))
    if method in ("nceer", "youd_2001"):
        return float(10 ** 2.24 / Mw ** 2.56)
    if method in ("boulanger_idriss_2014", "bi2014"):
        MSF_max = 1.8  # for sands
        return float(1 + (MSF_max - 1) * (8.64 * np.exp(-Mw / 4) - 1.325))
    raise ValueError(
        "method must be 'idriss_1999', 'nceer', or 'boulanger_idriss_2014'.")


# -----------------------------------------------------------------
# 5. Factor of safety against liquefaction
# -----------------------------------------------------------------
def liquefaction_fos(
    CSR: float, CRR: float, Mw: float = 7.5,
    MSF: float = None, K_sigma: float = 1.0, K_alpha: float = 1.0,
    MSF_method: str = "idriss_1999",
) -> dict:
    """Factor of safety against liquefaction triggering.

        FS_L = (CRR * MSF * K_sigma * K_alpha) / CSR

    Parameters
    ----------
    CSR : float
        Cyclic stress ratio from earthquake loading.
    CRR : float
        Cyclic resistance ratio at Mw=7.5.
    Mw : float
        Moment magnitude (for MSF if not provided).
    MSF : float, optional
        Pre-computed magnitude scaling factor.
    K_sigma : float
        Overburden correction (default 1).
    K_alpha : float
        Static shear stress correction (default 1).

    Returns
    -------
    dict
        ``{'FS': ..., 'MSF': ..., 'liquefies': bool}``.

    Reference
    ---------
    Youd et al. (2001) Eq. 26; Idriss & Boulanger (2008).
    """
    check_positive(CSR, "CSR")
    check_positive(CRR, "CRR")
    if MSF is None:
        MSF = magnitude_scaling_factor(Mw, method=MSF_method)
    FS = CRR * MSF * K_sigma * K_alpha / CSR
    return {
        "FS": float(FS), "MSF": float(MSF),
        "K_sigma": float(K_sigma), "K_alpha": float(K_alpha),
        "liquefies": bool(FS < 1.0),
    }
