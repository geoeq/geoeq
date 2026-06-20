"""Tests for ge.dynamics."""

import numpy as np
import pytest

from geoeq.dynamics import (
    gmax, gmax_hardin, modulus_reduction, damping_ratio,
    depth_reduction, liquefaction_csr, liquefaction_crr,
    liquefaction_fos, magnitude_scaling_factor, response_spectrum,
)


# --- Gmax -------------------------------------------------------

def test_gmax_from_Vs():
    # gamma=20 kN/m^3, Vs=200 m/s -> rho = 2038 kg/m^3 -> G = 2038*40000 = 81 MPa
    G = gmax(Vs=200, gamma=20)
    assert G == pytest.approx(81550, rel=0.02)


def test_gmax_hardin_round_sand():
    G = gmax_hardin(e=0.6, sigma_m_eff=100, soil_type="round_grained")
    assert G > 0


def test_gmax_hardin_OCR_increases_G():
    G1 = gmax_hardin(e=0.6, sigma_m_eff=100, OCR=1, PI=40)
    G4 = gmax_hardin(e=0.6, sigma_m_eff=100, OCR=4, PI=40)
    assert G4 > G1


def test_gmax_hardin_bad_soil_type():
    with pytest.raises(ValueError):
        gmax_hardin(0.6, 100, soil_type="weird")


# --- modulus reduction -----------------------------------------

def test_modulus_reduction_at_zero_strain():
    # Very small strain -> ~1.0
    assert modulus_reduction(strain=1e-7, PI=15) > 0.99


def test_modulus_reduction_decreases_with_strain():
    G1 = modulus_reduction(strain=1e-5, PI=15)
    G2 = modulus_reduction(strain=1e-3, PI=15)
    G3 = modulus_reduction(strain=1e-2, PI=15)
    assert G1 > G2 > G3


def test_modulus_reduction_PI_effect():
    # Higher PI -> less reduction at a given strain.
    G_low = modulus_reduction(strain=1e-3, PI=0)
    G_high = modulus_reduction(strain=1e-3, PI=100)
    assert G_high > G_low


# --- damping ---------------------------------------------------

def test_damping_at_small_strain_is_small():
    D = damping_ratio(strain=1e-7, PI=15)
    assert 0 < D < 0.05  # < 5%


def test_damping_increases_with_strain():
    D1 = damping_ratio(strain=1e-5, PI=15)
    D2 = damping_ratio(strain=1e-3, PI=15)
    assert D2 > D1


# --- depth reduction rd ----------------------------------------

def test_rd_at_surface():
    # rd ~ 1 at the surface.
    rd = depth_reduction(z=0, method="idriss_1999")
    assert rd == pytest.approx(1.0, abs=0.05)


def test_rd_decreases_with_depth():
    rd_shallow = depth_reduction(z=2, method="idriss_1999")
    rd_deep = depth_reduction(z=15, method="idriss_1999")
    assert rd_deep < rd_shallow


def test_rd_liao_whitman_bilinear():
    # At z=5 m: rd = 1 - 0.00765*5 = 0.962
    rd = depth_reduction(z=5, method="liao_whitman_1986")
    assert rd == pytest.approx(0.962, abs=0.01)


def test_rd_invalid_method():
    with pytest.raises(ValueError):
        depth_reduction(z=5, method="bogus")


# --- CSR -------------------------------------------------------

def test_csr_classic_example():
    # Seed & Idriss example: amax=0.2g, sigma_v=100, sigma_v_eff=50, z=5
    # rd at z=5, M=7.5 ~ 0.97; CSR = 0.65 * 0.2 * (100/50) * 0.97 ~ 0.252
    res = liquefaction_csr(amax=0.2, sigma_v=100, sigma_v_eff=50,
                            z=5, Mw=7.5)
    assert 0.2 < res["CSR"] < 0.3


def test_csr_amax_required_as_g_fraction():
    # API contract: amax must be passed as fraction of g.
    res = liquefaction_csr(amax=0.2, sigma_v=100, sigma_v_eff=50, z=5)
    assert res["amax_g"] == pytest.approx(0.2)


def test_csr_provide_rd_directly():
    res = liquefaction_csr(amax=0.2, sigma_v=100, sigma_v_eff=50, rd=0.95)
    expected = 0.65 * 0.2 * 2 * 0.95
    assert res["CSR"] == pytest.approx(expected)


# --- CRR -------------------------------------------------------

def test_crr_youd_low_N():
    # N160cs = 10 -> CRR ~ 0.11 (textbook value).
    res = liquefaction_crr(N160cs=10, method="youd_2001")
    assert 0.10 < res["CRR"] < 0.13


def test_crr_youd_high_N_safe():
    # N160cs >= 30 -> CRR capped (non-liquefiable).
    res = liquefaction_crr(N160cs=30, method="youd_2001")
    assert res["CRR"] >= 0.5


def test_crr_idriss_boulanger():
    res = liquefaction_crr(N160cs=15, method="idriss_boulanger_2008")
    assert 0.05 < res["CRR"] < 0.25


def test_crr_cpt():
    res = liquefaction_crr(qc1Ncs=100, method="idriss_boulanger_2008_cpt")
    assert res["CRR"] > 0


def test_crr_vs():
    res = liquefaction_crr(Vs1=150, method="andrus_stokoe_2000")
    assert res["CRR"] > 0


def test_crr_invalid_method():
    with pytest.raises(ValueError):
        liquefaction_crr(N160cs=15, method="bogus")


# --- MSF -------------------------------------------------------

def test_msf_at_7_5_close_to_1():
    msf = magnitude_scaling_factor(7.5, method="idriss_1999")
    assert msf == pytest.approx(1.0, abs=0.1)


def test_msf_at_low_magnitude_larger():
    msf_small = magnitude_scaling_factor(5.5)
    msf_large = magnitude_scaling_factor(8.5)
    assert msf_small > msf_large


def test_msf_nceer():
    # Youd et al. (2001) Table for M=7.5 should be ~1.0
    msf = magnitude_scaling_factor(7.5, method="nceer")
    assert msf == pytest.approx(1.0, abs=0.05)


# --- FoS -------------------------------------------------------

def test_fos_safe():
    res = liquefaction_fos(CSR=0.10, CRR=0.30, Mw=7.5)
    assert res["FS"] > 1
    assert not res["liquefies"]


def test_fos_liquefies():
    res = liquefaction_fos(CSR=0.30, CRR=0.10, Mw=7.5)
    assert res["FS"] < 1
    assert res["liquefies"]


def test_fos_includes_MSF():
    # At small magnitudes MSF > 1 increases FS.
    fs1 = liquefaction_fos(CSR=0.15, CRR=0.15, Mw=7.5)["FS"]
    fs2 = liquefaction_fos(CSR=0.15, CRR=0.15, Mw=5.5)["FS"]
    assert fs2 > fs1


# --- response spectrum ----------------------------------------

def test_response_spectrum_returns_figure():
    T = np.linspace(0.01, 5, 100)
    Sa = 0.5 * np.exp(-T)
    fig = response_spectrum(T, Sa, site_class="D")
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close(fig)
