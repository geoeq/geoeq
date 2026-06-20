"""Tests for ge.io -- CSV, AGS, GEF readers and CPT container."""

import textwrap
from pathlib import Path

import numpy as np
import pytest

from geoeq.io import read_csv, read_ags, read_gef, CPT


# --- CSV reader ---------------------------------------------------

def test_read_csv_basic(tmp_path):
    p = tmp_path / "data.csv"
    p.write_text("depth,qc,fs\n0.5,2.0,0.05\n1.0,3.0,0.06\n")
    res = read_csv(p)
    assert "depth" in res["data"]
    np.testing.assert_allclose(res["data"]["depth"], [0.5, 1.0])
    np.testing.assert_allclose(res["data"]["qc"], [2.0, 3.0])
    assert res["n_rows"] == 2


def test_read_csv_with_units_row(tmp_path):
    p = tmp_path / "u.csv"
    p.write_text("depth,qc\nm,MPa\n0.5,2.0\n1.0,3.0\n")
    res = read_csv(p, units_row=True)
    assert res["units"] == {"depth": "m", "qc": "MPa"}
    np.testing.assert_allclose(res["data"]["qc"], [2.0, 3.0])


def test_read_csv_auto_skip_metadata(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("Project: Mehmet\nSite: Tower 1\ndepth,qc\n0.5,2.0\n")
    res = read_csv(p)
    assert "depth" in res["data"]


# --- AGS reader ---------------------------------------------------

AGS_SAMPLE = textwrap.dedent("""\
    "GROUP","LOCA"
    "HEADING","LOCA_ID","LOCA_NATE","LOCA_NATN"
    "UNIT","","m","m"
    "TYPE","ID","2DP","2DP"
    "DATA","BH1","100.00","200.00"
    "DATA","BH2","150.00","250.00"
    "GROUP","GEOL"
    "HEADING","LOCA_ID","GEOL_TOP","GEOL_DESC"
    "UNIT","","m",""
    "TYPE","X","2DP","X"
    "DATA","BH1","0.00","Topsoil"
    "DATA","BH1","2.00","Sand"
""")


def test_read_ags_two_groups(tmp_path):
    p = tmp_path / "test.ags"
    p.write_text(AGS_SAMPLE)
    groups = read_ags(p)
    assert "LOCA" in groups
    assert "GEOL" in groups
    assert len(groups["LOCA"]["data"]) == 2
    assert groups["LOCA"]["data"][0]["LOCA_ID"] == "BH1"
    assert len(groups["GEOL"]["data"]) == 2


# --- GEF reader ---------------------------------------------------

GEF_SAMPLE = textwrap.dedent("""\
    #GEFID= 1, 1, 0
    #PROCEDURECODE= GEF-CPT-Report
    #COMPANYID= TestCorp
    #COLUMN= 3
    #COLUMNINFO= 1, m, depth, 11
    #COLUMNINFO= 2, MPa, qc, 1
    #COLUMNINFO= 3, MPa, fs, 2
    #COLUMNSEPARATOR= ;
    #COLUMNVOID= 1, -9999
    #EOH=
    0.50; 2.10; 0.05
    1.00; 3.20; 0.08
    1.50; -9999; 0.10
""")


def test_read_gef_basic(tmp_path):
    p = tmp_path / "t.gef"
    p.write_text(GEF_SAMPLE)
    res = read_gef(p)
    assert "depth" in res["data"]
    assert "qc" in res["data"]
    assert "fs" in res["data"]
    np.testing.assert_allclose(res["data"]["depth"], [0.5, 1.0, 1.5])
    np.testing.assert_allclose(res["data"]["qc"][:2], [2.10, 3.20])
    assert np.isnan(res["data"]["qc"][2])  # -9999 -> NaN
    assert res["units"][1] == "MPa"


# --- CPT container ----------------------------------------------

def test_cpt_construct_directly():
    c = CPT(depth=[0.5, 1.0, 1.5], qc=[1, 2, 3], fs=[0.1, 0.2, 0.3])
    assert len(c) == 3
    assert c.title == "CPT"
    np.testing.assert_allclose(c.qc, [1, 2, 3])


def test_cpt_from_gef(tmp_path):
    p = tmp_path / "t.gef"
    p.write_text(GEF_SAMPLE)
    c = CPT.from_gef(p)
    assert len(c) == 3
    assert c.title == "t"
    np.testing.assert_allclose(c.depth, [0.5, 1.0, 1.5])
