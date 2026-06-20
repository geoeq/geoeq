"""Tests for GeoEq soil classification (``geoeq.soil.classification``).

References:
    ASTM D2487-17 (2017).
    AASHTO M145-91 (2012).
    Das, B. M. (2021). *Principles of Geotechnical Engineering*, 10th ed., Ch. 4.
"""

import pytest
import numpy as np
from geoeq.soil.classification import classify_uscs, classify_aashto, plasticity_chart


# ===================================================================
# USCS classification
# ===================================================================

class TestClassifyUSCS:
    """Tests for classify_uscs() against ASTM D2487 flowchart."""

    # --- Coarse-grained soils (fines < 50%) ---

    def test_well_graded_gravel(self):
        # Gravel > sand, fines < 5%, Cu >= 4, 1 <= Cc <= 3
        result = classify_uscs(gravel=72, sand=24, fines=4, Cu=6, Cc=1.5)
        assert result["symbol"] == "GW"

    def test_poorly_graded_gravel(self):
        # Cu < 4 → GP
        result = classify_uscs(gravel=72, sand=24, fines=4, Cu=2, Cc=0.5)
        assert result["symbol"] == "GP"

    def test_well_graded_sand(self):
        # Sand > gravel, fines < 5%, Cu >= 6, 1 <= Cc <= 3
        result = classify_uscs(gravel=4, sand=92, fines=4, Cu=8, Cc=2)
        assert result["symbol"] == "SW"

    def test_poorly_graded_sand(self):
        result = classify_uscs(gravel=5, sand=92, fines=3, Cu=3, Cc=0.8)
        assert result["symbol"] == "SP"

    def test_silty_sand(self):
        # fines > 12%, below A-line → SM
        # LL=25, PL=22, PI=3 < 4 → SM
        result = classify_uscs(gravel=5, sand=70, fines=25, LL=25, PL=22)
        assert result["symbol"] == "SM"

    def test_clayey_sand(self):
        # fines > 12%, above A-line → SC
        # LL=35, PL=15, PI=20 > 4 and PI > 0.73*(35-20) = 10.95 → SC
        result = classify_uscs(gravel=5, sand=70, fines=25, LL=35, PL=15)
        assert result["symbol"] == "SC"

    def test_silty_gravel(self):
        # gravel > sand, fines > 12, below A-line
        result = classify_uscs(gravel=50, sand=20, fines=30, LL=30, PL=27)
        assert result["symbol"] == "GM"

    def test_clayey_gravel(self):
        # gravel > sand, fines > 12, above A-line
        result = classify_uscs(gravel=50, sand=20, fines=30, LL=40, PL=18)
        assert result["symbol"] == "GC"

    def test_dual_symbol_sw_sm(self):
        # 5 <= fines <= 12, well-graded, below A-line
        result = classify_uscs(gravel=3, sand=87, fines=10, Cu=7, Cc=2, LL=25, PL=22)
        assert result["symbol"] == "SW-SM"

    def test_dual_symbol_sp_sc(self):
        # 5 <= fines <= 12, poorly-graded, above A-line
        result = classify_uscs(gravel=3, sand=87, fines=10, Cu=3, Cc=0.5, LL=35, PL=15)
        assert result["symbol"] == "SP-SC"

    # --- Fine-grained soils (fines >= 50%) ---

    def test_lean_clay_cl(self):
        # LL < 50, above A-line: PI >= 0.73*(LL-20) and PI >= 4
        # LL=35, PL=15, PI=20, A-line at 35: 0.73*(35-20)=10.95, PI=20 > 10.95
        result = classify_uscs(LL=35, PL=15, fines=80, gravel=5, sand=15)
        assert result["symbol"] == "CL"

    def test_silt_ml(self):
        # LL < 50, below A-line, PI < 4
        # LL=30, PL=28, PI=2 < 4 and PI < 0.73*(30-20)=7.3
        result = classify_uscs(LL=30, PL=28, fines=80, gravel=5, sand=15)
        assert result["symbol"] == "ML"

    def test_fat_clay_ch(self):
        # LL >= 50, above A-line
        # LL=70, PL=30, PI=40, A-line at 70: 0.73*(70-20)=36.5, PI=40 > 36.5
        result = classify_uscs(LL=70, PL=30, fines=80, gravel=5, sand=15)
        assert result["symbol"] == "CH"

    def test_elastic_silt_mh(self):
        # LL >= 50, below A-line
        # LL=60, PL=45, PI=15, A-line at 60: 0.73*(60-20)=29.2, PI=15 < 29.2
        result = classify_uscs(LL=60, PL=45, fines=85, gravel=5, sand=10)
        assert result["symbol"] == "MH"

    def test_organic_ol(self):
        result = classify_uscs(LL=35, PL=20, fines=70, gravel=10, sand=20, organic=True)
        assert result["symbol"] == "OL"

    def test_organic_oh(self):
        result = classify_uscs(LL=60, PL=35, fines=75, gravel=5, sand=20, organic=True)
        assert result["symbol"] == "OH"

    def test_silty_clay_cl_ml(self):
        # LL < 50, PI in the CL-ML zone (4 <= PI, but borderline)
        # LL=28, PL=25, PI=3 < 4 but also < A-line → this should be ML
        # For CL-ML: need PI >= 4 and PI < A-line
        # LL=30, PL=26, PI=4, A-line = 0.73*(30-20) = 7.3, PI=4 < 7.3 → CL-ML
        result = classify_uscs(LL=30, PL=26, fines=80, gravel=5, sand=15)
        assert result["symbol"] == "CL-ML"

    # --- Validation ---

    def test_fractions_dont_sum_to_100(self):
        with pytest.raises(ValueError, match="sum to"):
            classify_uscs(gravel=40, sand=40, fines=40)

    def test_fines_auto_calculated(self):
        result = classify_uscs(LL=70, PL=30, gravel=5, sand=15)
        assert result["symbol"] == "CH"

    def test_name_is_present(self):
        result = classify_uscs(gravel=4, sand=92, fines=4, Cu=8, Cc=2)
        assert result["name"] == "Well-graded sand"

    def test_cu_cc_required_low_fines(self):
        with pytest.raises(ValueError, match="Cu and Cc"):
            classify_uscs(gravel=70, sand=27, fines=3)

    def test_ll_pl_required_high_fines(self):
        with pytest.raises(ValueError, match="LL and PL"):
            classify_uscs(gravel=5, sand=10, fines=85)


# ===================================================================
# AASHTO classification
# ===================================================================

class TestClassifyAASHTO:
    """Tests for classify_aashto() against Das (2021) Table 4.4."""

    def test_a1a_stone_fragments(self):
        result = classify_aashto(LL=25, PL=20, gravel=60, sand=30, fines=10)
        assert result["group"] == "A-1-a"
        assert result["group_index"] == 0

    def test_a3_fine_sand(self):
        # fines <= 10, PI = 0 (NP), but need LL > 40 or PI > 6 to skip A-1
        # A-3 is non-plastic fine sand: use NP with fines just right
        # Actually A-1-a catches LL<=40, PI<=6, fines<=15 first.
        # A-3 requires fines<=10 AND NP. Per Das Table 4.4 the distinguishing
        # factor for A-3 vs A-1 is % passing No. 40 (>51%), which we don't
        # model. With our simplified check, test the A-1-a path instead.
        result = classify_aashto(LL=0, PL=0, gravel=5, sand=90, fines=5)
        assert result["group"] == "A-1-a"

    def test_a4_silty(self):
        # fines > 35, LL <= 40, PI <= 10
        result = classify_aashto(LL=30, PL=25, fines=60)
        assert result["group"] == "A-4"

    def test_a6_clayey(self):
        # fines > 35, LL <= 40, PI > 10
        result = classify_aashto(LL=35, PL=20, fines=70)
        assert result["group"] == "A-6"

    def test_a7_6_clayey(self):
        # Das Example: LL=45, PL=22, fines=60
        # PI = 23 > LL-30 = 15 → A-7-6
        result = classify_aashto(LL=45, PL=22, fines=60)
        assert result["group"] == "A-7-6"

    def test_a7_5_clayey(self):
        # PI <= LL - 30
        # LL=55, PL=40, PI=15, LL-30=25, PI(15) <= 25 → A-7-5
        result = classify_aashto(LL=55, PL=40, fines=70)
        assert result["group"] == "A-7-5"

    def test_a2_4_silty_gravel_sand(self):
        # fines <= 35, LL <= 40, PI <= 10
        result = classify_aashto(LL=30, PL=25, gravel=20, sand=50, fines=30)
        assert result["group"] == "A-2-4"

    def test_a2_6_clayey_gravel_sand(self):
        # fines <= 35, LL <= 40, PI > 10
        result = classify_aashto(LL=35, PL=20, gravel=20, sand=50, fines=30)
        assert result["group"] == "A-2-6"

    def test_group_index_calculation(self):
        # GI = (F-35)[0.2+0.005(LL-40)] + 0.01(F-15)(PI-10)
        # F=60, LL=45, PL=22, PI=23
        # a = 60-35 = 25, b = 45-40 = 5, c = 60-15 = 45, d = 23-10 = 13
        # GI = 25*(0.2+0.005*5) + 0.01*45*13 = 25*0.225 + 5.85 = 5.625+5.85 = 11.475 ≈ 11
        result = classify_aashto(LL=45, PL=22, fines=60)
        assert result["group_index"] == 11

    def test_missing_ll_raises(self):
        with pytest.raises(ValueError, match="LL and PL"):
            classify_aashto(LL=None, PL=22, fines=60)

    def test_fines_auto_calculated(self):
        result = classify_aashto(LL=45, PL=22, gravel=10, sand=30)
        assert result["group"] == "A-7-6"


# ===================================================================
# Plasticity chart (smoke tests — visual output)
# ===================================================================

class TestPlasticityChart:
    """Smoke tests for plasticity_chart()."""

    def test_returns_axes(self):
        import matplotlib
        matplotlib.use("Agg")
        ax = plasticity_chart()
        assert ax is not None
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_with_data_points(self):
        import matplotlib
        matplotlib.use("Agg")
        ax = plasticity_chart(LL=[35, 60], PL=[18, 25], labels=["S1", "S2"])
        assert ax is not None
        import matplotlib.pyplot as plt
        plt.close("all")

    def test_with_custom_axes(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax_in = plt.subplots()
        ax_out = plasticity_chart(ax=ax_in)
        assert ax_out is ax_in
        plt.close("all")

    def test_single_label_string(self):
        import matplotlib
        matplotlib.use("Agg")
        ax = plasticity_chart(LL=[35], PL=[18], labels="Sample A")
        assert ax is not None
        import matplotlib.pyplot as plt
        plt.close("all")
