"""Tests for geoeq.site.field_cbr — DCP and field CBR."""

import numpy as np
import pytest
from geoeq.site.field_cbr import dcp_cbr, field_cbr_test


class TestDcpCbr:
    def test_webster(self):
        # CBR = 10^(1.12 - 0.39*log10(10)) = 10^(1.12-0.39) = 10^0.73
        expected = 10**(1.12 - 0.39 * np.log10(10))
        assert dcp_cbr(10, method="webster") == pytest.approx(expected, rel=1e-3)

    def test_trl(self):
        expected = 10**(2.48 - 1.057 * np.log10(10))
        assert dcp_cbr(10, method="trl") == pytest.approx(expected, rel=1e-3)

    def test_kleyn(self):
        expected = 10**(2.62 - 1.27 * np.log10(10))
        assert dcp_cbr(10, method="kleyn") == pytest.approx(expected, rel=1e-3)

    def test_lower_dcpi_higher_cbr(self):
        cbr1 = dcp_cbr(5)
        cbr2 = dcp_cbr(50)
        assert cbr1 > cbr2

    def test_array_input(self):
        dcpi = np.array([5, 10, 20, 50])
        result = dcp_cbr(dcpi)
        assert len(result) == 4
        # Should be decreasing
        assert all(result[i] > result[i + 1] for i in range(3))

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            dcp_cbr(10, method="bogus")


class TestFieldCbrTest:
    def test_basic(self):
        pen = [0, 0.64, 1.27, 1.91, 2.54, 3.81, 5.08, 7.62, 10.16, 12.70]
        load = [0, 0.5, 1.2, 2.5, 4.0, 6.5, 9.0, 13.0, 16.0, 18.0]
        res = field_cbr_test(pen, load)
        assert res["CBR"] > 0
        assert res["CBR_2_54"] == pytest.approx(4.0 / 13.24 * 100, rel=1e-2)
        assert res["CBR_5_08"] == pytest.approx(9.0 / 19.96 * 100, rel=1e-2)

    def test_governing_cbr(self):
        pen = [0, 0.64, 1.27, 1.91, 2.54, 3.81, 5.08, 7.62, 10.16, 12.70]
        load = [0, 0.5, 1.2, 2.5, 4.0, 6.5, 9.0, 13.0, 16.0, 18.0]
        res = field_cbr_test(pen, load)
        assert res["CBR"] == max(res["CBR_2_54"], res["CBR_5_08"])

    def test_mismatched_lengths(self):
        with pytest.raises(ValueError):
            field_cbr_test([0, 1, 2], [0, 1])
