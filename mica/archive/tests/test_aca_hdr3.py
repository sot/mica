# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Basic functionality and regression tests for ACA hdr3 (diagnostic) telemetry.
"""

import os

import numpy as np
import pytest
from astropy.table import Table

from mica.archive import aca_hdr3

has_l0_2010_archive = os.path.exists(
    os.path.join(aca_hdr3.aca_l0.CONFIG["data_root"], "2010")
)

has_l0_2026_archive = os.path.exists(
    os.path.join(aca_hdr3.aca_l0.CONFIG["data_root"], "2026")
)


@pytest.mark.skipif("not has_l0_2010_archive", reason="Test requires 2010 L0 archive")
def test_MSIDset():
    """
    Read all available MSIDs into a single MSIDset.  Use the empirically determined
    lengths as regression tests.
    """
    msids = [hdr3["msid"] for hdr3 in aca_hdr3.HDR3_DEF.values() if "value" in hdr3]
    msids = sorted(msids)

    # Read all MSIDs as a set. This is from old times when there was a mix of 6x6 and
    # 8x8 data so some MSIDs are only sampled for a shorter time.
    dat = aca_hdr3.MSIDset(msids, "2010:001:12:00:00", "2010:003:12:00:00")
    lengths = {msid: len(dat[msid].vals) for msid in msids}
    assert lengths == {
        "aca_temp": 10679,
        "ad_15v_ps": 40991,
        "ad_27v_ps": 40991,
        "ad_5v_ps": 40528,
        "ad_achhs_therm": 40514,
        "ad_achohs_therm": 40514,
        "ad_analog_gnd": 40991,
        "ad_converter_therm": 40991,
        "ad_lc_therm": 40514,
        "ad_m15v_ps": 40991,
        "ad_pmhs_therm": 40514,
        "ad_pmohs_therm": 40514,
        "ad_smhs_therm": 40991,
        "ad_smohs_therm": 40514,
        "avg_bkg": 10731,
        "ccd_det_therm": 40528,
        "ccd_molyb_therm_1": 40528,
        "ccd_molyb_therm_2": 40528,
        "ccd_setpoint": 10679,
        "ccd_temp": 10760,
        "dac": 10679,
        "zero_off16_quad_a": 10710,
        "zero_off16_quad_b": 10710,
        "zero_off16_quad_c": 10710,
        "zero_off16_quad_d": 10267,
        "zero_off32_quad_a": 10267,
        "zero_off32_quad_b": 10267,
        "zero_off32_quad_c_lsb": 10731,
        "zero_off32_quad_c_msb": 10267,
        "zero_off32_quad_d": 10731,
    }


def test_two_byte_sum():
    bytes0 = np.array([0x00, 0xF0, 0x0F, 0xFF, 0x00], dtype=np.uint8)
    bytes1 = np.array([0x00, 0x0F, 0xF0, 0xFF, 0xFF], dtype=np.uint8)
    exp = np.array([0, -4081, 4080, -1, 255])

    slot_data = Table([bytes0, bytes1], names=["byte0", "byte1"])
    vals = aca_hdr3.two_byte_sum(["byte0", "byte1"])(slot_data)
    assert vals.dtype == np.int16
    assert np.all(vals == exp)

    scale = 2.5
    vals = aca_hdr3.two_byte_sum(["byte0", "byte1"], scale=scale)(slot_data)
    assert vals.dtype == np.float64
    assert np.all(vals == scale * exp)


@pytest.mark.skipif("not has_l0_2026_archive", reason="Test requires 2026 L0 archive")
def test_quad_offset_16_vs_32():
    """Test that 16-bit and 32-bit quadrant offsets agree to within 1.0 A/D count.

    The 32-bit value is a higher-resolution version of the same measurement, so the
    16-bit rounded value (interpolated to the 32-bit timestamps) should differ by less
    than 1 LSB (1 A/D converter count).
    """
    start, stop = "2026:021:02:00:00", "2026:021:03:00:00"
    for quad in ("a", "b", "c", "d"):
        dat16 = aca_hdr3.MSID(f"zero_off16_quad_{quad}", start, stop)
        dat32 = aca_hdr3.MSID(f"zero_off32_quad_{quad}", start, stop)

        # Interpolate the 16-bit values onto the 32-bit timestamps
        vals16_interp = np.interp(dat32.times, dat16.times, dat16.vals)

        diff = np.abs(vals16_interp - dat32.vals)
        assert np.all(diff < 1.0), (
            f"quad {quad}: max diff {diff.max():.4f} exceeds 1.0 A/D count"
        )


# ---------------------------------------------------------------------------
# Unit tests for _fuzzy_join_times
# ---------------------------------------------------------------------------


def test_fuzzy_join_times_exact_match():
    """Identical arrays: every element matches its counterpart exactly."""
    msb = np.array([10.0, 20.0, 30.0])
    lsb = np.array([10.0, 20.0, 30.0])
    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb)
    np.testing.assert_array_equal(idx_msb, [0, 1, 2])
    np.testing.assert_array_equal(idx_lsb, [0, 1, 2])


def test_fuzzy_join_times_lsb_slightly_ahead():
    """LSB times slightly ahead of MSB (searchsorted lands on the correct bin directly)."""
    msb = np.array([10.0, 20.0, 30.0])
    lsb = np.array([10.5, 20.5, 30.5])
    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb)
    np.testing.assert_array_equal(idx_msb, [0, 1, 2])
    np.testing.assert_array_equal(idx_lsb, [0, 1, 2])


def test_fuzzy_join_times_lsb_slightly_behind():
    """LSB times slightly behind MSB (exercises the prev-candidate-is-closer branch)."""
    msb = np.array([10.0, 20.0, 30.0])
    lsb = np.array([9.5, 19.5, 29.5])
    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb)
    np.testing.assert_array_equal(idx_msb, [0, 1, 2])
    np.testing.assert_array_equal(idx_lsb, [0, 1, 2])


def test_fuzzy_join_times_gap_excludes_entry():
    """An MSB entry with no LSB neighbor within tolerance is dropped."""
    msb = np.array([10.0, 20.0, 30.0])
    lsb = np.array([10.0, 25.0, 30.0])  # lsb[1] is 5 s from msb[1] > tol
    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb)
    np.testing.assert_array_equal(idx_msb, [0, 2])
    np.testing.assert_array_equal(idx_lsb, [0, 2])


def test_fuzzy_join_times_tolerance_boundary():
    """Entry at exactly tol is included; entry just over tol is excluded."""
    msb = np.array([10.0, 20.0, 30.0])
    lsb = np.array(
        [10.0, 22.0, 33.0]
    )  # msb[1]-lsb[1]=2.0 (in), msb[2]-lsb[2]=3.0 (out)
    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb, tol=2.0)
    np.testing.assert_array_equal(idx_msb, [0, 1])
    np.testing.assert_array_equal(idx_lsb, [0, 1])


def test_fuzzy_join_times_uniform_offset():
    """MSB and LSB staggered by 1 s; all pairs match with tol=1.1, none with tol=0.9."""
    msb = np.array([0.0, 2.0, 4.0, 6.0, 8.0])
    lsb = np.array([1.0, 3.0, 5.0, 7.0, 9.0])
    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb, tol=1.1)
    np.testing.assert_array_equal(idx_msb, [0, 1, 2, 3, 4])
    np.testing.assert_array_equal(idx_lsb, [0, 1, 2, 3, 4])

    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb, tol=0.9)
    assert len(idx_msb) == 0
    assert len(idx_lsb) == 0


def test_fuzzy_join_times_unequal_lengths():
    """MSB and LSB have different lengths with a gap; multiple MSB can match same LSB.

    msb = [0, 2, 4, 6, 8, 12, 14], lsb = [1, 3, 7.06, 9.04, 10.8, 13.15], tol = 1.1
      msb[0]=0  -> lsb[0]=1     (diff 1.0)
      msb[1]=2  -> lsb[1]=3     (diff 1.0)
      msb[2]=4  -> lsb[1]=3     (diff 1.0; lsb[2]=7.06 is farther)
      msb[3]=6  -> lsb[2]=7.06  (diff 1.06)
      msb[4]=8  -> lsb[2]=7.06  (diff 0.94; lsb[3]=9.04 is farther)
      msb[5]=12 -> lsb[5]=13.15 (diff 1.15 > tol; no match)
      msb[6]=14 -> lsb[5]=13.15 (diff 0.85)
    """
    msb = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 12.0, 14.0])
    lsb = np.array([1.0, 3.0, 7.06, 9.04, 10.8, 13.15])
    idx_lsb, idx_msb = aca_hdr3._fuzzy_join_times(msb, lsb, tol=1.1)
    np.testing.assert_array_equal(idx_msb, [0, 1, 2, 3, 4, 6])
    np.testing.assert_array_equal(idx_lsb, [0, 1, 1, 2, 2, 5])
