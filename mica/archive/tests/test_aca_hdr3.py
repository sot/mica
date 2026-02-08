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
    }


def test_two_byte_sum():
    bytes0 = np.ma.array([0x00, 0xF0, 0x0F, 0xFF, 0xFF], dtype=np.uint8)
    bytes1 = np.ma.array([0x00, 0x0F, 0xF0, 0xFF, 0xFF], dtype=np.uint8)
    bytes0[-1] = np.ma.masked
    bytes1[-1] = np.ma.masked

    # Original code prior to PR #315
    out1 = (
        (bytes0.astype("int") >> 7) * (-1 * 65535)
        + (bytes0.astype("int") << 8)
        + (bytes1.astype("int"))
    )
    assert np.all(out1 == np.ma.array([0, -4080, 4080, 0, 0], mask=[0, 0, 0, 0, 1]))

    # New code in PR #315
    slot_data = Table([bytes0, bytes1], names=["byte0", "byte1"])
    ints16 = aca_hdr3.two_byte_sum(["byte0", "byte1"])(slot_data)

    assert np.all(ints16 == np.ma.array([0, -4081, 4080, -1, 0], mask=[0, 0, 0, 0, 1]))
