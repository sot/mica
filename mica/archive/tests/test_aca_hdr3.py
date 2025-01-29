# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Basic functionality and regression tests for ACA hdr3 (diagnostic) telemetry.
"""

import os

import numpy as np
import pytest

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

    # Read all MSIDs as a set
    dat = aca_hdr3.MSIDset(msids, "2010:001:12:00:00", "2010:003:12:00:00")

    val_lengths = np.array([len(dat[msid].vals) for msid in msids])
    time_lengths = np.array([len(dat[msid].times) for msid in msids])
    assert np.all(val_lengths == time_lengths)
    assert np.all(val_lengths == 44432)

    for msid in msids:
        dat[msid].filter_bad()
    val_lengths = np.array([len(dat[msid].vals) for msid in msids])
    time_lengths = np.array([len(dat[msid].times) for msid in msids])
    assert np.all(val_lengths == time_lengths)
    assert np.all(
        val_lengths
        == [
            10679,
            40991,
            40991,
            40528,
            40514,
            40514,
            40991,
            40991,
            40514,
            40991,
            40514,
            40514,
            40991,
            40514,
            10731,
            40528,
            40528,
            40528,
            10679,
            10760,
            10679,
        ]
    )


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
    bytes8_2xN = np.ma.vstack([bytes0, bytes1], dtype=np.uint8)
    bytes8 = bytes8_2xN.transpose().flatten().copy()
    ints16 = np.ma.array(bytes8.data.view(">i2"), mask=bytes8.mask[::2])

    assert np.all(ints16 == np.ma.array([0, -4081, 4080, -1, 0], mask=[0, 0, 0, 0, 1]))
