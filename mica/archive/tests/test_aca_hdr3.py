"""
Basic functionality and regression tests for ACA hdr3 (diagnostic) telemetry.
"""

import numpy as np

from .. import aca_hdr3


def test_MSIDset():
    """
    Read all available MSIDs into a single MSIDset.  Use the empirically determined
    lengths as regression tests.
    """
    msids = [hdr3['msid'] for hdr3 in aca_hdr3.HDR3_DEF.values() if 'value' in hdr3]

    # Read all MSIDs as a set
    dat = aca_hdr3.MSIDset(msids, '2010:001', '2010:003')

    val_lengths = np.array([len(dat[msid].vals) for msid in msids])
    time_lengths = np.array([len(dat[msid].times) for msid in msids])
    assert np.all(val_lengths == time_lengths)
    assert np.all(val_lengths == 44432)

    for msid in msids:
        dat[msid].filter_bad()
    val_lengths = np.array([len(dat[msid].vals) for msid in msids])
    time_lengths = np.array([len(dat[msid].times) for msid in msids])
    assert np.all(val_lengths == time_lengths)
    assert np.all(val_lengths == [40528, 44432, 44432, 10679, 44432, 10760, 10731, 44432,
                                  44432, 10679, 44432, 44432, 44432, 44432, 40991, 44432,
                                  40991, 44432, 44432, 40991, 10679])
