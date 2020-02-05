# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tempfile
import os
import pytest

from .. import update_acq_stats as acq_stats

HAS_OBSPAR_ARCHIVE = os.path.exists(
        acq_stats.mica.archive.obspar.CONFIG['data_root'])


@pytest.mark.skipif('not HAS_OBSPAR_ARCHIVE', reason='Test requires mica obspars')
def test_calc_stats():
    acq_stats.calc_stats(17210)
    acq_stats.calc_stats(15175)
    acq_stats.calc_stats(4911)
    acq_stats.calc_stats(19386)


@pytest.mark.skipif('not HAS_OBSPAR_ARCHIVE', reason='Test requires mica obspars')
def test_make_acq_stats():
    """
    Save the acq stats for one obsid into a newly-created table
    """
    # Get a temporary file, but then delete it, because _save_acq_stats will only
    # make a new table if the supplied file doesn't exist
    fh, fn = tempfile.mkstemp(suffix='.h5')
    os.unlink(fn)
    acq_stats.table_file = fn
    obsid = 20001
    obsid_info, acq, star_info, catalog, temp = acq_stats.calc_stats(obsid)
    t = acq_stats.table_acq_stats(obsid_info, acq, star_info, catalog, temp)
    acq_stats._save_acq_stats(t)
    os.unlink(fn)
