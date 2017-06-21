import tempfile
import os

from .. import acq_stats


def test_calc_stats():
    acq_stats.calc_stats(17210)


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
