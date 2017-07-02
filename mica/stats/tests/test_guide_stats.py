import tempfile
import os

from .. import guide_stats


def test_calc_stats():
    guide_stats.calc_stats(17210)


def test_make_gui_stats():
    """
    Save the guide stats for one obsid into a newly-created table
    """
    # Get a temporary file, but then delete it, because _save_acq_stats will only
    # make a new table if the supplied file doesn't exist
    fh, fn = tempfile.mkstemp(suffix='.h5')
    os.unlink(fn)
    guide_stats.TABLE_FILE = fn
    obsid = 20001
    obsid_info, gui, star_info, catalog, temp = guide_stats.calc_stats(obsid)
    t = guide_stats.table_gui_stats(obsid_info, gui, star_info, catalog, temp)
    guide_stats._save_gui_stats(t)
    os.unlink(fn)
