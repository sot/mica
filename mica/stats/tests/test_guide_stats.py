# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tempfile
import os
import numpy as np
import pytest
from pathlib import Path

from .. import guide_stats
from .. import update_guide_stats

HAS_GS_TABLE = os.path.exists(guide_stats.TABLE_FILE)


@pytest.mark.skipif("not HAS_GS_TABLE", reason="Test requires guide stats table")
def test_read_stats():
    stats = guide_stats.get_stats()
    slot = stats[(stats["obsid"] == 5438) & (stats["slot"] == 5)][0]
    np.isclose(slot["dy_std"], 0.16008819321078668)
    np.isclose(slot["dz_std"], 0.23807435722775061)
    # Fetch the stats for the slot 5 BOT by id
    single = guide_stats.get_star_stats(839386400)
    np.isclose(single[single["obsid"] == 5438][0]["dy_std"], 0.16008819321078668)


HAS_OBSPAR_ARCHIVE = (
    Path(update_guide_stats.mica.archive.obspar.CONFIG["data_root"]) / "00"
).exists()


@pytest.mark.skipif("not HAS_OBSPAR_ARCHIVE", reason="Test requires mica obspars")
def test_calc_stats():
    update_guide_stats.calc_stats(17210)


@pytest.mark.skipif("not HAS_OBSPAR_ARCHIVE", reason="Test requires mica obspars")
def test_calc_stats_with_bright_trans():
    s = update_guide_stats.calc_stats(17472)
    # Assert that the std on the slot 7 residuals are reasonable
    # even in this obsid that had a transition to BRIT
    assert s[1][7]["dr_std"] < 1


@pytest.mark.skipif("not HAS_OBSPAR_ARCHIVE", reason="Test requires mica obspars")
@pytest.mark.filterwarnings("ignore: object name")
def test_make_gui_stats():
    """
    Save the guide stats for one obsid into a newly-created table
    """
    # Get a temporary file, but then delete it, because _save_acq_stats will only
    # make a new table if the supplied file doesn't exist
    fh, fn = tempfile.mkstemp(suffix=".h5")
    os.close(fh)
    os.unlink(fn)
    obsid = 20001
    obsid_info, gui, star_info, catalog, temp = update_guide_stats.calc_stats(obsid)
    t = update_guide_stats.table_gui_stats(obsid_info, gui, star_info, catalog, temp)
    update_guide_stats._save_gui_stats(t, table_file=fn)
    os.unlink(fn)
