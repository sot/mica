# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tempfile
import os
import shutil
import pytest

from .. import report

try:
    import Ska.DBI
    with Ska.DBI.DBI(server='sqlsao', dbi='sybase', user='aca_ops', database='axafocat') as db:
        assert db.conn._is_connected == 1
        HAS_SYBASE_ACCESS = True
except:
    HAS_SYBASE_ACCESS = False


HAS_SC_ARCHIVE = os.path.exists(report.starcheck.FILES['data_root'])


@pytest.mark.skipif('not HAS_SYBASE_ACCESS', reason='Report test requires Sybase/OCAT access')
@pytest.mark.skipif('not HAS_SC_ARCHIVE', reason='Report test requires mica starcheck archive')
def test_write_reports():
    """
    Make a report and database
    """
    tempdir = tempfile.mkdtemp()
    # Get a temporary file, but then delete it, because report.py will only
    # make a new table if the supplied file doesn't exist
    fh, fn = tempfile.mkstemp(dir=tempdir, suffix='.db3')
    os.unlink(fn)
    report.REPORT_ROOT = tempdir
    report.REPORT_SERVER = fn
    for obsid in [20001, 15175, 54778]:
        report.main(obsid)
    os.unlink(fn)
    shutil.rmtree(tempdir)
