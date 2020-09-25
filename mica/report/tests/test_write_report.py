# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tempfile
import os
import shutil
import pytest
from warnings import warn

from testr.test_helper import on_head_network, has_sybase

from .. import report

try:
    import Ska.DBI
    user = os.environ.get('USER') or os.environ.get('LOGNAME')
    with Ska.DBI.DBI(server='sqlsao', dbi='sybase', user=user, database='axafvv') as db:
        HAS_SYBASE_ACCESS = True
except Exception:
    if on_head_network() and not has_sybase():
        warn("On HEAD but no sybase access. Run ska_envs or define SYBASE/SYBASE_OCS")
    HAS_SYBASE_ACCESS = False


HAS_SC_ARCHIVE = os.path.exists(report.starcheck.FILES['data_root'])


@pytest.mark.skipif('not HAS_SYBASE_ACCESS', reason='Report test requires Sybase VV access')
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
