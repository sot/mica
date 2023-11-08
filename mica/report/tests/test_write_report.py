# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tempfile
import os
from pathlib import Path
import getpass
import shutil
import pytest
from warnings import warn

from testr.test_helper import on_head_network, has_sybase

from .. import report

user = getpass.getuser()

try:
    import Ska.DBI
    with Ska.DBI.DBI(server='sqlsao', dbi='sybase', user=user, database='axafvv') as db:
        HAS_SYBASE_ACCESS = True
except Exception:
    HAS_SYBASE_ACCESS = False

    # If the user should have access, warn about the issue.
    if (on_head_network() and not has_sybase() and
            Path(os.environ['SKA'], 'data', 'aspect_authorization',
                 f'sqlsao-axafvv-{user}').exists()):
        warn("On HEAD but no sybase access. Run test from production environment.")


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
    for obsid in [20001, 15175, 54778, 44077]:
        report.main(obsid)
    os.unlink(fn)
    shutil.rmtree(tempdir)
