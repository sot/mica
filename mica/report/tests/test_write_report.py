# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tempfile
import os
import shutil

from .. import report

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
