# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tempfile
import os
from pathlib import Path
import shutil
import pytest
from warnings import warn

from .. import update_ocat_target_table, ocat_target_table

def test_write_read_table():
    """
    Make a report and database
    """
    tempdir = tempfile.mkdtemp()
    file = Path(tempdir, 'target_table.h5')
    update_ocat_target_table.update_table(datafile=file)
    dat = ocat_target_table.get_ocat_target_table(datafile=file)
    ok = dat['obsid'] == 2121
    dat['target_name'][ok][0] == "MCG-5-23-16"
    if file.exists():
        os.unlink(file)
        shutil.rmtree(tempdir)
