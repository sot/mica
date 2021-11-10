# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from .. import update_ocat_target_table, ocat_target_table


def test_write_read_table(tmp_path):
    """
    Make a report and database
    """
    file = tmp_path / 'target_table.h5'
    update_ocat_target_table.update_table(datafile=file,
                                          filter="obsid=2115-2125")
    dat = ocat_target_table.get_ocat_target_table(datafile=file)
    ok = dat['obsid'] == 2121
    assert np.all(dat['target_name'][ok] == "MCG-5-23-16")
    dat = ocat_target_table.get_ocat_target_table(datafile=file,
                                                  read_where='obsid==2115')
    assert dat['target_name'] == 'LBQS 2350-0045A'
