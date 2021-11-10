# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import tables
from astropy.table import Table
from pathlib import Path
from mica.common import MICA_ARCHIVE

OCAT_TARGET_TABLE = Path(MICA_ARCHIVE) / 'ocat_target_table.h5'


def get_ocat_target_table(datafile=None, read_where=None):
    """
    Read the ocat target table h5 data product.

    :param datafile: h5 ocat target table data file.
                     Defaults to MICA_ARCHIVE/ocat_target_table.h5
    :param read_where: filter string to pass to tables read_where() to
                   limit returned results.  See
                   https://www.pytables.org/usersguide/condition_syntax.html

    :returns: astropy table of target table
    """
    if datafile is None:
        datafile = OCAT_TARGET_TABLE

    if read_where is None:
        dat = Table.read(datafile)
    else:
        with tables.open_file(datafile) as hdu:
            recs = hdu.root.root.read_where(read_where)
            dat = Table(recs)

    # Decode bytes to strings manually.  Fixed in numpy 1.20.
    # Eventually we will want just dat.convert_bytestring_to_unicode().
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'S':
            dat[name] = np.char.decode(col, 'utf-8')

    return dat
