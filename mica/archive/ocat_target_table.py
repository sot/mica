# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from astropy.table import Table
from pathlib import Path
from mica.common import MICA_ARCHIVE

OCAT_TARGET_TABLE = Path(MICA_ARCHIVE) / 'ocat_target_table.h5'

def get_ocat_target_table(datafile=None):
    """
    Read the ocat target table h5 data product.

    :param datafile: h5 ocat target table data file.
                     Defaults to MICA_ARCHIVE/ocat_target_table.h5
    :returns: astropy table of target table
    """
    if datafile is None:
        datafile = OCAT_TARGET_TABLE

    dat = Table.read(datafile)

    # Decode bytes to strings manually.  Perhaps fixed in newer versions of numpy.
    # Eventually we will want just dat.convert_bytestring_to_unicode().
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'S':
            dat[name] = np.char.decode(col, 'utf-8')

    return dat
