# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Update a mica Ska target table
"""
import argparse
from pathlib import Path
from astropy.table import Table
import requests
import numpy as np

from ska_helpers.retry import retry_call
from mica.common import MICA_ARCHIVE

URL_OCATDETAILS = "https://cda.harvard.edu/srservices/ocatDetails.do?format=text"

def get_options():
    parser = argparse.ArgumentParser(
        description="Update target table")
    parser.add_argument("--datafile",
                        default='target_table.h5')
    opt = parser.parse_args()
    return opt


def main():
    opt = get_options()
    datafile = Path(opt.datafile)
    update_table(datafile)


def update_table(datafile):
    table = get_web_table()
    if datafile.name.endswith('h5'):
        table.write(datafile, path='root', serialize_meta=True, overwrite=True)
    else:
        table.write(datafile, overwrite=True)


def get_web_table():
    """
    Fetch all observation details from the SRServices ocatDetails page

    "https://cda.harvard.edu/srservices/ocatDetails.do?format=text"

    :return: astropy table of the observation details
    """

    url = URL_OCATDETAILS

    try:
        resp = retry_call(requests.get, [url], {"timeout": 120},
                          tries=4, delay=1)
    except Exception:
        return None
    if not resp.ok:
        return None

    lines = resp.text.split("\n")

    # Convert the type line to standard RDB
    # First find the line that begins the column descriptions
    for i, line in enumerate(lines):
        if line.startswith('SEQ_NUM'):
            header_start = i
            break

    # The line with the lengths and types is next (header_start + 1)
    ctypes = lines[header_start + 1].split("\t")
    for j, t in enumerate(ctypes):
        if not t.endswith("N"):
            ctypes[j] = t + "S"
    lines[header_start + 1] = "\t".join(ctypes)

    dat = Table.read(lines, format='ascii.rdb')

    # Set encoding to utf-8
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'U':
            dat[name] = np.char.encode(col, 'utf-8')

    return dat


if __name__ == "__main__":
    main()
