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
    update_table(opt.datafile)


def update_table(datafile, filter=None):
    table = get_web_table(filter)
    datafile = Path(datafile)
    if datafile.suffix == '.h5':
        table.write(datafile, path='root', serialize_meta=True, overwrite=True)
    else:
        table.write(datafile, overwrite=True)


def get_web_table(filter=None):
    """
    Fetch all observation details from the SRServices ocatDetails page

    "https://cda.harvard.edu/srservices/ocatDetails.do?format=text"

    :filter: string to include in web query to limit results
    :return: astropy table of the observation details
    """

    url = URL_OCATDETAILS
    if filter is not None:
        url = f"{url}&{filter}"

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

    # Munge length descriptions back to standard RDB (from "20" to "20S" etc)
    # while leaving the numeric types alone ("20N" stayes "20N").
    for j, ctype in enumerate(ctypes):
        if not ctype.endswith("N"):
            ctypes[j] = ctype + "S"
    lines[header_start + 1] = "\t".join(ctypes)

    dat = Table.read(lines, format='ascii.rdb')

    # Set encoding to utf-8
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'U':
            dat[name] = np.char.encode(col, 'utf-8')

    # Lower-case all the columns
    lc_names = [name.lower() for name in dat.colnames]
    dat.rename_columns(dat.colnames, lc_names)

    return dat
