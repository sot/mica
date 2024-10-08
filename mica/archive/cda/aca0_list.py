# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Maintain a table with a list of of aca level 0 files in the Chandra Data Archive.

This module calls the CGI available from:

https://icxc.harvard.edu/dbtm/CDA/aspect_fetch.html

fetches the list of aca0 files, and updates a full list of files in 'MICA_ARCHIVE/aca0/cda_aca0.h5'.
That list should be of all available aca level 0 files in the Chandra Data Archive.  When called
from the update_aca_l0.py update script in the mica cron job, this module attempts to use the end
of the complete data in cda_aca0.h5 to form the query for new files, so only recent file names
are fetched.

The list of files in cda_aca0.h5 is compared against the list of files that have been archived in
MICA_ARCHIVE so that none are missed.
"""

import argparse
import logging
import os
import re
import time
import urllib

import numpy as np
import tables
from astropy.table import Table
from Chandra.Time import DateTime
from ska_helpers.retry import retry_call

from mica.common import MICA_ARCHIVE

CONFIG = dict(
    data_root=os.path.join(MICA_ARCHIVE, "aca0"),
    cda_fetch_url="https://icxc.harvard.edu/dbtm/CDA/cgi/aspect_fetch.cgi",
    cda_table="cda_aca0.h5",
)

logger = logging.getLogger("CDA ACA0 file list update")
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def get_options():
    parser = argparse.ArgumentParser(
        description="Update table of list of CDA ingested ACA0 files"
    )
    defaults = dict(CONFIG)
    parser.set_defaults(**defaults)
    parser.add_argument("--data-root", help="parent directory for all data")
    parser.add_argument("--cda-fetch-url", help="URL for CDA CGI")
    opt = parser.parse_args()
    return opt


def make_data_table(lines):
    files = Table.read(
        lines,
        format="ascii.csv",
        guess=False,
        names=["filename", "status", "ingest_time"],
    )
    # replace variable spaces with single spaces and then strptime
    ingest_dates = [
        time.strftime("%Y:%j:%H:%M:%S.000", time.strptime(t, "%b %d %Y %I:%M%p"))
        for t in [" ".join(f.split()) for f in files["ingest_time"]]
    ]
    files["ingest_date"] = ingest_dates
    file_re = re.compile(r"acaf(\d{9,})N(\d{3})_(\d)_img0.fits(\.gz)?")
    files["aca_ingest"] = np.repeat(DateTime().date, len(ingest_dates))
    files["filetime"] = [int(file_re.search(f).group(1)) for f in files["filename"]]
    files["version"] = [int(file_re.search(f).group(2)) for f in files["filename"]]
    files.sort(["aca_ingest", "ingest_date", "filename"])
    return files


TABLE_DTYPE = np.dtype(
    [
        ("filename", "S29"),
        ("status", "S1"),
        ("ingest_time", "S19"),
        ("ingest_date", "S21"),
        ("aca_ingest", "S21"),
        ("filetime", "<i8"),
        ("version", "<i8"),
    ]
)


def make_table_from_scratch(table_file, cda_fetch_url, start="2020:320:12:00:00"):
    logger.info("Fetching new CDA list from %s" % start)
    ct_start = DateTime(start)
    query = "?tstart={:02d}-{:02d}-{:04d}&pattern=acaimgc%%25&submit=Search".format(
        ct_start.mon, ct_start.day, ct_start.year
    )
    url = cda_fetch_url + query
    logger.info("URL for fetch {}".format(url))
    new_lines = urllib.request.urlopen(url).read().decode().splitlines()
    files = make_data_table(new_lines)
    logger.info("Creating new table at %s" % table_file)
    h5f = tables.open_file(
        table_file, "a", filters=tables.Filters(complevel=5, complib="zlib")
    )
    desc, bo = tables.table.descr_from_dtype(TABLE_DTYPE)
    tbl = h5f.create_table("/", "data", desc)
    for file in files:
        row = tbl.row
        row["filename"] = file["filename"]
        row["status"] = file["status"]
        row["ingest_time"] = file["ingest_time"]
        row["ingest_date"] = file["ingest_date"]
        row["aca_ingest"] = file["aca_ingest"]
        row["filetime"] = file["filetime"]
        row.append()
    tbl.flush()
    h5f.close()


def update_cda_table(data_root=None, cda_table=None, cda_fetch_url=None):
    if data_root is None:
        data_root = CONFIG["data_root"]
    if cda_table is None:
        cda_table = CONFIG["cda_table"]
    if cda_fetch_url is None:
        cda_fetch_url = CONFIG["cda_fetch_url"]

    if not os.path.exists(data_root):
        os.makedirs(data_root)
    table_file = os.path.join(data_root, cda_table)

    if not os.path.exists(table_file):
        make_table_from_scratch(table_file, cda_fetch_url)
        return

    # Load the existing table
    with tables.open_file(table_file, "r") as h5f:
        tbl = h5f.get_node("/", "data")
        cda_files = Table(tbl[:])
        lastdate = DateTime(cda_files[-1]["ingest_date"])

    # Get new data
    logger.info("Fetching new CDA list from %s" % lastdate.date)
    query = "?tstart={:02d}-{:02d}-{:04d}&pattern=acaimgc%%25&submit=Search".format(
        lastdate.mon, lastdate.day, lastdate.year
    )
    url = cda_fetch_url + query
    logger.info("URL for fetch {}".format(url))
    try:
        resp = retry_call(urllib.request.urlopen, [url], tries=5, delay=5)
        new_lines = resp.read().decode().splitlines()
    except urllib.error.URLError as err:
        logger.info(err)
        return
    files = make_data_table(new_lines)
    match = np.flatnonzero(cda_files["filename"] == files[0]["filename"])
    if len(match) == 0:
        raise ValueError("no overlap")
    match_last_idx = match[-1]
    i_diff = 0
    for have_entry, new_entry in zip(cda_files[match_last_idx:], files):
        if have_entry["filename"] != new_entry["filename"]:
            break
        i_diff += 1

    if i_diff < len(files):
        with tables.open_file(table_file, "a") as h5f:
            tbl = h5f.get_node("/", "data")
            logger.info(
                "Updating %s with %d new rows" % (table_file, len(files[i_diff:]))
            )
            for file in files[i_diff:]:
                logger.debug(file)
                row = tbl.row
                row["filename"] = file["filename"]
                row["status"] = file["status"]
                row["ingest_time"] = file["ingest_time"]
                row["ingest_date"] = file["ingest_date"]
                row["aca_ingest"] = file["aca_ingest"]
                row["filetime"] = file["filetime"]
                row.append()
            tbl.flush()


def main():
    opt = get_options()
    kwargs = vars(opt)
    update_cda_table(**kwargs)


if __name__ == "__main__":
    main()
