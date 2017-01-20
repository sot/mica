import os
import re
import tables
import asciitable
import numpy as np
from six.moves import urllib, zip

import time
from Chandra.Time import DateTime
import Ska.Numpy
import logging
import argparse

from mica.common import MICA_ARCHIVE

CONFIG = dict(data_root=os.path.join(MICA_ARCHIVE, 'aca0'),
              cda_fetch_url='https://icxc.harvard.edu/dbtm/CDA/cgi/aspect_fetch.cgi',
              cda_table='cda_aca0.h5')

logger = logging.getLogger('CDA ACA0 file list update')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def get_options():
    parser = argparse.ArgumentParser(
        description="Update table of list of CDA ingested ACA0 files")
    defaults = dict(CONFIG)
    parser.set_defaults(**defaults)
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--cda-fetch-url",
                        help="URL for CDA CGI")
    opt = parser.parse_args()
    return opt


def make_data_table(lines):
    files = asciitable.read(lines,
                            names=['filename', 'status', 'ingest_time'])
    # replace variable spaces with single spaces and then strptime
    ingest_dates = [time.strftime(
            "%Y:%j:%H:%M:%S.000",
            time.strptime(
                t, '%b %d %Y %I:%M%p'))
                    for t in
                    [' '.join(f.split())
                     for f in files['ingest_time']]]
    file_re = re.compile(r'acaf(\d{9,})N(\d{3})_(\d)_img0.fits(\.gz)?')
    filetimes = [int(file_re.search(f).group(1)) for f in files['filename']]
    versions = [int(file_re.search(f).group(2)) for f in files['filename']]
    now_dates = np.repeat(DateTime().date, len(ingest_dates))
    files = Ska.Numpy.add_column(files, 'ingest_date', ingest_dates)
    files = Ska.Numpy.add_column(files, 'aca_ingest', now_dates)
    files = Ska.Numpy.add_column(files, 'filetime', filetimes)
    files = Ska.Numpy.add_column(files, 'version', versions)
    files.sort(order=['aca_ingest', 'ingest_date', 'filename'])
    return files


def make_table_from_scratch(table_file, cda_fetch_url, start='2015:001'):
    logger.info("Fetching new CDA list from %s" % start)
    ct_start = DateTime(start)
    query = ("?tstart={:02d}-{:02d}-{:04d}&pattern=acaimgc%%25&submit=Search".format(
            ct_start.mon, ct_start.day, ct_start.year))
    url = cda_fetch_url + query
    new_lines = urllib.request.urlopen(url).readlines()
    files = make_data_table(new_lines)
    logger.info("Creating new table at %s" % table_file)
    h5f = tables.openFile(table_file, 'a',
                          filters=tables.Filters(complevel=5, complib='zlib'))
    desc, bo = tables.table.descr_from_dtype(files[0].dtype)
    tbl = h5f.createTable('/', 'data', desc)
    tbl.append(files)
    tbl.flush()
    h5f.close()


def update_cda_table(data_root=None,
                     cda_table=None,
                     cda_fetch_url=None):
    if data_root is None:
        data_root = CONFIG['data_root']
    if cda_table is None:
        cda_table = CONFIG['cda_table']
    if cda_fetch_url is None:
        cda_fetch_url = CONFIG['cda_fetch_url']

    if not os.path.exists(data_root):
        os.makedirs(data_root)
    table_file = os.path.join(data_root, cda_table)
    if not os.path.exists(table_file):
        make_table_from_scratch(table_file, cda_fetch_url)

    h5f = tables.openFile(table_file, 'a')
    tbl = h5f.getNode('/', 'data')
    cda_files = tbl[:]
    lastdate = DateTime(cda_files[-1]['ingest_date'])

    logger.info("Fetching new CDA list from %s" % lastdate.date)
    query = ("?tstart={:02d}-{:02d}-{:04d}&pattern=acaimgc%%25&submit=Search".format(
            lastdate.mon, lastdate.day, lastdate.year))
    url = cda_fetch_url + query
    new_lines = urllib.request.urlopen(url).readlines()
    files = make_data_table(new_lines)
    match = np.flatnonzero(cda_files['filename'] == files[0]['filename'])
    if len(match) == 0:
        raise ValueError("no overlap")
    match_last_idx = match[-1]

    i_diff = 0
    for have_entry, new_entry in zip(cda_files[match_last_idx:], files):
        if have_entry['filename'] != new_entry['filename']:
            break
        i_diff += 1

    if i_diff < len(files):
        logger.info("Updating %s with %d new rows"
                    % (table_file, len(files[i_diff:])))
        for file in files[i_diff:]:
            logger.debug(file)
            row = tbl.row
            row['filename'] = file['filename']
            row['status'] = file['status']
            row['ingest_time'] = file['ingest_time']
            row['ingest_date'] = file['ingest_date']
            row['aca_ingest'] = file['aca_ingest']
            row['filetime'] = file['filetime']
            row.append()
        tbl.flush()
    h5f.close()


def main():
    opt = get_options()
    kwargs = vars(opt)
    update_cda_table(**kwargs)


if __name__ == '__main__':
    main()
