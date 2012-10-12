import os
import tables
import asciitable
import numpy as np
import urllib
import mx.DateTime
from itertools import izip
from Chandra.Time import DateTime
import logging
import argparse

config = dict(data_root='/data/aca/archive/aca0',
              cda_fetch_url='https://icxc.harvard.edu/dbtm/CDA/cgi/aspect_fetch.cgi',
              cda_table='cda_aca0.h5')

logger = logging.getLogger('CDA ACA0 file list update')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

def get_options():
    parser = argparse.ArgumentParser(
        description="Update table of list of CDA ingested ACA0 files")
    defaults = dict(config)
    parser.set_defaults(**defaults)
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--cda-fetch-url",
                        help="URL for CDA CGI")
    opt = parser.parse_args()
    return opt

def make_table_from_scratch(table_file, cda_fetch_url, start='1999:001'):
    logger.info("Fetching new data from %s" % start)
    startmx = DateTime(start).mxDateTime
    query = ("?tstart=%s&pattern=aca%%25&submit=Search"
             % startmx.strftime("%m-%d-%Y"))
    url = cda_fetch_url + query
    new_lines = urllib.urlopen(url).readlines()
    files = asciitable.read(new_lines,
                            names=['filename', 'status', 'ingest_time'])
    logger.info("Creating new table at %s" % table_file)
    h5f = tables.openFile(table_file, 'a',
                          filters=tables.Filters(complevel=5, complib='zlib'))
    desc, bo = tables.table.descr_from_dtype(files[0].dtype)
    tbl = h5f.createTable('/', 'data', desc)
    tbl.append(files)
    tbl.flush()
    h5f.close()


def update_cda_table(data_root=config['data_root'],
         cda_table=config['cda_table'],
         cda_fetch_url=config['cda_fetch_url']):

    table_file = os.path.join(data_root, cda_table)
    if not os.path.exists(table_file):
        make_table_from_scratch(table_file, cda_fetch_url)
    
    h5f = tables.openFile(table_file, 'a')
    tbl = h5f.getNode('/', 'data')
    cda_files = tbl[:]
    lasttime = cda_files[-1]['ingest_time']

    mxlast = mx.DateTime.strptime(lasttime, "%b%t%d%t%Y%t%I:%M%p")
    logger.info("Fetching new data from %s" 
                % mxlast.strftime("%m-%d-%Y"))
    query = ("?tstart=%s&pattern=aca%%25&submit=Search"
                 % mxlast.strftime("%m-%d-%Y"))
    url = cda_fetch_url + query
    new_lines = urllib.urlopen(url).readlines()
    arc_files = asciitable.read(new_lines,
                                names = ['filename', 'status', 'ingest_time'])
    match_first = np.flatnonzero(cda_files['filename'] == arc_files[0]['filename'])
    if len(match_first) == 0:
        raise ValueError("no overlap")
    if len(match_first) > 1:
        raise ValueError("duplicate files")
    match_first_idx = match_first[0]

    i_diff = 0
    for have_entry, new_entry in izip(cda_files[match_first_idx:], arc_files):
        if have_entry['filename'] != new_entry['filename']:
            break
        i_diff += 1

    if i_diff < len(arc_files):
        logger.info("Updating %s with %d new rows" 
                    % (table_file, len(arc_files[i_diff:])))
        for file in arc_files[i_diff:]:
            logger.debug(file)
            row = tbl.row
            row['filename'] = file['filename']
            row['status'] = file['status']
            row['ingest_time'] = file['ingest_time']
            row.append()
        tbl.flush()
    h5f.close()

def main():
    opt = get_options()
    kwargs = vars(opt)
    update_cda_table(**kwargs)


if __name__ == '__main__':
    main()
