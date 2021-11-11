# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pathlib import Path
import re
import warnings

import requests
import numpy as np
import tables
from astropy.table import Table

from mica.common import MICA_ARCHIVE

__all__ = ['get_ocat_target_table', 'get_cda_archive_file_list', 'get_cda_prop_abstract',
           'get_cda_ocat']

OCAT_TARGET_TABLE = Path(MICA_ARCHIVE) / 'ocat_target_table.h5'
URL_CDA_SERVICES = "https://cda.harvard.edu/srservices"
CDA_SERVICES = {
    'prop_abstract': 'propAbstract',
    'ocat_summary': 'ocatList',
    'ocat_details': 'ocatDetails',
    'archive_file_list': 'archiveFileList'}


def html_to_text(html):
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(html, features='lxml')
    text = soup.get_text()
    text = re.sub(r'\n+', '\n', text)
    return text


def clean_text(text):
    out = text.encode('ascii', errors='ignore').decode('ascii')
    out = out.replace('\n', ' ').replace('\r', '').strip()
    return out


def get_cda_archive_file_list(obsid, detector, level, dataset='flight', **params):
    """
    Get list of archive files for given ``obsid``, ``detector``, ``level``, and ``dataset``.

    Other parameters can be ``subdetector``, ``filetype``, ``filename``, and ``obi``.

    Note: this may not be working for level 0 products.

    Examples::

       >>> get_cda_archive_file_list(obsid=2365, detector='pcad',
       ...                           subdetector='aca', level=1, obi=2)
       <Table length=27>
               Filename            Filesize      Timestamp
                   str30               int64          str19
       ------------------------------ -------- -------------------
       pcadf126690624N007_asol1.fits  7300800 2021-04-09 08:04:29
       pcadf02365_002N001_asol1.fits  4728960 2021-04-09 08:04:30
                               ...      ...                 ...
       pcadf126695890N007_adat61.fits  1293120 2021-04-09 08:04:28
       pcadf126695890N007_adat71.fits  1293120 2021-04-09 08:04:28

       >>> get_cda_archive_file_list(obsid=400, detector='acis', level=2, filetype='evt2')
       <Table length=1>
               Filename         Filesize      Timestamp
                str24            int64          str19
       ------------------------ -------- -------------------
       acisf00400N007_evt2.fits  4619520 2011-07-08 13:52:57

    :param obsid: int, str
        Observation ID
    :param detector: str
        Detector name (e.g. 'pcad', 'acis')
    :param level: int, float, str
        Level name (e.g. 0, 0.5, 1, 1.5, 2, 3)
    :param dataset: str
        Dataset name (default='flight')
    :param **params: dict
        Additional parameters to filter query (subdetector, filetype, obi, filename)
    :return: astropy Table
        Table of archive files
    """
    params['dataset'] = dataset
    params['detector'] = detector
    params['level'] = level
    params['obsid'] = obsid

    text = get_cda_service_text('archive_file_list', **params)
    dat = Table.read(text.splitlines(), format='ascii.basic', delimiter='\t', guess=False)
    filesize = [int(x.replace(',', '')) for x in dat['Filesize']]
    dat['Filesize'] = filesize

    return dat


def get_cda_prop_abstract(obsid=None, propnum=None, timeout=30):
    """Get a proposal abstract from the CDA services.

    One of ``obsid`` or ``propnum`` must be provided.

    :param obsid: int, str
        Observation ID
    :param propnum: int, str
        Proposal number
    :param timeout: int, float
        Timeout in seconds for the request
    :returns: dict
        Dictionary of proposal abstract
    """
    params = {}
    if obsid is not None:
        params['obsid'] = obsid
    if propnum is not None:
        params['propnum'] = propnum

    if not params:
        raise ValueError('must provide obsid or propnum')

    html = get_cda_service_text('prop_abstract', timeout=timeout, **params)
    text = html_to_text(html)

    delims = ['Proposal Title',
              'Proposal Number',
              'Principal Investigator',
              'Abstract',
              '']
    out = {}
    for delim0, delim1 in zip(delims[:-1], delims[1:]):
        name = '_'.join(word.lower() for word in delim0.split())
        print(rf'{delim0}:(.+){delim1}:')
        if match := re.search(rf'{delim0}:(.+){delim1}:', text, re.DOTALL):
            out[name] = clean_text(match.group(1))
        else:
            warnings.warn(f'failed to find {delim0} in result')

    return out


def get_cda_ocat(query_type='details', timeout=30, **params):
    """
    Get either the Ocat summary or details from the CDA services.

    If ``obsid`` is provided and corresponds to an integer,
    then the returned result is a ``dict``.  Otherwise, the returned result is
    an ``astropy.table.Table``.

    :param query_type: str
        Ocat query type ('details' or 'summary')
    :param **params: dict
        Parameters to filter query
    :param timeout: int, float
        Timeout in seconds for the request
    :return: astropy Table or dict of the observation details
    """
    if query_type not in ['details', 'summary']:
        raise ValueError(f'unknown query_type {query_type!r}, '
                         'must be one of details or summary')

    params['format'] = 'text'
    html = get_cda_service_text(f'ocat_{query_type}', **params)
    dat = get_table_from_cda_rdb_text(html)

    # If obsid is a single integer then return the row as a dict.
    if 'obsid' in params:
        try:
            int(params['obsid'])
        except (ValueError, TypeError):
            pass
        else:
            if len(dat) == 1:
                # Maybe multi-obi obsids early in the mission could have multiple rows?
                dat = dict(dat[0])

    return dat


def get_cda_service_text(service, timeout=30, **params):
    """
    Fetch all observation details from one of the CDA SRService pages

    :param service: str
        Name of the service ('prop_abstract', 'ocat_summary', 'ocat_details', 'archive_file_list')
    :param timeout: int, float
        Timeout in seconds for the request
    :param **params: dict
        Additional parameters to pass to the service
    :return: str
        Returned text from the service
    """

    if service not in CDA_SERVICES:
        raise ValueError(f'unknown service {service!r}, must be one of {list(CDA_SERVICES)}')

    # Query the service and check for errors
    url = f'{URL_CDA_SERVICES}/{CDA_SERVICES[service]}.do'
    resp = requests.get(url, timeout=timeout, params=params)
    if not resp.ok:
        raise RuntimeError(f'got error {resp.status_code} for {resp.url}\n'
                           f'{html_to_text(resp.text)}')

    return resp.text


def get_table_from_cda_rdb_text(text):
    """Get astropy Table from the quasi-RDB text returned by the CDA services.

    :param text: str
        Text returned by the CDA services for a format='text' query
    :return: astropy Table
        Table of the returned data
    """
    lines = text.splitlines()

    # Convert the type line to standard RDB
    # First find the line that begins the column descriptions
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            header_start = i
            break

    # The line with the lengths and types is next (header_start + 1)
    ctypes = lines[header_start + 1].split("\t")

    # Munge length descriptions back to standard RDB (from "20" to "20S" etc)
    # while leaving the numeric types alone ("20N" stays "20N").
    for j, ctype in enumerate(ctypes):
        if not ctype.endswith("N"):
            ctypes[j] = ctype + "S"
    lines[header_start + 1] = "\t".join(ctypes)

    dat = Table.read(lines, format='ascii.rdb', guess=False)

    # Lower-case all the column names
    lc_names = [name.lower() for name in dat.colnames]
    dat.rename_columns(dat.colnames, lc_names)

    return dat


def main_update_ocat_hdf5():
    """
    Command line interface to write an Ocat HDF5 file.

    This overwrites the file from scratch each time.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Update target table")
    parser.add_argument("--datafile",
                        default='target_table.h5')
    opt = parser.parse_args()

    update_ocat_hdf5(opt.datafile)


def update_ocat_hdf5(datafile, **params):
    """Write HDF5 ``datafile`` with the Ocat "details" data.

    :param **params: dict
        Parameters to filter ``get_cda_ocat`` details query
    """
    table = get_cda_ocat(query_type='details', **params)
    table.write(datafile, path='data', serialize_meta=True, overwrite=True, format='hdf5')


def get_ocat_hdf5(datafile=None, read_where=None):
    """
    Read the Ocat target table from an HDF5 data file.

    :param datafile: str
        HDF5 Ocat target table data file.
        Defaults to MICA_ARCHIVE/ocat_target_table.h5
    :param read_where: str
        filter string to pass to tables read_where() to limit returned results.
        See https://www.pytables.org/usersguide/condition_syntax.html

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


