# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pathlib import Path
import re
import warnings

import requests
import numpy as np
import tables
from astropy.table import Table

from mica.common import MICA_ARCHIVE


__all__ = ['get_archive_file_list', 'get_proposal_abstract',
           'get_ocat_details', 'get_ocat_details_local', 'get_ocat_summary']

OCAT_TARGET_TABLE = Path(MICA_ARCHIVE) / 'ocat_target_table.h5'
URL_CDA_SERVICES = "https://cda.harvard.edu/srservices"
CDA_SERVICES = {
    'prop_abstract': 'propAbstract',
    'ocat_summary': 'ocatList',
    'ocat_details': 'ocatDetails',
    'archive_file_list': 'archiveFileList'}

PARAMETER_DOCS = """
    Search parameters::

        instrument=ACIS,ACIS-I,ACIS-S,HRC,HRC-I,HRC-S
        grating=NONE,LETG,HETG
        type=ER,TOO,CAL,GO,GTO,DDT
        cycle=00,01,02,03,04, ...
        category=SOLAR SYSTEM,
            NORMAL GALAXIES,
            STARS AND WD,
            WD BINARIES AND CV,
            BH AND NS BINARIES,
            NORMAL GALAXIES
            CLUSTERS OF GALAXIES,
            ACTIVE GALAXIES AND QUASARS,
            GALACTIC DIFFUSE EMISSION AND SURVEYS,
            EXTRAGALACTIC DIFFUSE EMISSION AND SURVEYS
        jointObservation= HST,XMM,Spitzer,NOAO,NRAO,NuSTAR,Suzaku,Swift,RXTE
        status= archived,observed,scheduled, unobserved,untriggered,canceled,deleted
        expMode= ACIS TE,ACIS CC,HRC Timing
        grid = 'is not null' or 'is null'

    Input coordinate specifications::

        inputCoordFrame=J2000 (other options: b1950, bxxxx, ec1950, ecxxxx, galactic)
        inputCoordEquinox=2000 (4 digit year)

    These parameters are single text entries::

        target: name will be resolved to RA, Dec
        piName: matches any part of PI name
        observer: matches any part of observer name
        propNum: proposal number
        propTitle: matches any part of proposal title

    These parameters form a cone search; if you use one you should use them all::

        lon
        lat
        radius (arcsec)

    These parameters form a box search; one lon & one lat are required.
    Open-ended ranges are allowed. (Such as lonMin=15 with no lonMax.)
    ::

        lonMin
        lonMax
        latMin
        latMax

    These parameters are range lists, where the range is indicated by a hyphen (-).
    Multiple ranges can be entered separated by commas::

        obsid  (eg. obsid=100,200-300,600-1000,1800)
        seqNum
        expTime
        appExpTime
        countRate

    These parameters are date range lists, where the range is
    indicated by a hyphen (/). Multiple ranges can be entered separated
    by commas. Valid dates are in one of the following formats:
    YYYY-MM-DD, YYYY-MM-DD hh:mm, or YYYY-MM-DD hh:mm:ss
    ::

        startDate
        releaseDate

    These specify how the data is displayed and ordered::

        outputCoordFrame=J2000 (other options: b1950, bxxxx, ec1950, ecxxxx, galactic)
        outputCoordEquinox=2000 (4 digit year)
        outputCoordUnits=decimal (other option: sexagesimal)
        sortColumn=ra (other options:
                dec,seqNum,instrument,grating,
                appExpTime,expTime,
                target,piName,observer,status,
                startDate,releaseDate,
                obsid,propNum,category,type,cycle)
        sortOrder=ascending (other option: descending)
        maxResults=#  (the number of results after which to stop displaying)
"""


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


def get_archive_file_list(obsid, detector, level, dataset='flight', **params):
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

    text = _get_cda_service_text('archive_file_list', **params)
    dat = Table.read(text.splitlines(), format='ascii.basic', delimiter='\t', guess=False)
    filesize = [int(x.replace(',', '')) for x in dat['Filesize']]
    dat['Filesize'] = filesize

    return dat


def get_proposal_abstract(obsid=None, propnum=None, timeout=30):
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

    html = _get_cda_service_text('prop_abstract', timeout=timeout, **params)
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


def get_ocat_summary(timeout=30, return_type='auto', **params):
    """
    Get the Ocat summaryfrom the CDA services.

    If ``return_type`` is 'auto', the return type is determined by the rules:
    - If ``obsid`` is provided AND the obsid corresponds to an integer
      AND the returned result has a single row THEN the return type is ``dict``.
    - Otherwise the return tuple is a ``Table``.
    Specify ``return_type='table'`` to always return a ``Table``.

    {PARAMETER_DOCS}

    :param timeout: int, float
        Timeout in seconds for the request
    :param return_type: str
        Return type (default='auto' => Table or dict)
    :param **params: dict
        Parameters passed to CDA web service
    :return: astropy Table or dict of the observation details
    """
    dat = _get_ocat_cda('ocat_summary', timeout, return_type, **params)
    return dat


get_ocat_summary.__doc__ = get_ocat_summary.__doc__.format(PARAMETER_DOCS=PARAMETER_DOCS)


def get_ocat_details(timeout=30, return_type='auto', **params):
    """
    Get the Ocat details from the CDA services.

    If ``return_type`` is 'auto', the return type is determined by the rules::

      If ``obsid`` is provided
        AND the obsid corresponds to an integer
        AND the returned result has a single row
      THEN the return type is ``dict``.
      ELSE the return tuple is a ``Table``.

    Specify ``return_type='table'`` to always return a ``Table``.

    {PARAMETER_DOCS}

    Special parameters that change the output table contents:
    - ``acisWindows='true'``: return ACIS windows details for a single obsid
    - ``rollReqs='true'``: return roll requirements for a single obsid
    - ``timeReqs='true'``: return time requirements for a single obsid

    :param timeout: int, float
        Timeout in seconds for the request
    :param return_type: str
        Return type (default='auto' => Table or dict)
    :param **params: dict
        Parameters passed to CDA web service
    :return: astropy Table or dict of the observation details
    """
    dat = _get_ocat_cda('ocat_details', timeout, return_type, **params)
    return dat


get_ocat_details.__doc__ = get_ocat_details.__doc__.format(PARAMETER_DOCS=PARAMETER_DOCS)


def _get_ocat_cda(service, timeout=30, return_type='auto', **params):
    """
    Get either the Ocat summary or details from the CDA services.

    :param timeout: int, float
        Timeout in seconds for the request
    :param **params: dict
        Parameters passed to CDA web service
    :return: astropy Table or dict of the observation details
    """
    if return_type not in ('auto', 'table'):
        raise ValueError(f"invalid return_type {return_type!r}, must be 'auto' or 'table'")

    params['format'] = 'text'

    # Default to decimal RA, Dec not sexagesimal
    if 'outputCoordUnits' not in params:
        params['outputCoordUnits'] = 'decimal'

    html = _get_cda_service_text(service, **params)
    dat = _get_table_or_dict_from_cda_rdb_text(html, return_type, params.get('obsid'))

    return dat


def _get_cda_service_text(service, timeout=30, **params):
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


def _get_table_or_dict_from_cda_rdb_text(text, return_type, obsid):
    """Get astropy Table or dict from the quasi-RDB text returned by the CDA services.

    :param text: str
        Text returned by the CDA services for a format='text' query
    :param return_type: str
        Return type (default='auto' => Table or dict)
    :param obsid: int, str, None
        Observation ID if provided
    :return: astropy Table, dict
        Table of the returned data, or dict if just one obsid selected
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

    # If obsid is a single integer then return the row as a dict.
    if return_type == 'auto' and obsid is not None:
        try:
            int(obsid)
        except (ValueError, TypeError):
            pass
        else:
            if len(dat) == 1:
                # Maybe multi-obi obsids early in the mission could have multiple rows?
                dat = dict(dat[0])

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
    dat = get_ocat_details(**params)

    # Encode unicode strings to bytes manually.  Fixed in numpy 1.20.
    # Eventually we will want just dat.convert_bytestring_to_unicode().
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'U':
            dat[name] = np.char.encode(col, 'utf-8')

    dat.write(datafile, path='data', serialize_meta=True, overwrite=True, format='hdf5')


def get_ocat_details_local(datafile=None, where=None):
    """
    Read the Ocat target table from a local data file.

    The local data file is assumed to be an HDF5 file that contains a copy of
    the Ocat details, typically updated by a cron job running on HEAD and
    potentially synced to the local host.

    :param datafile: str
        HDF5 Ocat target table data file.
        Defaults to MICA_ARCHIVE/ocat_target_table.h5
    :param where: str
        Filter string to pass to tables read_where() to limit returned results.
        See https://www.pytables.org/usersguide/condition_syntax.html

    :returns: astropy table of target table
    """
    if datafile is None:
        datafile = OCAT_TARGET_TABLE

    if where is None:
        dat = Table.read(datafile)
    else:
        with tables.open_file(datafile) as hdu:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=tables.NaturalNameWarning)
                dat = hdu.root.data.read_where(where)
            dat = Table(dat)

    # Decode bytes to strings manually.  Fixed in numpy 1.20.
    # Eventually we will want just dat.convert_bytestring_to_unicode().
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'S':
            dat[name] = np.char.decode(col, 'utf-8')

    return dat
