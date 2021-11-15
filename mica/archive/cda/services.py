# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pathlib import Path
import re
import warnings

import requests
import numpy as np
import tables
from astropy.table import Table, MaskedColumn
from astropy.coordinates import SkyCoord

from mica.common import MICA_ARCHIVE


__all__ = ['get_archive_file_list', 'get_proposal_abstract',
           'get_ocat_details_web', 'get_ocat_details_local', 'get_ocat_summary_web']

OCAT_TARGET_TABLE = Path(MICA_ARCHIVE) / 'ocat_target_table.h5'
URL_CDA_SERVICES = "https://cda.harvard.edu/srservices"
CDA_SERVICES = {
    'prop_abstract': 'propAbstract',
    'ocat_summary': 'ocatList',
    'ocat_details': 'ocatDetails',
    'archive_file_list': 'archiveFileList'}

# Units copied from https://github.com/jzuhone/pycda/blob/
# 5a4261328eab989bab91bed17f426ad17d876988/pycda/obscat.py#L38
OCAT_UNITS = {
    "app_exp": "ks",
    "count_rate": "s**-1",
    "est_cnt_rate": "s**-1",
    "evfil_lo": "keV",
    "evfil_ra": "keV",
    "exp_time": "ks",
    "f_time": "s",
    "forder_cnt_rate": "s**-1",
    "soe_roll": "degree",
    "x_sim": "mm",
    "y_off": "arcmin",
    "z_off": "arcmin",
    "z_sim": "mm",
}


RETURN_TYPE_DOCS = """If ``return_type='auto'`` the return type is determined by the rules:

    - If ``obsid`` is provided
      AND the obsid corresponds to an integer
      AND the returned result has a single row
      THEN the return type is ``dict``
      ELSE the return tuple is a ``Table``.

    If ``return_type='table'`` then always return a ``Table``."""

CDA_PARAM_DOCS = """Additional function args for CDA search parameters::

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

        target: matches any part of target name
        piName: matches any part of PI name
        observer: matches any part of observer name
        propNum: proposal number
        propTitle: matches any part of proposal title

    These parameters form a cone search; if you use one you should use them all::

        lon
        lat
        radius (arcmin, default=1.0)

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

COMMON_PARAM_DOCS = """:param target_name: str, optional
        Target name, used in SkyCoord.from_name() to define ``ra`` and ``dec``
        if ``resolve_name`` is True, otherwise matches a substring of the
        table column ``target_name`` (ignoring spaces).
    :param resolve_name: bool, optional
        If True, use ``target_name`` to resolve ``ra`` and ``dec``.
    :param ra: float, optional
        Right Ascension in decimal degrees
    :param dec: float, optional
        Declination in decimal degrees
    :param radius: float, optional
        Search radius in arcmin (default=1.0)"""


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

       >>> get_archive_file_list(obsid=2365, detector='pcad',
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

       >>> get_archive_file_list(obsid=400, detector='acis', level=2, filetype='evt2')
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
    # Original Filesize has commas for the thousands like 11,233,456
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

    # Return value is a text string with these section header lines. Use them
    # to split the text into sections.
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


def _update_params_from_kwargs(params, obsid,
                               target_name, resolve_name,
                               ra, dec, radius):
    """Update params dict for CDA Ocat queries from specified keyword args.
    """
    if obsid is not None:
        params['obsid'] = obsid

    if ra is not None:
        params['ra'] = ra
    if dec is not None:
        params['dec'] = dec

    if target_name is not None:
        if resolve_name:
            coord = SkyCoord.from_name(target_name)
            params['ra'] = coord.ra.deg
            params['dec'] = coord.dec.deg
        else:
            # SR services API uses "target" to select substrings of target_name
            params['target'] = target_name

    # For any positional search include the radius
    if 'ra' in params and 'dec' in params:
        params['radius'] = radius

    return params


def get_ocat_summary_web(obsid=None, *,
                         target_name=None, resolve_name=False,
                         ra=None, dec=None, radius=1.0,
                         return_type='auto',
                         timeout=30, **params):
    """
    Get the Ocat summary from the CDA services.

    {RETURN_TYPE_DOCS}

    {CDA_PARAM_DOCS}

    :param obsid: int, str
        Observation ID or string with ObsId range or list of ObsIds
    {COMMON_PARAM_DOCS}
    :param timeout: int, float
        Timeout in seconds for the request
    :param return_type: str
        Return type (default='auto' => Table or dict)
    :param **params: dict
        Parameters passed to CDA web service
    :return: astropy Table or dict of the observation details
    """
    _update_params_from_kwargs(params, obsid,
                               target_name, resolve_name,
                               ra, dec, radius)
    dat = _get_ocat_web('ocat_summary', timeout, return_type, **params)
    return dat


get_ocat_summary_web.__doc__ = get_ocat_summary_web.__doc__.format(
    RETURN_TYPE_DOCS=RETURN_TYPE_DOCS,
    CDA_PARAM_DOCS=CDA_PARAM_DOCS,
    COMMON_PARAM_DOCS=COMMON_PARAM_DOCS)


def get_ocat_details_web(obsid=None, *,
                         target_name=None, resolve_name=False,
                         ra=None, dec=None, radius=1.0,
                         return_type='auto',
                         timeout=30, **params):
    """
    Get the Ocat details from the CDA services.

    {RETURN_TYPE_DOCS}

    {CDA_PARAM_DOCS}

    Special parameters that change the output table contents:
    - ``acisWindows='true'``: return ACIS windows details for a single obsid
    - ``rollReqs='true'``: return roll requirements for a single obsid
    - ``timeReqs='true'``: return time requirements for a single obsid

    :param obsid: int, str
        Observation ID or string with ObsId range or list of ObsIds
    {COMMON_PARAM_DOCS}
    :param timeout: int, float
        Timeout in seconds for the request
    :param return_type: str
        Return type (default='auto' => Table or dict)
    :param **params: dict
        Parameters passed to CDA web service
    :return: astropy Table or dict of the observation details
    """
    _update_params_from_kwargs(params, obsid,
                               target_name, resolve_name,
                               ra, dec, radius)
    dat = _get_ocat_web('ocat_details', timeout, return_type, **params)
    return dat


get_ocat_details_web.__doc__ = get_ocat_details_web.__doc__.format(
    RETURN_TYPE_DOCS=RETURN_TYPE_DOCS,
    CDA_PARAM_DOCS=CDA_PARAM_DOCS,
    COMMON_PARAM_DOCS=COMMON_PARAM_DOCS)


def _get_ocat_web(service, timeout=30, return_type='auto', **params):
    """
    Get either the Ocat summary or details from the CDA services.

    :param service: str
        Service name
    :param timeout: int, float
        Timeout in seconds for the request
    :param **params: dict
        Parameters passed to CDA web service
    :return: astropy Table or dict of the observation details
    """
    if return_type not in ('auto', 'table'):
        raise ValueError(f"invalid return_type {return_type!r}, must be 'auto' or 'table'")

    params['format'] = 'text'

    # Force RA, Dec in sexagesimal because decimal returns only 3 decimal digits
    # which is insufficient.
    params['outputCoordUnits'] = 'sexagesimal'

    text = _get_cda_service_text(service, **params)
    dat = _get_table_or_dict_from_cda_rdb_text(text, return_type, params.get('obsid'))

    if dat is None:
        # Query returned no rows. If a single obsid was specified with return_type
        # of 'auto' then we would have expected to return a dict, but instead
        # raise a ValueError. Otherwise we return an empty table with the right
        # column names.
        if return_type == 'auto' and _is_int(params.get('obsid')):
            raise ValueError(f"failed to find obsid {params['obsid']}")
        else:
            dat = _get_ocat_web(service, return_type='table', obsid=8000)[0:0]

    # Change RA, Dec to decimal
    sc = SkyCoord(dat['ra'], dat['dec'], unit='hr,deg')
    dat['ra'] = sc.ra.deg
    dat['dec'] = sc.dec.deg

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
    verbose = params.pop('verbose', False)
    resp = requests.get(url, timeout=timeout, params=params)
    if verbose:
        print(f'GET {resp.url}')

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
    :return: astropy Table, dict, None
        Table of the returned data, or dict if just one obsid selected, or None
        if the query returned no data.
    """
    lines = text.splitlines()

    # Convert the type line to standard RDB
    # First find the line that begins the column descriptions
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            header_start = i
            break
    else:
        return None

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

    # Apply units to the columns
    for name, col in dat.columns.items():
        if name in OCAT_UNITS:
            col.info.unit = OCAT_UNITS[name]

    # Possibly get just the first row as a dict
    dat = _get_table_or_dict(return_type, obsid, dat)

    return dat


def _is_int(val):
    """Check if a value looks like an integet."""
    try:
        return int(val) == float(val)
    except (ValueError, TypeError):
        return False


def _get_table_or_dict(return_type, obsid, dat):
    # If obsid is a single integer and there was just one row then return the
    # row as a dict.
    if (return_type == 'auto'
            and _is_int(obsid)
            and len(dat) == 1):
        dat = dict(dat[0])
    return dat


def main_update_ocat_local():
    """
    Command line interface to write a local Ocat HDF5 file.

    This overwrites the file from scratch each time.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Update target table")
    parser.add_argument("--datafile",
                        default='target_table.h5')
    opt = parser.parse_args()

    update_ocat_local(opt.datafile)


def update_ocat_local(datafile, **params):
    """Write HDF5 ``datafile`` with the Ocat "details" data.

    :param **params: dict
        Parameters to filter ``get_ocat_details_web`` query
    """
    dat = get_ocat_details_web(**params)

    # Encode unicode strings to bytes manually.  Fixed in numpy 1.20.
    # Eventually we will want just dat.convert_bytestring_to_unicode().
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'U':
            dat[name] = np.char.encode(col, 'utf-8')

    dat.write(datafile, path='data', serialize_meta=True, overwrite=True, format='hdf5')


def get_ocat_details_local(obsid=None, *,
                           target_name=None, resolve_name=False,
                           ra=None, dec=None, radius=1.0,
                           return_type='auto',
                           datafile=None, where=None, **params):
    """
    Read the Ocat target table from a local data file.

    The local data file is assumed to be an HDF5 file that contains a copy of
    the Ocat details, typically updated by a cron job running on HEAD and
    potentially synced to the local host.

    {RETURN_TYPE_DOCS}

    :param obsid: int, optional
        Observation ID
    {COMMON_PARAM_DOCS}
    :param datafile: str, optional
        HDF5 Ocat target table data file.
        Defaults to MICA_ARCHIVE/ocat_target_table.h5
    :param where: str
        Filter string to pass to tables read_where() to limit returned results.
        See https://www.pytables.org/usersguide/condition_syntax.html
    :param **params: dict
        Additional filter criteria as ``<colname> == <value>`` key/value pairs.

    :returns: astropy Table or dict of Ocat details
    """
    where_parts = []  # build up bits of the where clause
    if where is not None:
        where_parts.append(where)

    if datafile is None:
        datafile = OCAT_TARGET_TABLE

    if obsid is not None:
        where_parts.append(f"obsid=={obsid}")

    if target_name is not None and resolve_name:
        coord = SkyCoord.from_name(target_name)
        ra = coord.ra.deg
        dec = coord.dec.deg

    if ra is not None and dec is not None:
        d2r = np.pi / 180.0  # Degrees to radians
        # Use great-circle distance to find targets within radius. This is
        # accurate enough for this application.
        where = (f'arccos(sin({ra * d2r})*sin(ra*{d2r}) + '
                 f'cos({ra * d2r})*cos(ra*{d2r})*cos({dec*d2r}-dec*d2r))'
                 f'< {radius / 60 * d2r}')
        where_parts.append(where)

    for col_name, value in params.items():
        where_parts.append(f'{col_name}=={value!r}')

    if where_parts:
        dat = _table_read_where(datafile, where_parts)
    else:
        dat = Table.read(datafile)

    # Decode bytes to strings manually.  Fixed in numpy 1.20.
    # Eventually we will want just dat.convert_bytestring_to_unicode().
    for name, col in dat.columns.items():
        if col.info.dtype.kind == 'S':
            dat[name] = np.char.decode(col, 'utf-8')

    # Match target_name as a substring of the table target_name column.
    if target_name is not None and not resolve_name:
        target_name = target_name.lower().replace(' ', '')
        target_names = np.char.lower(np.char.replace(dat['target_name'], ' ', ''))
        ok = np.char.find(target_names, target_name) != -1
        dat = dat[ok]

    # Apply units to the columns
    for name, col in dat.columns.items():
        if name in OCAT_UNITS:
            col.info.unit = OCAT_UNITS[name]

    # Possibly get just the first row as a dict
    dat = _get_table_or_dict(return_type, obsid, dat)

    return dat


def _table_read_where(datafile, where_parts):
    """Read HDF5 ``datafile`` using read_where() and ``where_parts``.
    """
    where = '&'.join(f'({where})' for where in where_parts)

    with tables.open_file(datafile) as h5:
        # PyTables is unhappy with all the column names that cannot be an
        # object attribute, so squelch that warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=tables.NaturalNameWarning)
            dat = h5.root.data.read_where(where)
    dat = Table(dat)

    # Manually make MaskedColumn's as needed since we are not using Table.read()
    # which handles this. This assumes a correspondence betweeen <name>.mask
    # and <name>, but this is always true for the Ocat target table.
    masked_names = [name for name in dat.colnames if name.endswith('.mask')]
    for masked_name in masked_names:
        name = masked_name[:-5]
        dat[name] = MaskedColumn(dat[name], mask=dat[masked_name])
        dat.remove_column(masked_name)
    return dat
