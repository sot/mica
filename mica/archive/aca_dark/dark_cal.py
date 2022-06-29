# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import os
from collections import OrderedDict
import json

import six
from six.moves import zip

import numpy as np
from astropy.io import fits, ascii
from astropy.table import Column
import pyyaks.context
from cxotime import CxoTime
from chandra_aca.aca_image import ACAImage

from chandra_aca.dark_model import dark_temp_scale
from chandra_aca.dark_model import DARK_SCALE_4C  # noqa
from mica.cache import lru_cache
from mica.common import MICA_ARCHIVE_PATH, MissingDataError
from . import file_defs

DARK_CAL = pyyaks.context.ContextDict('dark_cal')

MICA_FILES = pyyaks.context.ContextDict('update_mica_files',
                                        basedir=MICA_ARCHIVE_PATH)
MICA_FILES.update(file_defs.MICA_FILES)


def date_to_dark_id(date):
    """
    Convert ``date`` to the corresponding YYYYDOY format for a dark cal identifiers.

    :param date: any CxoTime compatible format
    :returns: dark id (YYYYDOY)
    """
    time = CxoTime(date)
    date_str = time.date
    if not time.shape:
        return date_str[:4] + date_str[5:8]
    date_str = np.atleast_1d(date_str)
    chars = date_str.view((str, 1)).reshape(-1, date_str.dtype.itemsize // 4)
    result = np.hstack([chars[:, :4], chars[:, 5:8]])
    result = np.frombuffer(result.tobytes(), dtype=(str, 7))
    return result.reshape(time.shape)


def dark_id_to_date(dark_id):
    """
    Convert ``dark_id`` (YYYYDOY) to the corresponding CxoTime 'date' format.

    :param date: dark id (YYYYDOY)
    :returns: str in CxoTime 'date' format
    """
    shape = np.shape(dark_id)
    if not shape:
        return '{}:{}'.format(dark_id[:4], dark_id[4:])
    b = np.atleast_1d(dark_id).view((str, 1)).reshape(-1, 7)
    sep = np.array([':'] * len(b))[:, None]
    b = np.hstack([b[:, :4], sep, b[:, 4:]])
    b = np.frombuffer(b.tobytes(), dtype=(str, 8))
    return b.reshape(shape)


@lru_cache()
def get_dark_cal_dirs(dark_cals_dir=MICA_FILES['dark_cals_dir'].abs):
    """
    Get an ordered dict of directory paths containing dark current calibration files,
    where the key is the dark cal identifier (YYYYDOY) and the value is the path.

    :param source: source of dark cal directories ('mica'|'ska')
    :returns: ordered dict of absolute directory paths
    """
    dark_cal_ids = sorted([fn for fn in os.listdir(dark_cals_dir)
                           if re.match(r'[12]\d{6}$', fn)])
    dark_cal_dirs = [os.path.join(dark_cals_dir, id_)
                     for id_ in dark_cal_ids]
    return OrderedDict(zip(dark_cal_ids, dark_cal_dirs))


def get_dark_cal_id(date, select='before'):
    """
    Return the dark calibration id corresponding to ``date``.

    If ``select`` is ``'before'`` (default) then use the first calibration which
    occurs before ``date``.  Other valid options are ``'after'`` and ``'nearest'``.

    :param date: date in any CxoTime format
    :param select: method to select dark cal (before|nearest|after)

    :returns: dark cal id string (YYYYDOY)
    """
    return _get_dark_cal_id_vector(date, select=select)


def _get_dark_cal_id_scalar(date, select='before'):
    dark_cals = get_dark_cal_dirs()
    dark_id = date_to_dark_id(date)

    # Special case if dark_id is exactly an existing dark cal then return that dark
    # cal regardless of the select method.
    if dark_id in dark_cals:
        return dark_id

    dark_cal_ids = list(dark_cals.keys())
    date_secs = CxoTime(date).secs
    dark_cal_secs = CxoTime(np.array([dark_id_to_date(id_) for id_ in dark_cal_ids])).secs

    if select == 'nearest':
        ii = np.argmin(np.abs(dark_cal_secs - date_secs))
    elif select in ('before', 'after'):
        ii = np.searchsorted(dark_cal_secs, date_secs)
        if select == 'before':
            ii -= 1
    else:
        raise ValueError('select arg must be one of "nearest", "before", or "after"')

    if ii < 0:
        earliest = CxoTime(dark_cal_secs[0]).date[:8]
        raise ValueError(
            f'trying to get last dark cal before {date}, '
            f'which is before the earliest dark cal on {earliest}'
        )

    try:
        out_dark_id = dark_cal_ids[ii]
    except IndexError:
        raise MissingDataError('No dark cal found {} {}'.format(select, date))

    return out_dark_id


_get_dark_cal_id_vector = np.vectorize(_get_dark_cal_id_scalar, excluded=['select'])


@DARK_CAL.cache
def _get_dark_cal_image_props(date, select='before', t_ccd_ref=None, aca_image=False,
                              allow_negative=False):
    """
    Return the dark calibration image (e-/s) nearest to ``date`` and the corresponding
    dark_props file.

    :param date: date in any CxoTime format
    :param select: method to select dark cal (before|nearest|after)
    :param t_ccd_ref: rescale dark map to temperature (degC, default=no scaling)
    :param aca_image: return an AcaImage instance
    :param allow_negative: allow negative values in raw dark map (default=False)

    :returns: 1024 x 1024 ndarray with dark cal image in e-/s, props dict
    """
    DARK_CAL['id'] = get_dark_cal_id(date, select)

    with fits.open(MICA_FILES['dark_image.fits'].abs, memmap=False) as hdus:
        dark = hdus[0].data
        # Recast as native byte ordering since FITS is typically not. This
        # statement is normally the same as dark.astype(np.float32).
        dark = dark.astype(dark.dtype.type)

    with open(MICA_FILES['dark_props.json'].abs, 'r') as fh:
        props = json.load(fh)

    # Change unicode to ascii at top level
    if six.PY2:
        keys = list(props.keys())
        asciiprops = {}
        for key in keys:
            asciiprops[str(key)] = props[key]
        props = asciiprops

    if t_ccd_ref is not None:
        # Scale factor to adjust data to an effective temperature of t_ccd_ref.
        # For t_ccds warmer than t_ccd_ref this scale factor is < 1, i.e. the
        # observed dark current is made smaller to match what it would be at the
        # lower reference temperature.
        t_ccd = props['ccd_temp']
        dark *= dark_temp_scale(t_ccd, t_ccd_ref)

    if not allow_negative:
        np.clip(dark, a_min=0, a_max=None, out=dark)

    if aca_image:
        dark = ACAImage(dark, row0=-512, col0=-512)

    return dark, props


def get_dark_cal_image(date, select='before', t_ccd_ref=None, aca_image=False,
                       allow_negative=False):
    """
    Return the dark calibration image (e-/s) nearest to ``date``.

    If ``select`` is ``'before'`` (default) then use the first calibration which
    occurs before ``date``.  Other valid options are ``'after'`` and ``'nearest'``.

    :param date: date in any CxoTime format
    :param select: method to select dark cal (before|nearest|after)
    :param t_ccd_ref: rescale dark map to temperature (degC, default=no scaling)
    :param aca_image: return an ACAImage instance instead of ndarray
    :param allow_negative: allow negative values in raw dark map (default=False)

    :returns: 1024 x 1024 ndarray with dark cal image in e-/s
    """
    dark, props = _get_dark_cal_image_props(date, select=select, t_ccd_ref=t_ccd_ref,
                                            aca_image=aca_image,
                                            allow_negative=allow_negative)
    return dark


def get_dark_cal_props(date, select='before', include_image=False, t_ccd_ref=None,
                       aca_image=False, allow_negative=False):
    """
    Return a dark calibration properties structure for ``date``

    If ``select`` is ``'before'`` (default) then use the first calibration which
    occurs before ``date``.  Other valid options are ``'after'`` and ``'nearest'``.

    If ``include_image`` is True then an additional column or key ``image`` is
    defined which contains the corresponding 1024x1024 dark cal image.

    :param date: date in any CxoTime format
    :param select: method to select dark cal (before|nearest|after)
    :param include_image: include the dark cal images in output (default=False)
    :param t_ccd_ref: rescale dark map to temperature (degC, default=no scaling)
    :param aca_image: return an ACAImage instance instead of ndarray
    :param allow_negative: allow negative values in raw dark map (default=False)

    :returns: dict of dark calibration properties
    """
    dark, props = _get_dark_cal_image_props(date, select=select, t_ccd_ref=None,
                                            aca_image=aca_image,
                                            allow_negative=allow_negative)

    if include_image:
        props['image'] = dark

    return props


def get_dark_cal_props_table(start=None, stop=None, include_image=False, as_table=True):
    """
    Return a table of dark calibration properties between ``start`` and ``stop``.

    If ``include_image`` is True then an additional column or key ``image`` is
    defined which contains the corresponding 1024x1024 dark cal image.

    If ``as_table`` is True (default) then the result is an astropy Table object.
    If False then a list of dicts is returned.  In this case the full contents
    of the properties file including replica properties is available.

    :param start: start time (default=beginning of mission)
    :param stop: stop time (default=now)
    :param include_image: include the dark cal images in output (default=False)
    :param as_table: return a Table instead of a list (default=True)

    :returns: astropy Table or list of dark calibration properties
    """
    start_id = date_to_dark_id('1999:001:12:00:00' if start is None else start)
    stop_id = date_to_dark_id(stop)
    dark_dirs = [dark_id for dark_id in get_dark_cal_dirs()
                 if dark_id >= start_id and dark_id <= stop_id]

    # Get the list of properties structures
    props = [get_dark_cal_props(dark_id, include_image=include_image) for dark_id in dark_dirs]

    if as_table:
        # Get rid of a non-scalar data structures and collect col names
        names = []
        for prop in props:
            for key, val in prop.copy().items():
                if isinstance(val, (list, tuple, dict)):
                    del prop[key]
                elif key not in names:
                    names.append(key)

        lines = []
        lines.append(','.join(names))
        for prop in props:
            vals = [str(prop.get(name, '')) for name in names]
            lines.append(','.join(vals))
        table_props = ascii.read(lines, format='csv')

        if include_image:
            x = np.vstack([prop['image'][np.newaxis, :] for prop in props])
            images = Column(x)
            table_props['image'] = images
        props = table_props

    return props
