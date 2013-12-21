import re
import os
from collections import OrderedDict
import json

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
import pyyaks.context
from Chandra.Time import DateTime

from mica.common import MICA_ARCHIVE_PATH
from . import file_defs

DARK_CAL = pyyaks.context.ContextDict('dark_cal')

SKA_FILES = pyyaks.context.ContextDict('ska_files', basedir='/proj/sot/ska')
SKA_FILES.update(file_defs.SKA_FILES)

MICA_FILES = pyyaks.context.ContextDict('update_mica_files',
                                        basedir=os.path.join(MICA_ARCHIVE_PATH))
MICA_FILES.update(file_defs.MICA_FILES)

__all__ = ['get_dark_cal_image', 'get_dark_cal_properties', 'dark_temp_scale',
           'get_dark_cal_dirs']


def date_to_dark_id(date):
    date = DateTime(date).date
    return date[:4] + date[5:8]


def dark_id_to_date(dark_id):
    return '{}:{}'.format(dark_id[:4], dark_id[4:])


def dark_temp_scale(t_ccd, t_ccd_ref=-19.0):
    """
    Return the multiplicative scale factor to convert a CCD dark map from
    the actual temperature ``t_ccd` to the reference temperature ``t_ccd_ref``.

    Based on best global fit for dark current model in plot_predicted_warmpix.py.
    Previous value was 0.62 instead of 0.70.  This represents the change in
    dark current for each 4 degC decrease::

      >>> from mica.aca_dark import temp_scalefac
      >>> print temp_scalefac(t_ccd=-15, t_ccd_ref=-19)
      0.7

    :param t_ccd: actual temperature (degC)
    :param t_ccd_ref: reference temperature (degC, default=-19.0)

    :returns: scale factor
    """
    return np.exp(np.log(0.70) / 4.0 * (t_ccd - t_ccd_ref))


def get_dark_cal_dirs(source='mica'):
    """
    Get an ordered dict of directory paths containing dark current calibration files,
    where the key is the dark cal identifier (YYYYDOY) and the value is the path.

    :param source: source of dark cal directories ('mica'|'ska')
    :returns: ordered dict of absolute directory paths
    """
    files = {'ska': SKA_FILES, 'mica': MICA_FILES}[source]
    dark_cal_ids = sorted([fn for fn in os.listdir(files['dark_cals_dir'].abs)
                           if re.match(r'[12]\d{6}$', fn)])
    dark_cal_dirs = [os.path.join(files['dark_cals_dir'].abs, id_)
                     for id_ in dark_cal_ids]
    return OrderedDict(zip(dark_cal_ids, dark_cal_dirs))


def get_dark_cal_image(date, before_date=False, t_ccd_ref=None):
    """
    Return the dark calibration image (e-/s) nearest to ``date``.

    If ``before_date`` is True (default) then use the first calibration which
    occurs before ``date``.

    :param date: date in any DateTime format
    :param before_date: use first cal before date
    :param t_ccd_ref: rescale dark map to temperature (degC, default=no scaling)

    :returns: 1024 x 1024 ndarray with dark cal image in e-/s
    """

    hdus = fits.open(os.path.join('proj', 'sot', 'ska', 'data', 'aca_dark_cal', date, 'imd.fits'))
    dark = hdus[0].data
    hdus.close()

    if t_ccd_ref is not None:
        # Scale factor to adjust data to an effective temperature of t_ccd_ref.
        # For t_ccds warmer than t_ccd_ref this scale factor is < 1, i.e. the
        # observed dark current is made smaller to match what it would be at the
        # lower reference temperature.

        # t_ccd = ???
        # dark *= dark_temp_scale(t_ccd, t_ccd_ref)
        raise NotImplementedError('No code yet to get dark calibration temperature')

    return dark


def get_dark_cal_properties(start=None, stop=None, include_images=False, as_table=True):
    """
    Return a table of dark calibration properties between ``start`` and ``stop``.

    If ``include_images`` is True then an additional column or key ``image`` is
    defined which contains the corresponding 1024x1024 dark cal image.

    If ``as_table`` is True (default) then the result is an astropy Table object.
    If False then a list of dicts is returned.  In this case the full contents
    of the properties file including replica properties is available.

    :param start: start time (default=beginning of mission)
    :param stop: stop time (default=now)
    :param include_images: include the dark cal images in output (default=False)
    :param as_table: return a Table instead of a list (default=True)

    :returns: astropy Table or list of dark calibration properties
    """
    # TODO: use a context manager or decorator to cache DARK_CAL
    start_id = date_to_dark_id('1999:001' if start is None else start)
    stop_id = date_to_dark_id(stop)
    dark_dirs = [dark_id for dark_id in get_dark_cal_dirs()
                 if dark_id >= start_id and dark_id <= stop_id]
    props = []
    for dark_id in dark_dirs:
        DARK_CAL['id'] = dark_id
        with open(MICA_FILES['dark_props.json'].abs, 'r') as fh:
            prop = json.load(fh)
            # Change unicode to ascii at top level
            keys = prop.keys()
            asciiprop = {}
            for key in keys:
                asciiprop[str(key)] = prop[key]
            prop = asciiprop

        if include_images:
            hdus = fits.open(MICA_FILES['dark_image.fits'].abs)
            prop['image'] = hdus[0].data
            hdus.close()

        props.append(prop)

    if as_table:
        for prop in props:
            del prop['replicas']
        table_props = Table(props)
        if include_images:
            x = np.vstack([prop['image'][np.newaxis, :] for prop in props])
            images = Column(x)
            table_props['image'] = images
        props = table_props

    return props
