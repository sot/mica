import os

import numpy as np

from astropy.io import fits

__all__ = ['get_dark_cal_image', 'get_dark_cal_properties', 'dark_temp_scale']


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


def get_dark_cal_properties(start=None, stop=None, include_images=False):
    """
    Return a table of dark calibration properties between ``start`` and ``stop``.

    If ``include_images`` is True then an additional column ``dark_image`` is
    defined which contains the corresponding 1024x1024 dark cal image.

    :param start: Start time (default=beginning of mission)
    :param stop: Stop time (default=now)

    :returns: astropy Table of dark calibration properties
    """
    raise NotImplementedError()
