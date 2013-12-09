__all__ = ['get_dark_cal_image', 'get_dark_cal_properties']


def get_dark_cal_image(date, before_date=False):
    """
    Return the dark calibration image (e-/s) nearest to ``date``.

    If ``before_date`` is True (default) then use the first calibration which
    occurs before ``date``.

    :param date: date in any DateTime format
    :param before_date: use first cal before date

    :returns: 1024 x 1024 ndarray with dark cal image in e-/s
    """
    raise NotImplementedError()


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
