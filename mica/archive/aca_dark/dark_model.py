"""
Routines related to the dark current model and guide / acq success prediction.
"""


def get_dark_model(date, t_ccd):
    """
    Return the dark current model corresponding to ``date`` and ``t_ccd``.

    :param date: date in any DateTime format
    :param t_ccd: CCD temperature (deg C)

    :returns: TBD
    """
    raise NotImplementedError()


def get_acq_success(date, t_ccd, mag):
    """
    Return probability of acquisition success for given date, temperature and mag.

    Any of the inputs can be scalars or arrays, with the output being the result of
    the broadcasted dimension of the inputs.

    This is based on the dark model and acquisition success fitting presented
    in the State of the ACA 2013 (sot/state_of_aca/guide_acq_stats)

    :param date: Date(s) (scalar or np.ndarray)
    :param t_ccd: CD temperature(s) (degC, scalar or np.ndarray)
    :param mag: Star magnitude(s) (scalar or np.ndarray)

    :returns: Acquisition success probability(s)
    """
    raise NotImplementedError()


def get_guide_success(date, t_ccd, mag):
    """
    Return probability of guide (bad_trak) success for given date, temperature and mag.

    Any of the inputs can be scalars or arrays, with the output being the result of
    the broadcasted dimension of the inputs.

    This is based on the dark model and guide success fitting presented
    in the State of the ACA 2013 (sot/state_of_aca/guide_acq_stats)

    :param date: Date(s) (scalar or np.ndarray)
    :param t_ccd: CD temperature(s) (degC, scalar or np.ndarray)
    :param mag: Star magnitude(s) (scalar or np.ndarray)

    :returns: Guide success probability(s)
    """
    raise NotImplementedError()
