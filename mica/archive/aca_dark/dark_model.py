"""
Routines related to the dark current model and guide / acq success prediction.
"""

from itertools import izip
import numpy as np
from numpy import exp, log, arange

import Ska.Numpy
from Chandra.Time import DateTime

# Define a common fixed binning of dark current distribution
import darkbins

# Some constants and globals.  Done this way to support sherpa fitting.
# Needs to be re-worked to be nicer.

# Fixed gaussian for smoothing the broken power law
dx = 0.1
sigma = 0.30                            # Gaussian sigma in log space
xg = arange(-2.5 * sigma, 2.5 * sigma, dx, dtype=float)
yg = exp(-0.5 * (xg / sigma) ** 2)
yg /= np.sum(yg)

NPIX = 1024 ** 2

# Fixed
xbins = darkbins.bins
xall = darkbins.bin_centers
imin = 0
imax = len(xall)

# scale and offset fit of polynomial to acq failures in log space
acq_fit = {
    'scale': (-0.491, 0.990, 0.185),
    'offset': (0.280, 0.999, -1.489),
    }

warm_threshold = 100.

def get_dark_model(date, t_ccd):
    """
    Return the dark current model corresponding to ``date`` and ``t_ccd``.

    :param date: date in any DateTime format
    :param t_ccd: CCD temperature (deg C)

    :returns: TBD
    """
    raise NotImplementedError()


def get_dark_hist(date, t_ccd):
    """
    Return the dark current histogram corresponding to ``date`` and ``t_ccd``.

    :param date: date in any DateTime format
    :param t_ccd: CCD temperature (deg C)

    :returns: bin_centers, bins, darkhist
    """
    pars = get_sbp_pars(date)
    x = darkbins.bin_centers
    y = smooth_broken_pow(pars, x)

    scale = temp_scalefac(t_ccd)
    xbins = darkbins.bins * scale
    x = x * scale
    return x, xbins, y


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
    try:
        date = DateTime(date).secs
        zeros = np.atleast_1d(np.zeros_like(date + mag + t_ccd))
        dates = date + zeros
        t_ccds = t_ccd + zeros
        mags = mag + zeros
    except:
        raise ValueError("Incompatible input shapes for 'date', 't_ccd', 'mag'")

    warm_fracs = []
    for sdate, stemp  in izip(dates, t_ccds):
        warm_frac = get_warm_fracs(warm_threshold,
                                   date=sdate, T_ccd=stemp)
        warm_fracs.append(warm_frac)
    probs = acq_success_prob(mags, warm_fracs)

    if (np.array(date).ndim == 0 and np.array(t_ccd).ndim == 0
        and np.array(mag).ndim == 0):
           return probs[0]
    else:
        return probs



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




def acq_success_prob(mag, warm_frac, prob_offset=0):
    """
    Calculate probability of acquisition success for a star with ``mag``
    magnitude and a CCD warm fraction ``warm_frac``.  Uses the empirical relation:

       P_acq_success = offset(mag) + scale(mag) * warm_frac

    In ../guide_acq_success/plot_acq_success.py we find the best fit relation:

      log10(scale) = 0.185 + 0.990 * (mag - 10) + -0.491 * (mag - 10)**2
      log10(offset) = -1.489 + 0.888 * (mag - 10) + 0.280 * (mag - 10)**2
    """
    mag10 = mag - 10.0
    scale = 10. ** np.polyval(acq_fit['scale'], mag10)
    offset = 10. ** np.polyval(acq_fit['offset'], mag10)
    # these are actually fits on failure prob, so subtract from 1.00
    return 1.00 - (offset + scale * warm_frac - prob_offset)


def smooth_broken_pow(pars, x):
    """Smoothed broken power-law.  Pars are same as bpl1d (NOT + gaussian sigma):
    1: gamma1
    2: gamma2
    3: x_b (break point)
    4: x_r (normalization reference point)
    5: ampl1
    #   NOT 6: sigma (bins)"""
    (gamma1, gamma2, x_b, x_r, ampl1) = pars
    ampl2 = ampl1 * (x_b / x_r) ** (gamma2 - gamma1)
    ok = xall > x_b
    y = ampl1 * (xall / x_r) ** (-gamma1)
    y[ok] = ampl2 * (xall[ok] / x_r) ** (-gamma2)
    imin = np.searchsorted(xall, x[0] - 1e-3)
    imax = np.searchsorted(xall, x[-1] + 1e-3)
    return np.convolve(y, yg, mode='same')[imin:imax]


def temp_scalefac(T_ccd):
    """Return the multiplicative scale factor to convert a CCD dark map from
    the nominal -19C temperature to the temperature T.  Based on best global fit for
    dark current model in plot_predicted_warmpix.py.  Previous value was 0.62 instead
    of 0.70.
    """
    return exp(log(0.70) / 4.0 * (-19.0 - T_ccd))


def as_array(vals):
    if np.array(vals).ndim == 0:
        is_scalar = True
        vals = np.array([vals])
    else:
        is_scalar = False

    vals = np.array(vals)
    return vals, is_scalar


def get_sbp_pars(dates):
    """
    Return smooth broken powerlaw parameters at ``date``.  This is based on the
    sbp fits for the darkhist_peaknorm histograms, with parameters derived from
    by-hand inspection of fit trending.  See NOTES.
    """
    dates, is_scalar = as_array(dates)
    n_dates = len(dates)

    years = DateTime(dates).frac_year

    ones = np.ones(n_dates)
    g1 = 0.05 * ones
    g2 = 3.15 * ones
    x_r = 50.0 * ones

    ampl = (years - 2000.0) * 1390.2 + 1666.4

    bp_years = np.array([1999.0, 2000.9, 2003.5, 2007.0, 2007.01, 2011.5, 2011.51, 2013.7])
    bp_vals = np.array([125.000, 125.00, 110.00, 117.80, 111.50, 125.00, 108.80, 115.00])
    bp = Ska.Numpy.interpolate(bp_vals, bp_years, years, method='linear')

    if is_scalar:
        g1 = g1[0]
        g2 = g2[0]
        ampl = ampl[0]
        bp = bp[0]
        x_r = x_r[0]

    return g1, g2, bp, x_r, ampl


def get_warm_fracs(warm_threshold, date='2013:001', T_ccd=-19.0):
    x, xbins, y = get_dark_hist(date, T_ccd)
    # First get the full bins to right of warm_threshold
    ii = np.searchsorted(xbins, warm_threshold)
    warmpix = np.sum(y[ii:])
    lx = np.log(warm_threshold)
    lx0 = np.log(xbins[ii - 1])
    lx1 = np.log(xbins[ii])
    ly0 = np.log(y[ii - 1])
    ly1 = np.log(y[ii])
    m = (ly1 - ly0) / (lx1 - lx0)
    partial_bin = y[ii] * (lx1 ** m - lx ** m) / (lx1 ** m - lx0 ** m)
    # print ii, x[ii], xbins[ii - 1], xbins[ii], y[ii], partial_bin
    warmpix += partial_bin

    return warmpix / (1024.0 ** 2)
