import collections
import scipy.stats
import scipy.special
import numpy as np
from astropy.table import Table, vstack

from . import catalogs

import pandas as pd
import numba
from Chandra.Time import DateTime
from agasc import get_star
from cheta import fetch
from Quaternion import Quat
import Ska.quatutil
import mica
from mica.archive import aca_l0
from mica.archive.aca_dark.dark_cal import get_dark_cal_image
from chandra_aca.transform import count_rate_to_mag, pixels_to_yagzag
from cxotime import CxoTime

version = mica.__version__

MAX_MAG = 14
MASK = {
'mouse_bit': np.array([[ True, True, True, True, True, True, True, True],
                       [ True, True, False, False, False, False, True, True],
                       [ True, False, False, False, False, False, False, True],
                       [ True, False, False, False, False, False, False, True],
                       [ True, False, False, False, False, False, False, True],
                       [ True, False, False, False, False, False, False, True],
                       [ True, True, False, False, False, False, True, True],
                       [ True, True, True, True, True, True, True, True]])
}


EXCEPTION_MSG = {
    -1: 'Unknown',
    0: 'OK',
    1: 'No level 0 data',
    2: 'No telemetry data',
    3: 'Mismatch in telemetry between aca_l0 and cheta',
    4: 'Time mismatch between cheta and level0',
    5: 'Failed job'
}
EXCEPTION_CODES = collections.defaultdict(lambda: -1)
EXCEPTION_CODES.update({msg: code for code, msg in EXCEPTION_MSG.items() if code > 0})


class MagStatsException(Exception):
    def __init__(self, msg='', agasc_id=None, obsid=None, timeline_id=None):
        super().__init__(msg)
        self.error_code = EXCEPTION_CODES[msg]
        self.msg = msg
        self.agasc_id = agasc_id
        self.obsid = obsid
        self.timeline_id = timeline_id

    def __str__(self):
        return f'MagStatsException: {self.msg} (agasc_id: {self.agasc_id}, ' \
               f'obsid: {self.obsid}, timeline_id: {self.timeline_id})'

    def __iter__(self):
        yield 'error_code', self.error_code
        yield 'msg', self.msg
        yield 'agasc_id', self.agasc_id
        yield 'obsid', self.obsid
        yield 'timeline_id', self.timeline_id


def _magnitude_correction(time, mag_aca):
    """
    Apply a time-dependent correction to mag_aca.

    :param time: Chandra.Time.DateTime
    :param mag_aca: np.array
    :return: np.array
    """
    params = {"t_ref": "2011-01-01 12:00:00.000",
              "p": [0.005899340720522751,
                    0.12029019332761458,
                    -2.99386247406073e-10,
                    -6.9534637950633265,
                    0.7916261423307238]}

    q = params['p']
    t_ref = DateTime(params['t_ref'])
    dmag = (q[0] + (q[1] + q[2] * np.atleast_1d(time)) *
            np.exp(q[3] + q[4] * np.atleast_1d(mag_aca)))
    dmag[np.atleast_1d(time) < t_ref.secs] = 0
    return np.squeeze(dmag)


def _responsivity(time):
    """
    ACA magnitude response over time.

    This was estimated with bright stars that were observed more than a hundred times during the
    mission. More details in notebook:

    :param time: float
        Time in CXC seconds
    :return:
    """
    a, b, c = [3.24323540e-02, 5.35475665e+08, 8.44499654e+07]
    return - a*(1 + scipy.special.erf((time - b)/c))/2


def _droop_systematic_shift(magnitude):
    """
    Difference between the magnitude determined from DC-subtracted image telemetry and
    the catalog ACA magnitude.

    The magnitude shift is time-independent. It depends only on the catalog magnitude and is zero
    for bright stars. More details in notebook:

    :param magnitude: float
        Catalog ACA magnitude
    :return:
    """
    a, b = [11.24808388,  0.54867984]
    return np.exp((magnitude - a) / b)


def get_star_position(star, slot, telem):
    """
    Residuals for a given AGASC record at a given slot/time.

    :param star:
        Table Row of one AGASC entry
    :param slot: int
    :param times: np.array
        times in CXC seconds, only values exactly matching this time are returned. No interpolation.
    :return:
    """
    aca_misalign = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    R2A = 206264.81

    q = np.array([telem['AOATTQT1'],
                  telem['AOATTQT2'],
                  telem['AOATTQT3'],
                  telem['AOATTQT4']]).transpose()
    norm = np.sum(q**2, axis=1, keepdims=True)
    # I am just normalizing q, just in case.
    n = np.squeeze(np.sqrt(norm))
    q[n != 0] /= np.sqrt(norm)[n != 0]  # prevent warning when dividing by zero (it happens)
    q_att = Quat(q=q)
    Ts = q_att.transform

    star_pos_eci = Ska.quatutil.radec2eci(star['RA_PMCORR'], star['DEC_PMCORR'])
    d_aca = np.dot(np.dot(aca_misalign, Ts.transpose(0, 2, 1)),
                       star_pos_eci).transpose()
    yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
    zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A

    return {
        'yang_star': yag,
        'zang_star': zag,
    }


# this is in case one has to return empty telemetry
_telem_dtype = [('times', 'float64'),
                ('IMGSIZE', 'int32'),
                ('IMGROW0', 'int16'),
                ('IMGCOL0', 'int16'),
                ('IMGRAW', 'float32'),
                ('AOACASEQ', '<U4'),
                ('AOPCADMD', '<U4'),
                ('AOATTQT1', 'float64'),
                ('AOATTQT2', 'float64'),
                ('AOATTQT3', 'float64'),
                ('AOATTQT4', 'float64'),
                ('AOACIIR', '<U3'),
                ('AOACISP', '<U3'),
                ('AOACYAN', 'float64'),
                ('AOACZAN', 'float64'),
                ('AOACMAG', 'float32'),
                ('AOACFCT', '<U4'),
                ('mags_img', 'float64'),
                ('yang_img', 'float64'),
                ('zang_img', 'float64'),
                ('yang_star', 'float64'),
                ('zang_star', 'float64'),
                ('mags', 'float64'),
                ('dy', 'float64'),
                ('dz', 'float64'),
                ('dr', 'float64')]


def get_telemetry(obs):
    """
    Get all telemetry relevant for the mag_stats task.

    This gets:
    - AOACASEQ
    - AOPCADMD
    - AOACMAG (ACA estimated magnitude)
    - AOACIIR (ACA ionizing radiation flag)
    - AOACISP (ACA saturated pixel flag)

    MSIDs are renamed to remove the slot number.
    This assumes all MSIDs occur at the same times (they do)

    :param obs: astropy.table.Row
        It must have the following columns: 'agasc_id', 'mp_starcat_time', 'mag', 'slot'
    :return: dict
    """
    catalogs.load()
    dwell = catalogs.DWELLS_NP[catalogs.DWELLS_MAP[obs['mp_starcat_time']]]
    star = get_star(obs['agasc_id'], date=dwell['tstart'])
    start = dwell['tstart']
    stop = dwell['tstop']
    slot = obs['slot']

    # first we get slot data from mica and magnitudes from cheta and match them in time
    # to match them in time, we assume they come in steps of 1.025 seconds, starting from the first
    # time sample.
    slot_data_cols = ['TIME', 'END_INTEG_TIME', 'IMGSIZE', 'IMGROW0', 'IMGCOL0', 'TEMPCCD', 'IMGRAW']
    slot_data = aca_l0.get_slot_data(start, stop, slot=obs['slot'],
                                     img_shape_8x8=True, columns=slot_data_cols)

    msid = fetch.MSID(f'AOACMAG{slot}', start, stop)
    if len(slot_data) == 0:
        raise MagStatsException('No level 0 data',
                                agasc_id=obs["agasc_id"],
                                obsid=obs["obsid"],
                                timeline_id=obs["timeline_id"])
    tmin = np.min([np.min(slot_data['END_INTEG_TIME']), np.min(msid.times)])
    t1 = np.round((msid.times - tmin)/1.025)
    t2 = np.round((slot_data['END_INTEG_TIME'].data - tmin)/1.025)
    c, i1, i2 = np.intersect1d(t1, t2, return_indices=True)
    times = msid.times[i1]

    # the following line removes a couple of points at the edges. I have not checked why they differ
    slot_data = slot_data[i2]

    if len(times) == 0:
        # the intersection was null.
        raise MagStatsException('Time mismatch between cheta and level0',
                                agasc_id=obs["agasc_id"],
                                obsid=obs["obsid"],
                                timeline_id=obs["timeline_id"])

    # Now that we have the times, we get the rest of the MSIDs
    telem = {
        'times': times
    }
    telem.update({k: slot_data[k] for k in slot_data_cols[2:]})

    names = ['AOACASEQ', 'AOPCADMD',
             f'AOACIIR{slot}', f'AOACISP{slot}', f'AOACMAG{slot}', f'AOACFCT{slot}',
             f'AOACZAN{slot}', f'AOACYAN{slot}'] + [f'AOATTQT{i}' for i in range(1, 5)]
    msids = fetch.MSIDset(names, times[0] - 4, times[-1] + 4)
    # the following just works...
    t = np.in1d(msids[names[0]].times, times)
    telem.update({n: msids[n].vals[t] for n in names})
    if len(telem['AOACASEQ']) != len(telem['IMGSIZE']):
        raise MagStatsException(
            "Mismatch in telemetry between aca_l0 and cheta",
            agasc_id=obs['agasc_id'], obsid=obs['obsid'], timeline_id=obs['timeline_id']
        )
    for name in ['AOACIIR', 'AOACISP', 'AOACYAN', 'AOACZAN', 'AOACMAG', 'AOACFCT']:
        telem[name] = telem[f'{name}{slot}']
        del telem[f'{name}{slot}']
    for name in ['AOACIIR', 'AOACISP']:
        telem[name] = np.char.rstrip(telem[name])
    ok = (telem['AOACASEQ'] == 'KALM') & (telem[f'AOACIIR'] == 'OK') & \
         (telem[f'AOACISP'] == 'OK') & (telem['AOPCADMD'] == 'NPNT') & \
         (telem['AOACFCT'] == 'TRAK')

    # etc...
    telem.update(get_mag_from_img(slot_data, start, ok))
    telem.update(get_star_position(star=star, slot=obs['slot'], telem=telem))

    droop_shift = _droop_systematic_shift(star['MAG_ACA'])
    responsivity = _responsivity(start)
    telem['mags'] = telem['mags_img'] - responsivity - droop_shift
    telem['mags'][~ok] = 0.
    telem['ok'] = ok

    telem['dy'] = np.ones(len(ok)) * np.inf
    telem['dz'] = np.ones(len(ok)) * np.inf
    telem['dr'] = np.ones(len(ok)) * np.inf
    if np.any(ok):
        yang = telem['yang_img'] - telem['yang_star']
        zang = telem['zang_img'] - telem['zang_star']
        y25, y50, y75 = np.quantile(yang[ok], [0.25, 0.5, 0.75])
        z25, z50, z75 = np.quantile(zang[ok], [0.25, 0.5, 0.75])
        centroid_outlier = ((yang > y75 + 3 * (y75 - y25)) |
                            (yang < y25 - 3 * (y75 - y25)) |
                            (zang > z75 + 3 * (z75 - z25)) |
                            (zang < z25 - 3 * (z75 - z25)))

        telem['dy'] = yang - np.mean(yang[ok & ~centroid_outlier])
        telem['dz'] = zang - np.mean(zang[ok & ~centroid_outlier])
        telem['dr'] = (telem['dy'] ** 2 + telem['dz'] ** 2) ** .5

    return telem


def astropy_table_cache(name, dir):
    def decorator_cache(func):
        import inspect
        signature = inspect.signature(func)
        from functools import wraps
        @wraps(func)
        def wrapper(*args, name=name, **kwargs):
            import os
            name += '::'
            s_args = signature.bind(*args, **kwargs).arguments
            arg_str = '-'.join(
                ['{a}:{v}'.format(a=a, v=s_args[a]) for a in s_args])
            filename = f'{name}{arg_str}.fits'
            filename = os.path.join(dir, filename)
            if os.path.exists(filename):
                # Table.read returns a big-endian array, and as_array converts to native endianness
                return Table(Table.read(filename).as_array())
            else:
                result = func(*args, **kwargs)
                result.write(filename)
                return result
        return wrapper
    return decorator_cache


#@astropy_table_cache(name='telem', dir='/Users/javierg/SAO/mica/mag_stats_cache/telem')
def get_telemetry_by_agasc_id(agasc_id, obsid=None, ignore_exceptions=False):
    """
    Get all telemetry relevant for the mag_stats task for a given AGASC ID.

    :param agasc_id: int
    :param obsid: int (optional)
    :return: dict
    """
    catalogs.load()
    if obsid is None:
        obs = catalogs.STARS_OBS[
            (catalogs.STARS_OBS['agasc_id'] == agasc_id)]
    else:
        obs = catalogs.STARS_OBS[(catalogs.STARS_OBS['agasc_id'] == agasc_id) &
                                 (catalogs.STARS_OBS['obsid'] == obsid)]
    if len(obs) > 1:
        obs = obs.loc['mp_starcat_time', sorted(obs['mp_starcat_time'])]
    telem = []
    for i, o in enumerate(obs):
        try:
            t = Table(get_telemetry(o))
            t['obsid'] = o['obsid']
            t['agasc_id'] = agasc_id
            telem.append(t)
        except Exception as e:
            import sys
            import traceback
            print(f'{agasc_id} failed', e)
            exc_type, exc_value, exc_traceback = sys.exc_info()
            trace = traceback.extract_tb(exc_traceback)
            print(f'{exc_type.__name__} {exc_value}')
            for step in trace:
                print(f'  in {step.filename}:{step.lineno}/{step.name}:')
                print(f'    {step.line}')
            if not ignore_exceptions:
                raise
    return vstack(telem)


def add_obsid_info(telem, obs_stats):
    obs_stats['obsid_ok'] = (
        (obs_stats['n'] > 10) &
        (obs_stats['f_track'] > 0.3) &
        (obs_stats['dr_star'] < 10) &
        (obs_stats['lf_variability_100s'] < 1)
    )
    obs_stats['comments'] = np.zeros(len(obs_stats), dtype='<U80')

    telem = vstack(telem)
    telem['obsid_ok'] = True
    telem['obsid_outlier'] = False

    for s in obs_stats:
        obsid = s['obsid']
        o = (telem['obsid'] == obsid)
        telem['obsid_ok'][o] = np.ones(np.sum(o), dtype=bool) * s['obsid_ok']
        if (np.any(telem['ok'][o]) and s['f_track'] > 0 and
                np.isfinite(s['q75']) and np.isfinite(s['q25'])):
            iqr = s['q75'] - s['q25']
            telem['obsid_outlier'][o] = (telem[o]['ok'] & (iqr > 0) &
                                         (telem[o]['mags'] < s['q25'] - 1.5*iqr) |
                                         (telem[o]['mags'] > s['q75'] + 1.5*iqr))
    return telem


@numba.jit(nopython=True)
def staggered_aca_slice(array_in, array_out, row, col):
    for i in np.arange(len(row)):
        if row[i]+8 < 1024 and col[i]+8 < 1024:
            array_out[i] = array_in[row[i]:row[i]+8, col[i]:col[i]+8]


def get_mag_from_img(slot_data, t_start, ok=True):
    dark_cal = get_dark_cal_image(t_start, 'nearest',
                                  t_ccd_ref=np.mean(slot_data['TEMPCCD'] - 273.16),
                                  aca_image=False)

    # all images will be 8x8, with a centered mask, imgrow will always be the one of the 8x8 corner.
    imgrow_8x8 = np.where(slot_data['IMGSIZE'] == 8,
                          slot_data['IMGROW0'],
                          slot_data['IMGROW0'] - 1
                          )
    imgcol_8x8 = np.where(slot_data['IMGSIZE'] == 8,
                          slot_data['IMGCOL0'],
                          slot_data['IMGCOL0'] - 1
                          )

    # subtract closest dark cal
    dark = np.zeros([len(slot_data), 8, 8], dtype=np.float64)
    staggered_aca_slice(dark_cal.astype(float), dark, 512 + imgrow_8x8, 512 + imgcol_8x8)
    img_sub = slot_data['IMGRAW'] - dark * 1.696 / 5
    img_sub.mask *= MASK['mouse_bit']

    # calculate magnitude
    mag = np.ones(len(slot_data)) * MAX_MAG
    counts = np.ma.sum(np.ma.sum(img_sub, axis=1), axis=1)
    m = ok & np.isfinite(counts) & (counts > 0)
    mag[m] = count_rate_to_mag(counts[m] * 5 / 1.7)
    # this extra step is to investigate the background scale
    dark = np.ma.array(dark * 1.696 / 5, mask=img_sub.mask)
    img_raw = np.ma.array(slot_data['IMGRAW'], mask=img_sub.mask)
    dark_count = np.ma.sum(np.ma.sum(dark, axis=1), axis=1)
    img_count = np.ma.sum(np.ma.sum(img_raw, axis=1), axis=1)

    # centroids
    yag = np.zeros(len(slot_data))
    zag = np.zeros(len(slot_data))
    pixel_center = np.arange(8) + 0.5
    projected_image = np.ma.sum(slot_data['IMGRAW'], axis=1)
    col = np.ma.sum(pixel_center * projected_image, axis=1) / np.ma.sum(projected_image, axis=1)
    projected_image = np.ma.sum(slot_data['IMGRAW'], axis=2)
    row = np.ma.sum(pixel_center * projected_image, axis=1) / np.ma.sum(projected_image, axis=1)

    y_pixel = row + imgrow_8x8
    z_pixel = col + imgcol_8x8
    yag[m], zag[m] = pixels_to_yagzag(y_pixel[m], z_pixel[m])
    return {
        'mags_img': mag,
        'yang_img': yag,
        'zang_img': zag,
        'counts_img': img_count,
        'counts_dark': dark_count
    }


def get_obsid_stats(obs, telem=None):
    """
    Get summary magnitude statistics for an observation.

    :param obs: astropy.table.Row
        a "star observation" row. From the join of starcheck catalog and starcat commands
        It must have the following columns: 'agasc_id', 'mp_starcat_time', 'mag', 'slot'
    :param telem: dict
        Dictionary with telemetry (output of get_telemetry)
    :return: dict
        dictionary with stats
    """
    catalogs.load()
    if telem is None:
        telem = get_telemetry(obs)

    star = get_star(obs['agasc_id'])
    dwell = catalogs.DWELLS_NP[catalogs.DWELLS_MAP[obs['mp_starcat_time']]]
    start = dwell['tstart']
    stop = dwell['tstop']

    stats = {k: obs[k] for k in
             ['agasc_id', 'obsid', 'slot', 'type', 'mp_starcat_time', 'timeline_id']}
    droop_shift = _droop_systematic_shift(star['MAG_ACA'])
    responsivity = _responsivity(start)
    stats.update({'tstart': start,
                  'tstop': stop,
                  'mag_correction': - responsivity - droop_shift,
                  'responsivity': responsivity,
                  'droop_shift': droop_shift,
                  'mag_aca': star['MAG_ACA'],
                  'mag_aca_err': star['MAG_ACA_ERR'] / 100,
                  'row': obs['row'],
                  'col': obs['col'],
                  })
    stats.update(calc_obsid_stats(telem))

    return stats


def calc_obsid_stats(telem):
    """
    Get summary magnitude statistics for an observation.

    :param telem: dict
        Dictionary with telemetry (output of get_telemetry)
    :return: dict
        dictionary with stats
    """
    times = telem['times']

    kalman = (telem['AOACASEQ'] == 'KALM') & (telem['AOPCADMD'] == 'NPNT')
    track = (telem[f'AOACIIR'] == 'OK') & (telem[f'AOACISP'] == 'OK') & \
            (telem['AOACFCT'] == 'TRAK')
    dr3 = (telem['dr'] < 3)
    dr5 = (telem['dr'] < 5)

    f_kalman = np.sum(kalman) / len(kalman)
    f_track = np.sum(kalman & track) / np.sum(kalman)
    f_3 = np.sum(kalman & track & dr3) / np.sum(kalman & track)
    f_5 = np.sum(kalman & track & dr5) / np.sum(kalman & track)

    ok = kalman & track & dr5
    f_ok = np.sum(ok) / len(ok)

    if np.any(ok):
        yang_mean = np.mean(telem['yang_img'][ok] - telem['yang_star'][ok])
        zang_mean = np.mean(telem['zang_img'][ok] - telem['zang_star'][ok])
        dr_star = np.sqrt(yang_mean**2 + zang_mean**2)
    else:
        dr_star = np.inf

    stats = {
        'aoacmag_mean': np.inf,
        'aoacmag_err': np.inf,
        'aoacmag_q25': np.inf,
        'aoacmag_median': np.inf,
        'aoacmag_q75': np.inf,
        'counts_img': np.inf,
        'counts_dark': np.inf,
        'f_kalman': f_kalman,
        'f_track': f_track,
        'f_dr5': f_5,
        'f_dr3': f_3,
        'f_ok': f_ok,
        'q25': np.inf,
        'median': np.inf,
        'q75': np.inf,
        'mean': np.inf,
        'mean_err': np.inf,
        'std': np.inf,
        'skew': np.inf,
        'kurt': np.inf,
        't_mean': np.inf,
        't_mean_err': np.inf,
        't_std': np.inf,
        't_skew': np.inf,
        't_kurt': np.inf,
        'n': len(telem['AOACMAG']),
        'n_ok': np.sum(ok),
        'outliers': -1,
        'lf_variability_100s': np.inf,
        'lf_variability_500s': np.inf,
        'lf_variability_1000s': np.inf,
        'tempccd': np.nan,
        'dr_star': dr_star,
    }
    if stats['n_ok'] < 10:
        return stats

    aoacmag_q25, aoacmag_q50, aoacmag_q75 = np.quantile(telem['AOACMAG'][ok], [0.25, 0.5, 0.75])

    mags = telem['mags']
    q25, q50, q75 = np.quantile(mags[ok], [0.25, 0.5, 0.75])
    iqr = q75 - q25
    outlier = ok & ((mags > q75 + 3 * iqr) | (mags < q25 - 3 * iqr))

    dt = np.mean(np.diff(times))
    s = pd.Series(mags[ok], index=times[ok])
    s_100s = s.rolling(window=int(100 / dt), min_periods=1, center=True).median().dropna()
    s_500s = s.rolling(window=int(500 / dt), center=True).median().dropna()
    s_1000s = s.rolling(window=int(1000 / dt), center=True).median().dropna()

    stats.update({
        'aoacmag_mean': np.mean(telem['AOACMAG'][ok]),
        'aoacmag_err': np.std(telem['AOACMAG'][ok]),
        'aoacmag_q25': aoacmag_q25,
        'aoacmag_median': aoacmag_q50,
        'aoacmag_q75': aoacmag_q75,
        'q25': q25,
        'median': q50,
        'q75': q75,
        'counts_img': np.mean(telem['counts_img'][ok]),
        'counts_dark': np.mean(telem['counts_dark'][ok]),
        'mean': np.mean(mags[ok]),
        'mean_err': scipy.stats.sem(mags[ok]),
        'std': np.std(mags[ok]),
        'skew': scipy.stats.skew(mags),
        'kurt': scipy.stats.kurtosis(mags),
        't_mean': np.mean(mags[ok & (~outlier)]),
        't_mean_err': scipy.stats.sem(mags[ok & (~outlier)]),
        't_std': np.std(mags[ok & (~outlier)]),
        't_skew': scipy.stats.skew(mags[ok & (~outlier)]),
        't_kurt': scipy.stats.kurtosis(mags[ok & (~outlier)]),
        'outliers': np.sum(outlier),
        'lf_variability_100s': np.max(s_100s) - np.min(s_100s),
        'lf_variability_500s': np.max(s_500s) - np.min(s_500s),
        'lf_variability_1000s': np.max(s_1000s) - np.min(s_1000s),
        'tempccd': np.mean(telem['TEMPCCD'][ok]) - 273.16,
    })
    return stats


def get_agasc_id_stats(agasc_id, excluded_observations={}, tstop=None):
    """
    Get summary magnitude statistics for an AGASC ID.

    :param agasc_id:
    :return: dict
        dictionary with stats
    """
    min_mag_obs_err = 0.3

    catalogs.load(tstop=tstop)
    # Get a table of every time the star has been observed
    idx0, idx1 = catalogs.STARS_OBS_MAP[agasc_id]
    star_obs = catalogs.STARS_OBS[idx0:idx1]
    if tstop is not None:
        star_obs[star_obs['mp_starcat_time'] <= tstop]
    if len(star_obs) > 1:
        star_obs = star_obs.loc['mp_starcat_time', sorted(star_obs['mp_starcat_time'])]

    failures = []
    all_telem = []
    stats = []
    for i, o in enumerate(star_obs):
        try:
            telem = Table(get_telemetry(o))
            all_telem.append(telem)
            stats.append(get_obsid_stats(o, telem={k: telem[k] for k in telem.colnames}))
            last_obs_time = CxoTime(o['mp_starcat_time']).cxcsec
        except MagStatsException as e:
            failures.append(dict(e))

    if len(all_telem) == 0:
        raise MagStatsException('No telemetry data', agasc_id=agasc_id)

    stats = Table(stats)
    n_obsids = len(star_obs)

    stats['obsid_ok'] = (
        (stats['n'] > 10) &
        (stats['f_track'] > 0.3) &
        (stats['dr_star'] < 10) &
        (stats['lf_variability_100s'] < 1)
    )
    stats['comments'] = np.zeros(len(stats), dtype='<U80')
    if excluded_observations:
        excluded_mask = np.in1d(stats['obsid'], np.fromiter(excluded_observations, dtype=int))
        stats['obsid_ok'] *= ~excluded_mask
        for i in np.argwhere(excluded_mask).T[0]:
            stats['comments'][i] = excluded_observations[stats[i]['obsid']]

    for s, t in zip(stats, all_telem):
        t['obsid_ok'] = np.ones_like(t['ok'], dtype=bool) * s['obsid_ok']
        t['obsid_outlier'] = np.zeros_like(t['ok'])
        if np.any(t['ok']) and s['f_track'] > 0 and s['obsid_ok']:
            iqr = s['q75'] - s['q25']
            t['obsid_outlier'] = (t['ok'] & (iqr > 0) &
                                  (t['mags'] < s['q25'] - 1.5*iqr) | (t['mags'] > s['q75'] + 1.5*iqr))
    all_telem = vstack([Table(t) for t in all_telem])

    mags = all_telem['mags']
    ok = all_telem['ok'] & all_telem['obsid_ok']

    f_ok = np.sum(ok)/len(ok)

    star = get_star(agasc_id, date=all_telem['times'][0])
    result = {
        'last_obs_time': last_obs_time,
        'agasc_id': agasc_id,
        'mag_aca': star['MAG_ACA'],
        'mag_aca_err': star['MAG_ACA_ERR']/100,
        'mag_obs': 0.,
        'mag_obs_err': min_mag_obs_err,
        'mag_obs_std': 0.,
        'color': star['COLOR1'],
        'n_obsids': n_obsids,
        'n_obsids_fail': len(failures),
        'n_obsids_ok': np.sum(stats['obsid_ok']),
        'n_no_track': np.sum(stats['f_ok'] < 0.3),
        'n': len(ok),
        'n_ok': np.sum(ok),
        'f_ok': f_ok,
        'median': 0,
        'sigma_minus': 0,
        'sigma_plus': 0,
        'mean': 0,
        'std': 0,
        'mag_weighted_mean': 0,
        'mag_weighted_std': 0,
        't_mean': 0,
        't_std': 0,
        'n_outlier': 0,
        't_mean_1': 0,
        't_std_1': 0,
        'n_outlier_1': 0,
        't_mean_2': 0,
        't_std_2': 0,
        'n_outlier_2': 0,
        # these are the criteria for including in supplement
        'selected_atol': False,
        'selected_rtol': False,
        'selected_mag_aca_err': False,
        'selected_color': False
    }

    for dr in [3, 5]:
        result.update({
            f't_mean_dr{dr}': 0,
            f't_std_dr{dr}': 0,
            f't_mean_dr{dr}_not': 0,
            f't_std_dr{dr}_not': 0,
            f'mean_dr{dr}': 0,
            f'std_dr{dr}': 0,
            f'f_dr{dr}': 0,
            f'n_dr{dr}': 0,
            f'n_dr{dr}_outliers': 0,
            f'median_dr{dr}': 0,
            f'sigma_minus_dr{dr}': 0,
            f'sigma_plus_dr{dr}': 0,
        })

    if result['n_ok'] < 10:
        return result, stats, failures

    sigma_minus, q25, median, q75, sigma_plus = np.quantile(mags[ok],
                                                            [0.158, 0.25, 0.5, 0.75, 0.842])
    iqr = q75 - q25
    outlier_1 = ok & ((mags > q75 + 1.5 * iqr) | (mags < q25 - 1.5 * iqr))
    outlier_2 = ok & ((mags > q75 + 3 * iqr) | (mags < q25 - 3 * iqr))
    outlier = all_telem['obsid_outlier']

    # combine measurements using a weighted mean
    obs_ok = stats['obsid_ok']
    min_std = max(0.1, stats[obs_ok]['std'].min())
    stats['w'] = np.nan
    stats['mean_corrected'] = np.nan
    stats['weighted_mean'] = np.nan
    stats['w'][obs_ok] = np.where(stats['std'][obs_ok] != 0,
                                  1. / stats['std'][obs_ok],
                                  1. / min_std)
    stats['mean_corrected'][obs_ok] = stats['t_mean'][obs_ok] + stats['mag_correction'][obs_ok]
    stats['weighted_mean'][obs_ok] = stats['mean_corrected'][obs_ok] * stats['w'][obs_ok]

    mag_weighted_mean = (stats[obs_ok]['weighted_mean'].sum() / stats[obs_ok]['w'].sum())
    mag_weighted_std = (
        np.sqrt(((stats[obs_ok]['mean'] - mag_weighted_mean)**2 * stats[obs_ok]['w']).sum() /
                stats[obs_ok]['w'].sum())
    )

    result.update({
        'agasc_id': agasc_id,
        'n': len(ok),
        'n_ok': np.sum(ok),
        'f_ok': f_ok,
        'median': median,
        'sigma_minus': sigma_minus,
        'sigma_plus': sigma_plus,
        'mean': np.mean(mags[ok]),
        'std': np.std(mags[ok]),
        'mag_weighted_mean': mag_weighted_mean,
        'mag_weighted_std': mag_weighted_std,
        't_mean': np.mean(mags[ok & (~outlier)]),
        't_std': np.std(mags[ok & (~outlier)]),
        'n_outlier': np.sum(ok & outlier),
        't_mean_1': np.mean(mags[ok & (~outlier_1)]),
        't_std_1': np.std(mags[ok & (~outlier_1)]),
        'n_outlier_1': np.sum(ok & outlier_1),
        't_mean_2': np.mean(mags[ok & (~outlier_2)]),
        't_std_2': np.std(mags[ok & (~outlier_2)]),
        'n_outlier_2': np.sum(ok & outlier_2),
    })

    for dr in [3, 5]:
        k = ok & (all_telem['dr'] < dr)
        k2 = ok & (all_telem['dr'] >= dr)
        if not np.any(k):
            continue
        sigma_minus, q25, median, q75, sigma_plus = np.quantile(mags[k],
                                                                [0.158, 0.25, 0.5, 0.75, 0.842])
        outlier = ok & all_telem['obsid_outlier']
        mag_not = np.nanmean(mags[k2 & (~outlier)]) if np.sum(k2 & (~outlier)) else np.nan
        std_not = np.nanstd(mags[k2 & (~outlier)]) if np.sum(k2 & (~outlier)) else np.nan
        result.update({
            f't_mean_dr{dr}': np.mean(mags[k & (~outlier)]),
            f't_std_dr{dr}': np.std(mags[k & (~outlier)]),
            f't_mean_dr{dr}_not': mag_not,
            f't_std_dr{dr}_not': std_not,
            f'mean_dr{dr}': np.mean(mags[k]),
            f'std_dr{dr}': np.std(mags[k]),
            f'f_dr{dr}': np.sum(k) / np.sum(ok),
            f'n_dr{dr}': np.sum(k),
            f'n_dr{dr}_outliers': np.sum(k & outlier),
            f'median_dr{dr}': median,
            f'sigma_minus_dr{dr}': sigma_minus,
            f'sigma_plus_dr{dr}': sigma_plus,
        })

    result.update({
        'mag_obs': result['t_mean_dr5'],
        'mag_obs_err': np.sqrt(result['t_std_dr5']**2 + min_mag_obs_err**2),
        'mag_obs_std': result['t_std_dr5'],
    })

    # these are the criteria for including in supplement
    result.update({
        'selected_atol': np.abs(result['mag_obs'] - result['mag_aca']) > 0.3,
        'selected_rtol': np.abs(result['mag_obs'] - result['mag_aca']) > 3 * result['mag_aca_err'],
        'selected_mag_aca_err': result['mag_aca_err'] > 0.2,
        'selected_color': (result['color'] == 1.5) | (np.isclose(result['color'], 0.7))
    })

    return result, stats, failures
