import scipy.stats
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

version = mica.__version__

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

    q_att = Quat(q=np.array([telem['AOATTQT1'],
                             telem['AOATTQT2'],
                             telem['AOATTQT3'],
                             telem['AOATTQT4']]).transpose())
    Ts = q_att.transform

    star_pos_eci = Ska.quatutil.radec2eci(star['RA_PMCORR'], star['DEC_PMCORR'])
    d_aca = np.dot(np.dot(aca_misalign, Ts.transpose(0, 2, 1)),
                   star_pos_eci).transpose()
    yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
    zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A

    return {
        'star_yag': yag,
        'star_zag': zag,
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
                ('mags_img', 'float64'),
                ('yag_img', 'float64'),
                ('zag_img', 'float64'),
                ('star_yag', 'float64'),
                ('star_zag', 'float64'),
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
    dwell = catalogs.DWELLS_NP[catalogs.DWELLS_MAP[obs['mp_starcat_time']]]
    star = get_star(obs['agasc_id'], date=dwell['tstart'])
    start = dwell['tstart']
    stop = dwell['tstop']
    slot = obs['slot']

    # first we get slot data from mica and magnitudes from cheta and match them in time
    slot_data_cols = ['TIME', 'END_INTEG_TIME', 'IMGSIZE', 'IMGROW0', 'IMGCOL0', 'IMGRAW']
    slot_data = aca_l0.get_slot_data(start, stop, slot=obs['slot'],
                                     img_shape_8x8=True, columns=slot_data_cols)

    msid = fetch.MSID(f'AOACMAG{slot}', start, stop)
    t1 = np.round(msid.times, 3)
    t2 = np.round(slot_data['END_INTEG_TIME'], 3)
    c, i1, i2 = np.intersect1d(t1, t2, return_indices=True)
    times = msid.times[i1]

    # the following line removes a couple of points at the edges. I have not checked why they differ
    slot_data = slot_data[i2]

    if len(times) == 0:
        # the intersection was null.
        return {k: np.array([], dtype=v) for k, v in _telem_dtype}
    # Now that we have the times, we get the rest of the MSIDs
    telem = {
        'times': times
    }
    telem.update({k: slot_data[k] for k in slot_data_cols[2:]})

    names = ['AOACASEQ', 'AOPCADMD', f'AOACIIR{slot}', f'AOACISP{slot}', f'AOACMAG{slot}',
             f'AOACZAN{slot}', f'AOACYAN{slot}'] + [f'AOATTQT{i}' for i in range(1, 5)]
    msids = fetch.MSIDset(names, times[0] - 4, times[-1] + 4)
    # the following just works...
    t = np.in1d(msids[names[0]].times, times)
    telem.update({n: msids[n].vals[t] for n in names})
    for name in ['AOACIIR', 'AOACISP', 'AOACYAN', 'AOACZAN', 'AOACMAG']:
        telem[name] = telem[f'{name}{slot}']
        del telem[f'{name}{slot}']

    # etc...
    telem.update(get_mag_from_img(slot_data, start))
    telem.update(get_star_position(star=star, slot=obs['slot'], telem=telem))

    telem['mags'] = telem['mags_img']
    #dmag = _magnitude_correction(start, obs['mag'])
    #telem['mags'] = telem['AOACMAG']
    #telem['mags_corrected'] = telem['AOACMAG'] - dmag

    telem['dy'] = telem['yag_img'] - telem['star_yag']
    telem['dz'] = telem['zag_img'] - telem['star_zag']
    #telem['dy'] = telem[f'AOACYAN{slot}'] - telem['star_yag']
    #telem['dz'] = telem[f'AOACZAN{slot}'] - telem['star_zag']
    # cheating here and ignoring spherical trig
    telem['dr'] = (telem['dy'] ** 2 + telem['dz'] ** 2) ** .5

    return telem


def get_telemetry_by_agasc_id(agasc_id, obsid=None):
    """
    Get all telemetry relevant for the mag_stats task for a given AGASC ID.

    :param agasc_id: int
    :param obsid: int (optional)
    :return: dict
    """
    if obsid is None:
        obs = catalogs.STARS_OBS[
            (catalogs.STARS_OBS['agasc_id'] == agasc_id)]
    else:
        obs = catalogs.STARS_OBS[(catalogs.STARS_OBS['agasc_id'] == agasc_id) &
                                 (catalogs.STARS_OBS['obsid'] == obsid)]
    telem = [Table(get_telemetry(o)) for o in obs]
    for i, obsid in enumerate(obs['obsid']):
        telem[i]['obsid'] = obsid
    return vstack(telem)


@numba.jit(nopython=True)
def staggered_aca_slice(array_in, array_out, row, col):
    for i in np.arange(len(row)):
        if row[i]+8 < 1024 and col[i]+8 < 1024:
            array_out[i] = array_in[row[i]:row[i]+8, col[i]:col[i]+8]


def get_mag_from_img(slot_data, t_start):
    dark_cal = get_dark_cal_image(t_start, 'nearest', t_ccd_ref=-11.2, aca_image=False)

    # all images will be 8x8, with a centered mask, imgrow will always be the one of the 8x8 corner.
    imgrow_8x8 = np.where(slot_data['IMGSIZE'] == 8,
                          slot_data['IMGROW0'],
                          slot_data['IMGROW0'] - 1
                          )
    imgcol_8x8 = np.where(slot_data['IMGSIZE'] == 8,
                          slot_data['IMGCOL0'],
                          slot_data['IMGCOL0'] - 1
                          )
    # a useful mask to select entries when it is tracking
    m = slot_data['IMGSIZE'] > 4

    # subtract closest dark cal
    dark = np.zeros([len(slot_data), 8, 8], dtype=np.float64)
    staggered_aca_slice(dark_cal.astype(float), dark, 512 + imgrow_8x8, 512 + imgcol_8x8)
    img_sub = slot_data['IMGRAW'] - dark * 1.696 / 5
    img_sub.mask *= MASK['mouse_bit']

    # calculate magnitude
    mag = np.zeros(len(slot_data))
    counts = np.ma.sum(np.ma.sum(img_sub, axis=1), axis=1)
    mag[m] = count_rate_to_mag(counts[m] * 5 / 1.7)

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
        'yag_img': yag,
        'zag_img': zag
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
    if telem is None:
        telem = get_telemetry(obs)

    dwell = catalogs.DWELLS_NP[catalogs.DWELLS_MAP[obs['mp_starcat_time']]]
    start = dwell['tstart']
    stop = dwell['tstop']

    stats = {k: obs[k] for k in obs.colnames}
    stats.update({'tstart': start,
                  'tstop': stop})
    stats.update(calc_obsid_stats(telem))

    dmag = _magnitude_correction(start, obs['mag'])
    stats['mag_correction'] = dmag[()]
    stats['mean_corrected'] = stats['mean'] - dmag

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
    mags = telem['mags']

    track = (telem['AOACASEQ'] == 'KALM') & (telem[f'AOACIIR'] == 'OK ') & \
            (telem[f'AOACISP'] == 'OK ') & (telem['AOPCADMD'] == 'NPNT') & \
            (telem['IMGSIZE'] > 4)
    dr3 = (telem['dr'] < 3)
    dr5 = (telem['dr'] < 5)

    f_track = sum(track) / len(track)
    f_3 = sum(track & dr3) / len(track)
    f_5 = sum(track & dr5) / len(track)
    if sum(track & dr3):
        f_14 = sum(track & dr3 & (mags >= 13.9)) / sum(track & dr3)
    else:
        f_14 = np.nan

    ok = track & dr3 & (mags < 13.9)
    f_ok = sum(ok) / len(ok)

    stats = {
        'f_track': f_track,
        'f_dr5': f_5,
        'f_dr3': f_3,
        'f_14': f_14,
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
        'n': len(mags),
        'n_ok': sum(ok),
        'outliers': -1,
        'lf_variability_100s': np.inf,
        'lf_variability_500s': np.inf,
        'lf_variability_1000s': np.inf,
    }
    if sum(ok) < 10:
        return stats

    q25, q50, q75 = np.quantile(mags[ok], [0.25, 0.5, 0.75])
    iqr = q75 - q25
    outlier = ok & ((mags > q75 + 3 * iqr) | (mags < q25 - 3 * iqr))

    dt = np.mean(np.diff(times))
    s = pd.Series(mags[ok], index=times[ok])
    s_100s = s.rolling(window=int(100 / dt), min_periods=1, center=True).median().dropna()
    s_500s = s.rolling(window=int(500 / dt), center=True).median().dropna()
    s_1000s = s.rolling(window=int(1000 / dt), center=True).median().dropna()

    stats.update({
        'q25': q25,
        'median': q50,
        'q75': q75,
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
        'outliers': sum(outlier),
        'lf_variability_100s': np.max(s_100s) - np.min(s_100s),
        'lf_variability_500s': np.max(s_500s) - np.min(s_500s),
        'lf_variability_1000s': np.max(s_1000s) - np.min(s_1000s),
    })
    return stats


def get_agasc_id_stats(agasc_id):
    """
    Get summary magnitude statistics for an AGASC ID.

    :param agasc_id:
    :return: dict
        dictionary with stats
    """
    # Get a table of every time the star has been observed
    idx0, idx1 = catalogs.STARS_OBS_MAP[agasc_id]
    star_obs = catalogs.STARS_OBS[idx0:idx1]

    all_telem = []
    for obs in star_obs:
        telem = get_telemetry(obs)
        telem['obsid'] = np.ones_like(telem['times']) * obs['obsid']
        telem['agasc_id'] = np.ones_like(telem['times']) * agasc_id
        if len(telem['times']):
            all_telem.append(telem)

    if len(all_telem) == 0:
        raise Exception(f'No telemetry data for agasc_id {agasc_id}')

    stats = Table([get_obsid_stats(obs, telem=telem) for obs, telem in zip(star_obs, all_telem)])
    n_obsids = len(stats)
    stats = stats[stats['f_dr3'] > 0]

    stats['obsid_ok'] = (
        (stats['n'] > 10) &
        (stats['f_ok'] > 0.3) &
        (stats['lf_variability_100s'] < 1) &
        ((stats['q75'] - stats['q25'] == 0) | (stats['outliers'] < 0.05 * stats['n']))
    )

    all_telem = vstack([Table(t) for t in all_telem])

    mags = all_telem['mags']
    ok = (all_telem['AOACASEQ'] == 'KALM') & (all_telem['AOACIIR'] == 'OK ') & \
         (all_telem['AOACISP'] == 'OK ') & (all_telem['AOPCADMD'] == 'NPNT') & \
        np.in1d(all_telem['obsid'], stats[stats['obsid_ok']]['obsid'])

    f_ok = sum(ok)/len(ok)

    ok *= (mags < 13.9)

    star = get_star(agasc_id, date=all_telem['times'][0])
    result = {
        'agasc_id': agasc_id,
        'mag_aca': star['MAG_ACA'],
        'mag_aca_err': star['MAG_ACA_ERR']/100,
        'color': star['COLOR1'],
        'n_obsids': n_obsids,
        'n_obsids_ok': sum(stats['obsid_ok']),
        'n': len(ok),
        'n_ok': sum(ok),
        'f_ok': f_ok,  # f_ok does not count samples with mag >= 13.9
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
        't_mean_2': 0,
        't_std_2': 0,
        'n_outlier_2': 0,
    }
    for dr in [3, 5]:
        result.update({
            f't_mean_dr{dr}': 0,
            f't_std_dr{dr}': 0,
            f'mean_dr{dr}': 0,
            f'std_dr{dr}': 0,
            f'f_dr{dr}': 0,
            f'n_dr{dr}': 0,
            f'n_dr{dr}_outliers': 0,
            f'median_dr{dr}': 0,
            f'sigma_minus_dr{dr}': 0,
            f'sigma_plus_dr{dr}': 0,
        })

    if sum(ok) < 10:
        return result, stats

    sigma_minus, q25, median, q75, sigma_plus = np.quantile(mags[ok],
                                                            [0.158, 0.25, 0.5, 0.75, 0.842])
    iqr = q75 - q25
    outlier_1 = ok & ((mags > q75 + 1.5 * iqr) | (mags < q25 - 1.5 * iqr))
    outlier_2 = ok & ((mags > q75 + 3 * iqr) | (mags < q25 - 3 * iqr))

    # combine measurements using a weighted mean
    min_std = max(0.1, stats['std'].min())
    stats['w'] = np.where(stats['std'] != 0,
                          1. / stats['std'],
                          1. / min_std)
    stats['weighted_mean'] = stats['mean_corrected'] * stats['w']

    mag_weighted_mean = (stats['weighted_mean'].sum() / stats['w'].sum())
    mag_weighted_std = (
        np.sqrt(((stats['mean'] - mag_weighted_mean)**2 * stats['w']).sum() / stats['w'].sum())
    )

    result.update({
        'agasc_id': agasc_id,
        'n_obsids': len(stats),
        'n_obsids_ok': sum(stats['obsid_ok']),
        'n': len(ok),
        'n_ok': sum(ok),
        'f_ok': f_ok,  # f_ok does not count samples with mag >= 13.9
        'median': median,
        'sigma_minus': sigma_minus,
        'sigma_plus': sigma_plus,
        'mean': np.mean(mags[ok]),
        'std': np.std(mags[ok]),
        'mag_weighted_mean': mag_weighted_mean,
        'mag_weighted_std': mag_weighted_std,
        't_mean': np.mean(mags[ok & (~outlier_1)]),
        't_std': np.std(mags[ok & (~outlier_1)]),
        'n_outlier': sum(ok & outlier_1),
        't_mean_2': np.mean(mags[ok & (~outlier_2)]),
        't_std_2': np.std(mags[ok & (~outlier_2)]),
        'n_outlier_2': sum(ok & outlier_2),
    })

    for dr in [3, 5]:
        k = ok & (all_telem['dr'] < dr)
        sigma_minus, q25, median, q75, sigma_plus = np.quantile(mags[k],
                                                                [0.158, 0.25, 0.5, 0.75, 0.842])
        outlier = ok & ((mags > q75 + 1.5 * iqr) | (mags < q25 - 1.5 * iqr))
        result.update({
            f't_mean_dr{dr}': np.mean(mags[k & (~outlier)]),
            f't_std_dr{dr}': np.std(mags[k & (~outlier)]),
            f'mean_dr{dr}': np.mean(mags[k]),
            f'std_dr{dr}': np.std(mags[k]),
            f'f_dr{dr}': sum(k) / sum(ok),
            f'n_dr{dr}': sum(k),
            f'n_dr{dr}_outliers': sum(k & outlier),
            f'median_dr{dr}': median,
            f'sigma_minus_dr{dr}': sigma_minus,
            f'sigma_plus_dr{dr}': sigma_plus,
        })

    return result, stats
