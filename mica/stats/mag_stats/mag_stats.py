import scipy.stats
import numpy as np
from astropy.table import Table, vstack

import catalogs

import pandas as pd
from Chandra.Time import DateTime
from agasc import get_star
from cheta import fetch
from Quaternion import Quat
import Ska.quatutil


version = 'v0.0.12'


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


def get_residuals(star, slot, times):
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
    names = [f'AOATTQT{i}' for i in range(1, 5)] + [f'AOACZAN{slot}', f'AOACYAN{slot}']

    msids = fetch.MSIDset(names, times[0] - 4, times[-1] + 4)
    t = np.in1d(msids[names[0]].times, times)
    telem = {n: msids[n].vals[t] for n in names}

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
    dy = telem[f'AOACYAN{slot}'] - yag
    dz = telem[f'AOACZAN{slot}'] - zag

    # telem['yag'] = yag
    # telem['zag'] = zag
    telem['AOACYAN'] = telem[f'AOACYAN{slot}']
    telem['AOACZAN'] = telem[f'AOACZAN{slot}']

    del telem[f'AOACYAN{slot}']
    del telem[f'AOACZAN{slot}']

    telem['dy'] = dy
    telem['dz'] = dz
    # cheating here and ignoring spherical trig
    telem['dr'] = (dy ** 2 + dz ** 2) ** .5
    return telem


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
    start = dwell['tstart']
    stop = dwell['tstop']
    slot = obs['slot']

    dmag = _magnitude_correction(start, obs['mag'])

    msid = fetch.MSID(f'AOACMAG{slot}', start, stop)
    times = msid.times[::2]
    mags = msid.vals[::2]

    telem = {
        'times': times,
        'mags': mags,
        'mags_corrected': mags - dmag
    }

    names = ['AOACASEQ', 'AOPCADMD', f'AOACIIR{slot}', f'AOACISP{slot}']
    msids = fetch.MSIDset(names, times[0] - 4, times[-1] + 4)
    # the following works because the MSIDs are in the same group
    t = np.in1d(msids[names[0]].times, times)
    telem.update({n: msids[n].vals[t] for n in names})

    telem['AOACIIR'] = telem[f'AOACIIR{slot}']
    telem['AOACISP'] = telem[f'AOACISP{slot}']

    del telem[f'AOACIIR{slot}']
    del telem[f'AOACISP{slot}']

    star = get_star(obs['agasc_id'], date=times[0])
    telem.update(get_residuals(star=star, slot=obs['slot'], times=times))

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
            (telem[f'AOACISP'] == 'OK ') & (telem['AOPCADMD'] == 'NPNT')
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
    s_100s = s.rolling(window=int(100 / dt), center=True).median().dropna()
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
        all_telem.append(telem)

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

    mags = all_telem['mags_corrected']

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

    if sum(ok) < 100:
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
