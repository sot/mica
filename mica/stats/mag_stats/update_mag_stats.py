#!/usr/bin/env python
import os
import numpy as np
from mica.stats.mag_stats import catalogs, mag_stats
import pickle
import tables

np.seterr(all='ignore')


def level0_archive_time_range():
    import sqlite3
    import os
    db_file = os.path.expandvars('$SKA/data/mica/archive/aca0/archfiles.db3')
    with sqlite3.connect(db_file) as connection:
        cursor = connection.cursor()
        cursor.execute("select tstop from archfiles order by tstop desc limit 1")
        t_stop = cursor.fetchall()[0][0]
        cursor.execute("select tstop from archfiles order by tstart asc limit 1")
        t_start = cursor.fetchall()[0][0]
        return t_stop, t_start


def get_agasc_id_stats(agasc_ids):
    """
    Call mag_stats.get_agasc_id_stats for each AGASC ID

    :param agasc_ids: list
    :return: astropy.table.Table, astropy.table.Table, list
        obsid_stats, agasc_stats, fails
    """
    from mica.stats.mag_stats import mag_stats
    from astropy.table import Table, vstack

    fails = []
    obsid_stats = []
    agasc_stats = []
    obs_failures = []
    for i, agasc_id in enumerate(agasc_ids):
        try:
            agasc_stat, obsid_stat, obs_fail = mag_stats.get_agasc_id_stats(agasc_id=agasc_id)
            agasc_stats.append(agasc_stat)
            obsid_stats.append(obsid_stat)
            obs_failures.append(obs_fail)
        except Exception as e:
            fails.append((agasc_id, str(e)))

    obs_failures = [f for obs in obs_failures for f in obs]
    try:
        agasc_stats = Table(agasc_stats) if agasc_stats else None
        obsid_stats = vstack(obsid_stats) if obsid_stats else None
    except Exception as e:
        agasc_stats = None
        obsid_stats = None
        fails = [(agasc_id, f'Exception at bottom of get_agasc_id_stats: {str(e)}')]
    return obsid_stats, agasc_stats, fails, obs_failures


def update_mag_stats(obsid_stats, agasc_stats, fails, outdir='.'):
    """
    Update the mag_stats catalog.

    I currently save three files:
    - mag_stats_agasc_{mag_stats.version}.fits with stats for each AGASC ID
    - mag_stats_obsid_{mag_stats.version}.fits with stats for each OBSID
    - mag_stats_fails_{mag_stats.version}.pkl with a list of failures

    :param obsid_stats:
    :param agasc_stats:
    :param fails:
    :param outdir:
    :return:
    """
    if len(agasc_stats):
        filename = os.path.join(outdir, f'mag_stats_agasc_{mag_stats.version}.fits')
        if os.path.exists(filename):
            os.remove(filename)
        agasc_stats.write(filename)
    if len(obsid_stats):
        filename = os.path.join(outdir, f'mag_stats_obsid_{mag_stats.version}.fits')
        if os.path.exists(filename):
            os.remove(filename)
        obsid_stats.write(filename)
    if len(fails):
        filename = os.path.join(outdir, f'mag_stats_fails_{mag_stats.version}.pkl')
        with open(filename, 'wb') as out:
            pickle.dump(fails, out)


def update_mags_stats_supplement(agasc_stats, filename=None):
    """
    Update the magnitude table of the AGASC supplement.

    :param agasc_stats:
    :param filename:
    :return:
    """

    agasc_stats['outlier'] = \
        np.abs(agasc_stats['t_mean_dr3'] - agasc_stats['mag_aca']) > 3 * agasc_stats['mag_aca_err']
    agasc_stats['mag_aca'] = agasc_stats['t_mean_dr3']
    agasc_stats['mag_aca_err'] = agasc_stats['t_std_dr3']

    outliers = agasc_stats[
        (agasc_stats['color'] == 1.5) | (agasc_stats['color'] == 0.7) | agasc_stats['outlier']]
    outliers = outliers[['agasc_id', 'color', 'mag_aca', 'mag_aca_err']].as_array()

    if filename is not None:
        filename = f'agasc_supplement_{mag_stats.version}.h5'
    mode = 'r+' if os.path.exists(filename) else 'w'
    with tables.File(filename, mode) as h5:
        if 'mags' in h5.root:
            h5.remove_node('/mags')
        h5.create_table('/', 'mags', outliers)


def main():
    agasc_ids = sorted(np.unique(catalogs.STARS_OBS['agasc_id']))
    print(f'Will process {len(agasc_ids)} stars on {len(catalogs.STARS_OBS)} observations')
    obsid_stats, agasc_stats, fails, obs_failures = get_agasc_id_stats(agasc_ids)
    print(f'Got:\n'
          f'  {len(obsid_stats)} OBSIDs,'
          f'  {len(agasc_stats)} stars,'
          f'  {len(fails)} failed stars,'
          f'  {len(obs_failures)} failed observations')
    update_mag_stats(obsid_stats, agasc_stats, fails)


if __name__ == '__main__':
    main()
