#!/usr/bin/env python
import os
import pickle
import argparse
import yaml
import numpy as np
import tables
from astropy import table

from mica.stats.mag_stats import catalogs, mag_stats, mag_stats_report as msr
from cxotime import CxoTime


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


def get_agasc_id_stats(agasc_ids, excluded_observations={}):
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
    for i, agasc_id in enumerate(agasc_ids):
        try:
            agasc_stat, obsid_stat, obs_fail = \
                mag_stats.get_agasc_id_stats(agasc_id=agasc_id,
                                             excluded_observations=excluded_observations)
            agasc_stats.append(agasc_stat)
            obsid_stats.append(obsid_stat)
            fails += obs_fail
        except mag_stats.MagStatsException as e:
            fails.append(dict(e))
        except Exception as e:
            # transform Exception to MagStatsException for standard book keeping
            fails.append(dict(mag_stats.MagStatsException(agasc_id=agasc_id, msg=str(e))))

    try:
        agasc_stats = Table(agasc_stats) if agasc_stats else None
        obsid_stats = vstack(obsid_stats) if obsid_stats else None
    except Exception as e:
        agasc_stats = None
        obsid_stats = None
        # transform Exception to MagStatsException for standard book keeping
        fails.append(dict(mag_stats.MagStatsException(
            msg=f'Exception at end of get_agasc_id_stats: {str(e)}')))

    return obsid_stats, agasc_stats, fails


def _update_table(table_old, table_new, keys):
    # checking names, because actual types change upon saving in fits format
    assert table_old.as_array().dtype.names == table_new.as_array().dtype.names, \
        'Tables have different dtype'
    table_old = table_old.copy()
    new_row = np.ones(len(table_new), dtype=bool)
    _, i_new, i_old = np.intersect1d(table_new[keys].as_array(),
                                     table_old[keys].as_array(),
                                     return_indices=True)
    new_row[i_new] = False
    table_old[i_old] = table_new[i_new]
    return table.vstack([table_old, table_new[new_row]])


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
    if agasc_stats is not None and len(agasc_stats):
        filename = os.path.join(outdir, f'mag_stats_agasc_{mag_stats.version}.fits')
        if os.path.exists(filename):
            agasc_stats = _update_table(table.Table.read(filename), agasc_stats,
                                        keys=['agasc_id'])
            os.remove(filename)
        agasc_stats.write(filename)
    if obsid_stats is not None and len(obsid_stats):
        filename = os.path.join(outdir, f'mag_stats_obsid_{mag_stats.version}.fits')
        if os.path.exists(filename):
            obsid_stats = _update_table(table.Table.read(filename), obsid_stats,
                                        keys=['agasc_id', 'obsid', 'timeline_id'])
            os.remove(filename)
        obsid_stats.write(filename)
    if len(fails):
        filename = os.path.join(outdir, f'mag_stats_fails_{mag_stats.version}.pkl')
        with open(filename, 'wb') as out:
            pickle.dump(fails, out)


def update_supplement(agasc_stats, filename=None):
    """
    Update the magnitude table of the AGASC supplement.

    :param agasc_stats:
    :param filename:
    :return:
    """
    outliers_new = agasc_stats[
        (agasc_stats['n_obsids_ok'] > 0) &
        (agasc_stats['selected_atol'] |
         agasc_stats['selected_rtol'] |
         agasc_stats['selected_color'] |
         agasc_stats['selected_mag_aca_err'])
    ]
    outliers_new['mag_aca'] = outliers_new['t_mean_dr3']
    outliers_new['mag_aca_err'] = outliers_new['t_std_dr3']
    names = ['agasc_id', 'color', 'mag_aca', 'mag_aca_err', 'last_obs_time']
    outliers_new = outliers_new[names].as_array()

    if filename is None:
        filename = f'agasc_supplement_{mag_stats.version}.h5'
    if os.path.exists(filename):
        # I could do what follows directly in place, but the table is not that large.
        with tables.File(filename, 'r+') as h5:
            outliers_current = h5.root.mags[:]
            # find the indices of agasc_ids in both current and new lists
            _, i_new, i_cur = np.intersect1d(outliers_new['agasc_id'],
                                             outliers_current['agasc_id'],
                                             return_indices=True)
            current = outliers_current[i_cur]
            new = outliers_new[i_new]

            # from those, find the ones which differ in last observation time
            i_cur = i_cur[current['last_obs_time'] != new['last_obs_time']]
            i_new = i_new[current['last_obs_time'] != new['last_obs_time']]
            # overwrite current values with new values (and calculate diff to return)
            updated_stars = np.zeros(len(outliers_new[i_new]),
                                     dtype=[('agasc_id', np.int64),
                                            ('mag_aca', np.float64),
                                            ('mag_aca_err', np.float64)])
            updated_stars['mag_aca'] = (outliers_new[i_new]['mag_aca'] -
                                        outliers_current[i_cur]['mag_aca'])
            updated_stars['mag_aca_err'] = (outliers_new[i_new]['mag_aca_err'] -
                                            outliers_current[i_cur]['mag_aca_err'])
            updated_stars['agasc_id'] = outliers_new[i_new]['agasc_id']
            outliers_current[i_cur] = outliers_new[i_new]

            # find agasc_ids in new list but not in current list
            new_stars = ~np.in1d(outliers_new['agasc_id'], outliers_current['agasc_id'])
            # and add them to the current list
            outliers_current = np.concatenate([outliers_current, outliers_new[new_stars]])
            outliers = np.sort(outliers_current)

            new_stars = outliers_new[new_stars]['agasc_id']
    else:
        outliers = outliers_new
        new_stars = outliers_new['agasc_id']
        updated_stars = np.array([], dtype=[('agasc_id', np.int64), ('mag_aca', np.float64),
                                            ('mag_aca_err', np.float64)])

    mode = 'r+' if os.path.exists(filename) else 'w'
    with tables.File(filename, mode) as h5:
        if 'mags' in h5.root:
            h5.remove_node('/mags')
        h5.create_table('/', 'mags', outliers)

    return new_stars, updated_stars


def parser():
    parse = argparse.ArgumentParser()
    parse.add_argument('--agasc-id-file')
    parse.add_argument('--start')
    parse.add_argument('--stop')
    parse.add_argument('--excluded-observations')
    parse.add_argument('--report', action='store_true', default=False)
    return parse


def do(get_stats=get_agasc_id_stats):
    args, _ = parser().parse_known_args()
    catalogs.load(args.stop)
    if args.agasc_id_file:
        with open(args.agasc_id_file, 'r') as f:
            agasc_ids = [int(l.strip()) for l in f.readlines()]
            agasc_ids = np.intersect1d(agasc_ids, catalogs.STARS_OBS['agasc_id'])
    elif args.start:
        if not args.stop:
            args.stop = CxoTime.now().date
        else:
            args.stop = CxoTime(args.stop).date
        args.start = CxoTime(args.start).date
        obs_in_time = ((catalogs.STARS_OBS['mp_starcat_time'] >= args.start) &
                       (catalogs.STARS_OBS['mp_starcat_time'] <= args.stop))
        agasc_ids = sorted(catalogs.STARS_OBS[obs_in_time]['agasc_id'])
    else:
        agasc_ids = sorted(catalogs.STARS_OBS['agasc_id'])
    agasc_ids = np.unique(agasc_ids)
    stars_obs = catalogs.STARS_OBS[np.in1d(catalogs.STARS_OBS['agasc_id'], agasc_ids)]

    if args.start is None:
        args.start = CxoTime(stars_obs['mp_starcat_time']).min().date
    if args.stop is None:
        args.stop = CxoTime(stars_obs['mp_starcat_time']).max().date

    excluded_observations = {}
    if args.excluded_observations:
        with open(args.excluded_observations) as f:
            exclude = yaml.load(f, Loader=yaml.FullLoader)
            for k in exclude:
                for obsid in exclude[k]:
                    if obsid not in excluded_observations:
                        excluded_observations[obsid] = ''
                    excluded_observations[obsid] += k

    print(f'Will process {len(agasc_ids)} stars on {len(stars_obs)} observations')
    obsid_stats, agasc_stats, fails = \
        get_stats(agasc_ids, excluded_observations=excluded_observations)

    failed_global = [f for f in fails if not f['agasc_id'] and not f['obsid']]
    failed_stars = [f for f in fails if f['agasc_id'] and not f['obsid']]
    failed_obs = [f for f in fails if f['obsid']]
    print(f'Got:\n'
          f'  {len(obsid_stats)} OBSIDs,'
          f'  {len(agasc_stats)} stars,'
          f'  {len(failed_stars)} failed stars,'
          f'  {len(failed_obs)} failed observations,'
          f'  {len(failed_global)} global errors')
    if len(agasc_stats):
        update_mag_stats(obsid_stats, agasc_stats, fails)
        new_stars, updated_stars = update_supplement(agasc_stats)

        print(f'{len(new_stars)} new stars, {len(updated_stars)} updated stars')
        if args.report:
            print("making report")
            msr.multi_star_html_report(agasc_stats, obsid_stats, new_stars, updated_stars,
                                       fails=fails, report_date=CxoTime.now().date,
                                       tstart=args.start, tstop=args.stop)

def main():
    do()


if __name__ == '__main__':
    main()
