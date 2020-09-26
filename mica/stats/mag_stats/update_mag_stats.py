#!/usr/bin/env python
import os
import pickle
import yaml
import numpy as np
import tables
import datetime
from functools import partial
from multiprocessing import Pool
from astropy import table

from mica.stats.mag_stats import catalogs, mag_stats, mag_stats_report as msr
from cxotime import CxoTime
from astropy import time, units as u


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
        return CxoTime(t_stop).date, CxoTime(t_start).date


def get_agasc_id_stats(agasc_ids, obs_status_override={}, tstop=None):
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
                                             obs_status_override=obs_status_override,
                                             tstop=tstop)
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


def get_agasc_id_stats_pool(agasc_ids, obs_status_override={}, batch_size=100, tstop=None):
    """
    Call update_mag_stats.get_agasc_id_stats multiple times using a multiprocessing.Pool

    :param agasc_ids: list
    :param batch_size: int
    :return: astropy.table.Table, astropy.table.Table, list
        obsid_stats, agasc_stats, fails, failed_jobs
    """
    import time
    from astropy.table import vstack, Table

    fmt = '%Y-%m-%d %H:%M'
    jobs = []
    n = len(agasc_ids)
    args = []
    progress = 0
    finished = 0
    for i in range(0, n, batch_size):
        args.append(agasc_ids[i:i + batch_size])
    with Pool() as pool:
        for arg in args:
            jobs.append(pool.apply_async(get_agasc_id_stats,
                                         [arg, obs_status_override, tstop]))
        start = datetime.datetime.now()
        now = None
        while finished < len(jobs):
            finished = sum([f.ready() for f in jobs])
            if now is None or 100*finished/len(jobs) - progress > 0.02:
                now = datetime.datetime.now()
                if finished == 0:
                    eta = ''
                else:
                    dt1 = (now - start).total_seconds()
                    dt = datetime.timedelta(seconds=(len(jobs)-finished) * dt1 / finished)
                    eta = f'ETA: {(now + dt).strftime(fmt)}'
                progress = 100*finished/len(jobs)
                print(f'{progress:6.2f}% at {now.strftime(fmt)}, {eta}')
            time.sleep(1)
    fails = []
    failed_agasc_ids = [i for arg, job in zip(args, jobs) if not job.successful() for i in arg]
    for agasc_id in failed_agasc_ids:
        fails.append(dict(mag_stats.MagStatsException(agasc_id=agasc_id, msg='Failed job')))

    results = [job.get() for job in jobs if job.successful()]

    # TODO: make sure nothing is skipped here
    obsid_stats = [r[0] for r in results if r[0] is not None]
    agasc_stats = [r[1] for r in results if r[1] is not None]
    obsid_stats = vstack(obsid_stats) if obsid_stats else Table()
    agasc_stats = vstack(agasc_stats) if agasc_stats else Table()
    fails += sum([r[2] for r in results], [])

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


def update_supplement(agasc_stats, filename, obs_status=None, include_all=True):
    """
    Update the magnitude table of the AGASC supplement.

    :param agasc_stats:
    :param filename:
    :return:
    """
    if include_all:
        outliers_new = agasc_stats[
            (agasc_stats['n_obsids_ok'] > 0)
        ]
    else:
        outliers_new = agasc_stats[
            (agasc_stats['n_obsids_ok'] > 0) &
            (agasc_stats['selected_atol'] |
             agasc_stats['selected_rtol'] |
             agasc_stats['selected_color'] |
             agasc_stats['selected_mag_aca_err'])
        ]
    outliers_new['mag_aca'] = outliers_new['mag_obs']
    outliers_new['mag_aca_err'] = outliers_new['mag_obs_err']
    names = ['agasc_id', 'color', 'mag_aca', 'mag_aca_err', 'last_obs_time']
    outliers_new = outliers_new[names].as_array()

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

    t = list(zip(*[[oi, ai, obs_status[(oi, ai)]['ok'], obs_status[(oi, ai)]['comments']]
                   for oi, ai in obs_status]))
    if t:
        obs_status = table.Table(t, names=['obsid', 'agasc_id', 'ok', 'comments']).as_array()
    else:
        obs_status = table.Table(dtype=[('obsid', int),
                                        ('agasc_id', int),
                                        ('ok', bool),
                                        ('comments', '<U80')]).as_array()

    mode = 'r+' if os.path.exists(filename) else 'w'
    with tables.File(filename, mode) as h5:
        if 'mags' in h5.root:
            h5.remove_node('/mags')
        h5.create_table('/', 'mags', outliers)
        if 'obs_status' in h5.root:
            h5.remove_node('/obs_status')
        h5.create_table('/', 'obs_status', obs_status)

    return new_stars, updated_stars


def do(args):
    filename = f'agasc_supplement_{mag_stats.version}.h5'

    if args.multi_process:
        get_stats = partial(get_agasc_id_stats_pool, batch_size=10)
    else:
        get_stats = get_agasc_id_stats
    catalogs.load(args.stop)

    # first, get list of AGASC IDs from file, from start/stop or take all observations.
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

    # default values for start/stop cover all the observations
    if args.start is None:
        args.start = CxoTime(stars_obs['mp_starcat_time']).min().date
    if args.stop is None:
        args.stop = CxoTime(stars_obs['mp_starcat_time']).max().date

    # exclude/include an ad-hoc list of observations
    obs_status_override = {}
    if args.obs_status_override:
        with open(args.obs_status_override) as f:
            obs_status_override = yaml.load(f, Loader=yaml.FullLoader)
        for obs in obs_status_override:
            if 'agasc_id' not in obs_status_override[obs]:
                obs_status_override[obs]['agasc_id'] = list(
                    catalogs.STARS_OBS[catalogs.STARS_OBS['obsid'] == obs]['agasc_id'])
        obs_status_override = {
            (obs, agasc_id): {
                'ok': obs_status_override[obs]['ok'],
                'comments': obs_status_override[obs]['comments']
            }
            for obs in obs_status_override for agasc_id in obs_status_override[obs]['agasc_id']}

    # if supplement exists, get the latest observation for each agasc_id in stars_obs,
    # find the ones already in the supplement, and drop the ones for which supplement.last_obs_time
    # is equal or larger than stars_obs.mp_starcat_time
    obs_status = {}
    if os.path.exists(filename):
        with tables.File(filename, 'r+') as h5:
            outliers_current = h5.root.mags[:]
            times = stars_obs[['agasc_id', 'mp_starcat_time']].group_by(
                'agasc_id').groups.aggregate(lambda d: np.max(CxoTime(d)).date)
            times = table.join(times,
                               table.Table(outliers_current[['agasc_id', 'last_obs_time']]),
                               join_type='left')
            if hasattr(times['last_obs_time'], 'mask'):
                # the mask can be there if there are new stars
                update = (times['last_obs_time'].mask | (
                          (~times['last_obs_time'].mask) &
                          (CxoTime(times['mp_starcat_time']).cxcsec > times['last_obs_time']).data))
            else:
                update = (CxoTime(times['mp_starcat_time']).cxcsec > times['last_obs_time'])

            # also update stars passed in the obs_status_override
            override_agasc_ids = np.in1d(times['agasc_id'], [ai for _, ai in obs_status_override])
            if np.any(override_agasc_ids):
                update += override_agasc_ids
                print(f'Not skipping {np.sum(override_agasc_ids)} stars overridden explicitly')

            stars_obs = stars_obs[np.in1d(stars_obs['agasc_id'], times[update]['agasc_id'])]
            agasc_ids = np.sort(np.unique(stars_obs['agasc_id']))
            print(f'Skipping {len(update) - np.sum(update)} stars (already in the supplement)')

            if 'obs' in h5.root:
                obs_status = {
                    (r['obsid'], r['agasc_id']): {'ok': r['ok'], 'comments': r['comments']}
                    for r in h5.root.obs_status[:]
                }

    obs_status.update(obs_status_override)

    # do the processing
    print(f'Will process {len(agasc_ids)} stars on {len(stars_obs)} observations')
    obsid_stats, agasc_stats, fails = \
        get_stats(agasc_ids, tstop=args.stop, obs_status_override=obs_status)

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
        new_stars, updated_stars = update_supplement(agasc_stats,
                                                     filename=filename, obs_status=obs_status)

        print(f'  {len(new_stars)} new stars, {len(updated_stars)} updated stars')
        if args.report:
            now = datetime.datetime.now()
            print(f"making report at {now}")
            sections = [{
                'id': 'new_stars',
                'title': 'New Stars',
                'stars': new_stars
            }, {
                'id': 'updated_stars',
                'title': 'Updated Stars',
                'stars': updated_stars['agasc_id'] if len(updated_stars) else []
            }]

            week = time.TimeDelta(7 * u.day)
            t = CxoTime(args.stop)
            nav_links = {
                'previous': f'../{(t - week).date[:8]}/index.html',
                'up': '..',
                'next': f'../{(t + week).date[:8]}/index.html'
            }
            report = msr.MagStatsReport(agasc_stats, obsid_stats,
                                        directory=os.path.join('weekly_reports', f'{t.date[:8]}'))
            report.multi_star_html(filename='index.html',
                                   sections=sections,
                                   updated_stars=updated_stars,
                                   fails=fails,
                                   report_date=CxoTime.now().date,
                                   tstart=args.start,
                                   tstop=args.stop,
                                   nav_links=nav_links,
                                   include_all_stars=True)
    now = datetime.datetime.now()
    print(f"done at {now}")
