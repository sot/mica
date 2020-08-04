#!/usr/bin/env python
import numpy as np
import datetime
from multiprocessing import Pool
from mica.stats.mag_stats import update_mag_stats, catalogs, mag_stats
from cxotime import CxoTime
import pickle
import argparse


def get_agasc_id_stats(agasc_ids, batch_size=100):
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
            jobs.append(pool.apply_async(update_mag_stats.get_agasc_id_stats, [arg]))
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
    failed_jobs = [arg for arg, job in zip(args, jobs) if not job.successful()]
    results = [job.get() for job in jobs if job.successful()]

    # TODO: make sure nothing is skipped here
    obsid_stats = [r[0] for r in results if r[0] is not None]
    agasc_stats = [r[1] for r in results if r[1] is not None]
    obsid_stats = vstack(obsid_stats) if obsid_stats else Table()
    agasc_stats = vstack(agasc_stats) if agasc_stats else Table()
    fails = sum([r[2] for r in results], [])
    obs_fails = sum([r[3] for r in results], [])

    return obsid_stats, agasc_stats, fails, failed_jobs, obs_fails


def parser():
    parse = argparse.ArgumentParser()
    parse.add_argument('--agasc-id-file')
    parse.add_argument('--start')
    parse.add_argument('--stop')
    return parse


def main():
    args = parser().parse_args()
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

    batch_size = 10
    print(f'Will process {len(agasc_ids)} stars'
          f' on {len(stars_obs)} observations')

    obsid_stats, agasc_stats, fails, fail_jobs, fail_obs = get_agasc_id_stats(agasc_ids, batch_size)

    print(f'Got:\n'
          f'  {len(obsid_stats)} OBSIDs,'
          f'  {len(agasc_stats)} stars,'
          f'  {len(fails)} failed stars,'
          f'  {len(fail_jobs)} failed jobs'
          f'  {len(fail_obs)} failed observations')
    filename = f'mag_stats_failed_jobs_{mag_stats.version}.pkl'
    with open(filename, 'wb') as out:
        pickle.dump(fail_jobs, out)
    if len(agasc_stats):
        update_mag_stats.update_mag_stats(obsid_stats, agasc_stats,
                                          {'fails': fails,
                                           'fail_jobs': fail_jobs,
                                           'fail_obs': fail_obs})

        new_stars, updated_stars = update_mag_stats.update_supplement(agasc_stats)


if __name__ == '__main__':
    import warnings
    warnings.simplefilter('ignore', UserWarning)
    main()
