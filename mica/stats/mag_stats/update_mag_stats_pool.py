#!/usr/bin/env python
import numpy as np
import datetime
from multiprocessing import Pool
from mica.stats.mag_stats import update_mag_stats, catalogs, mag_stats
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
        while finished < len(jobs):
            finished = sum([f.ready() for f in jobs])
            if 100*finished/len(jobs) - progress > 0.1:
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

    return obsid_stats, agasc_stats, fails, failed_jobs


def parser():
    parse = argparse.ArgumentParser()
    parse.add_argument('--agasc-id-file')
    return parse

def main():
    args = parser().parse_args()
    if args.agasc_id_file:
        with open(args.agasc_id_file, 'r') as f:
            agasc_ids = [int(l.strip()) for l in f.readlines()]
            agasc_ids = np.intersect1d(agasc_ids, catalogs.STARS_OBS['agasc_id'])
    else:
        agasc_ids = sorted(catalogs.STARS_OBS['agasc_id'])
    agasc_ids = np.unique(agasc_ids)
    stars_obs = catalogs.STARS_OBS[np.in1d(catalogs.STARS_OBS['agasc_id'], agasc_ids)]
    batch_size = 10
    print(f'Will process {len(agasc_ids)} stars'
          f' on {len(stars_obs)} observations')

    obsid_stats, agasc_stats, fails, fail_jobs = get_agasc_id_stats(agasc_ids, batch_size)

    print(f'Got:\n'
          f'  {len(obsid_stats)} OBSIDs,'
          f'  {len(agasc_stats)} stars,'
          f'  {len(fails)} failed stars,'
          f'  {len(fail_jobs)} failed jobs')
    filename = f'mag_stats_failed_jobs_{mag_stats.version}.pkl'
    with open(filename, 'wb') as out:
        pickle.dump(fail_jobs, out)
    if len(agasc_stats):
        update_mag_stats.update_mag_stats(obsid_stats, agasc_stats,
                                          {'fails': fails, 'fail_jobs': fail_jobs})



if __name__ == '__main__':
    main()