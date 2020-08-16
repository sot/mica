#!/usr/bin/env python
import datetime
from functools import partial
from multiprocessing import Pool
from mica.stats.mag_stats import update_mag_stats, mag_stats


def get_agasc_id_stats(agasc_ids, excluded_observations={}, batch_size=100):
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
    with Pool(processes=12) as pool:
        for arg in args:
            jobs.append(pool.apply_async(update_mag_stats.get_agasc_id_stats,
                                         [arg, excluded_observations]))
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


def main():
    get_stats = partial(get_agasc_id_stats, batch_size=10)
    update_mag_stats.do(get_stats)


if __name__ == '__main__':
    import warnings
    warnings.simplefilter('ignore', UserWarning)
    main()
