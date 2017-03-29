import os
import re
import logging
import argparse
from itertools import count
from operator import itemgetter
import numpy as np

from six.moves import zip
from Chandra.Time import DateTime
from Chandra.cmd_states import get_cmd_states
import Ska.Shell
import Ska.DBI

from mica.starcheck.starcheck_parser import read_starcheck
from mica.common import MICA_ARCHIVE

logger = logging.getLogger('starcheck ingest')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

DEFAULT_CONFIG = dict(starcheck_db=dict(dbi='sqlite',
                                        server=os.path.join(MICA_ARCHIVE, 'starcheck', 'starcheck.db3')),
                      mp_top_level='/data/mpcrit1/mplogs',
                      timelines_db=dict(dbi='sqlite',
                                        server=os.path.join(os.environ['SKA'], 'data', 'cmd_states', 'cmd_states.db3')))
FILES = dict(data_root=os.path.join(MICA_ARCHIVE, 'starcheck'),
             sql_def='starcheck.sql')



def ingest_obs(obs, obs_idx, sc_id, st, db, existing=None):
    if existing is not None:
        existing_obs = [eobs for eobs in existing
                        if eobs['obsid'] == int(obs['obsid'])]
        if len(existing_obs):
            logger.debug(
                "Skipping ingest of %s %d" % (st, int(obs['obsid'])))
            return
    obs_d = obs['obs']
    obs_d.update(dict(sc_id=sc_id,
                      obs_idx=obs_idx))
    if obs_d.get('sim_z_offset_steps') is not None:
        obs_d['sim_z_offset_mm'] = (
            float(obs['obs'].get('sim_z_offset_steps'))
            * 1000
            * 2.51432e-06)
    logger.debug("inserting obs %s from %s" % (obs['obsid'], st))
    db.insert(obs_d, 'starcheck_obs')
    for manvr in obs['manvrs']:
        logger.debug("inserting manvr at %s" % manvr['mp_targquat_time'])
        db.insert(
            dict(
                sc_id=sc_id,
                obsid=obs['obsid'],
                obs_idx=obs_idx,
                **manvr),
            'starcheck_manvr')
    for star in obs['catalog']:
        star_d = dict(
            sc_id=sc_id,
            obsid=obs['obsid'],
            obs_idx=obs_idx,
            mp_starcat_time=obs['obs']['mp_starcat_time'],
            **star)
        logger.debug("inserting %s idx of catalog" % star['idx'])
        db.insert(star_d, 'starcheck_catalog')
    logger.debug("inserting %d warnings" % len(obs['warnings']))
    for warn in obs['warnings']:
        warn_d = dict(
            sc_id=sc_id,
            obsid=obs['obsid'],
            obs_idx=obs_idx,
            **warn)
        db.insert(warn_d, 'starcheck_warnings')
    if obs['pred_ccd_temp'] is not None:
        db.insert(dict(sc_id=sc_id, obsid=obs['obsid'], pred_ccd_temp=obs['pred_ccd_temp']),
                  'starcheck_pred_temp')
    db.commit()


def get_options():
    parser = argparse.ArgumentParser(
        description="Update starcheck database")
    defaults = dict(DEFAULT_CONFIG)
    parser.set_defaults(**defaults)
    parser.add_argument("--dbi",
                        help="dbi for starcheck database")
    parser.add_argument("--server",
                        help="server or file for database")
    parser.add_argument("--mp-top-level",
                        help="top level SOT MP dir")
    parser.add_argument("--start",
                        help="update database with starcheck files after this start time")
    opt = parser.parse_args()
    return opt


def prune_dirs(dirs, regex):
    """
    Prune directories (in-place) that do not match ``regex``.
    (borrowed from kadi)
    """
    prunes = [x for x in dirs if not re.match(regex, x)]
    for prune in prunes:
        dirs.remove(prune)


def get_new_starcheck_files(rootdir, mtime=0):
    """
    Look for starcheck.txt files in a a given top-level SOT MP directory
    and return those with modification times at or after the optional supplied 'mtime'.
    """
    logger.info("Getting new starcheck.txt files from {}".format(rootdir))
    starchecks = []
    for root, dirs, files in os.walk(rootdir):
        root = root.rstrip('/')
        depth = len(root.split('/')) - 1
        if depth == 3:
            prune_dirs(dirs, r'\d{4}$')
        elif depth == 4:
            prune_dirs(dirs, r'[A-Z]{3}\d{4}$')
        elif depth == 5:
            prune_dirs(dirs, r'ofls[a-z]$')
        elif depth > 5:
            files = [x for x in files if re.match(r'starcheck\.txt$', x)]
            if len(files):
                starchecks.append(os.path.join(root, files[0]))
            while dirs:
                dirs.pop()
    starchecks_with_times = [{'file': st, 'mtime': os.path.getmtime(st)}
                             for st in starchecks if os.path.getmtime(st) >= mtime]
    starchecks_with_times = sorted(starchecks_with_times, key=itemgetter('mtime'))
    return starchecks_with_times


def update(config=None):
    if config is None:
        config = DEFAULT_CONFIG
    if config['starcheck_db']['dbi'] != 'sqlite':
        raise ValueError("Only sqlite DBI implemented")
    last_starcheck_mtime = 0
    if (not os.path.exists(config['starcheck_db']['server'])
            or os.stat(config['starcheck_db']['server']).st_size == 0):
        if not os.path.exists(os.path.dirname(config['starcheck_db']['server'])):
            os.makedirs(os.path.dirname(config['starcheck_db']['server']))
        db_sql = os.path.join('mica', 'starcheck', FILES['sql_def'])
        db_init_cmds = open(db_sql).read()
        db = Ska.DBI.DBI(dbi='sqlite', server=config['starcheck_db']['server'])
        db.execute(db_init_cmds)
    else:
        db = Ska.DBI.DBI(dbi='sqlite', server=config['starcheck_db']['server'])
        max_mtime = db.fetchone("select max(mtime) as mtime from starcheck_id")
        if max_mtime is not None:
            last_starcheck_mtime = max_mtime['mtime']
    # If a start time is explicitly requested, override 0 or last database value
    if 'start' in config:
        last_starcheck_mtime = DateTime(config['start']).unix
    starchecks_with_times = get_new_starcheck_files(config['mp_top_level'],
                                                    mtime=last_starcheck_mtime)

    for st_with_time in starchecks_with_times:
        st = st_with_time['file']
        logger.info("Attempting ingest of %s" % st)
        # get an existing id or insert a new one
        mp_dir = re.search(
            '\/data\/mpcrit1\/mplogs(\/\d{4}\/\w{7}/\w{5}\/)starcheck.txt',
            st).group(1)
        idcheck = db.fetchone(
            "select id from starcheck_id where dir = '%s'" % mp_dir)
        existing = None
        if idcheck is not None:
            sc_id = idcheck['id']
            existing = db.fetchall(
                "select * from starcheck_obs where sc_id = '%d'" % sc_id)
        else:
            sc_id = db.fetchone(
                "select max(id) + 1 as id from starcheck_id")['id'] or 0
            db.insert(dict(id=sc_id, dir=mp_dir, mtime=st_with_time['mtime']), 'starcheck_id')
            db.commit()

        starcheck = read_starcheck(st)
        for (obs, obs_idx) in zip(starcheck, count(0)):
            ingest_obs(obs, obs_idx, sc_id, st, db, existing=existing)


def main():
    """
    Command line interface to update starcheck database from files
    in MP directories
    """
    opt = get_options()
    update(config=vars(opt))

if __name__ == '__main__':
    main()
