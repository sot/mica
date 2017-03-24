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
             touch_file=os.path.join(MICA_ARCHIVE, 'starcheck', "starcheck_parser.touch"),
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
    opt = parser.parse_args()
    return opt


def update(config=None):
    if config is None:
        config = DEFAULT_CONFIG
    if config['starcheck_db']['dbi'] != 'sqlite':
        raise ValueError("Only sqlite DBI implemented")
    if (not os.path.exists(config['starcheck_db']['server'])
            or os.stat(config['starcheck_db']['server']).st_size == 0):
        if not os.path.exists(os.path.dirname(config['starcheck_db']['server'])):
            os.makedirs(os.path.dirname(config['starcheck_db']['server']))
        db_sql = os.path.join(os.environ['SKA_DATA'], 'mica', FILES['sql_def'])
        db_init_cmds = open(db_sql).read()
        db = Ska.DBI.DBI(dbi='sqlite', server=config['starcheck_db']['server'])
        db.execute(db_init_cmds)
        if not os.path.exists(FILES['touch_file']):
            # make a touch file with a time before starchecks
            Ska.Shell.bash("touch -t 199801010000 %s" % (FILES['touch_file']),
                           env={'MAILCHECK': -1})
    else:
        db = Ska.DBI.DBI(dbi='sqlite', server=config['starcheck_db']['server'])

    starchecks = Ska.Shell.bash(
        "find %s" % config['mp_top_level']
        + "/[12]???/[A-Z][A-Z][A-Z][0-9][0-9][0-9][0-9]/ofls?/ "
        + " -maxdepth 1 -wholename '*/ofls?/starcheck.txt' "
        + "-newer %s" % FILES['touch_file'],
        env={'MAILCHECK': -1})

    starchecks_with_times = [dict(file=st, mtime=os.path.getmtime(st))
                             for st in starchecks]
    for st_with_time in sorted(starchecks_with_times, key=itemgetter('mtime')):
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
            db.insert(dict(id=sc_id, dir=mp_dir), 'starcheck_id')
            db.commit()

        starcheck = read_starcheck(st)
        for (obs, obs_idx) in zip(starcheck, count(0)):
            ingest_obs(obs, obs_idx, sc_id, st, db, existing=existing)
        logger.info("Done with %s; updating touch file" % st)
        Ska.Shell.bash("touch -r %s %s" % (st, FILES['touch_file']),
                       env={'MAILCHECK': -1})


def main():
    """
    Command line interface to update starcheck database from files
    in MP directories
    """
    opt = get_options()
    update(config=vars(opt))

if __name__ == '__main__':
    main()
