import os
import sys
import re
import subprocess
import Ska.Shell
import json
import Ska.DBI
from itertools import izip, count
import logging
import argparse
from operator import itemgetter

logger = logging.getLogger('starcheck ingest')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

config = dict(data_root='/data/aca/archive/starcheck',
              touch_file="starcheck_parser.touch",
              dbi='sqlite',
              server='starcheck.db3',
              mp_top_level='/data/mpcrit1/mplogs')


def read_starcheck(starcheck_text_file):
    st = starcheck_text_file
    star_json, stderr = subprocess.Popen([
            'perl',
            '/proj/sot/ska/jeanproj/git/mica/mica/starcheck/starcheck_txt2json.pl',
            st],
                                         shell=False,
                                         stdout=subprocess.PIPE).communicate()
    if stderr:
        raise ValueError
    return json.loads(star_json)


def ingest_obs(obs, obs_idx, sc_id, st, db, existing=None):
    if 'obsid' not in obs['target']:
        if 'MP_STARCAT' in obs['times']:
            logger.warn(
                "No obsid for %s at %s" % (st, obs['times']['MP_STARCAT']))
        else:
            logger.warn(obs['warnings'])
        return
    if existing is not None:
        existing_obs = [eobs for eobs in existing
                        if eobs['obsid'] == int(obs['target']['obsid'])]
        if len(existing_obs):
            logger.debug(
                "Skipping ingest of %s %d" % (st, int(obs['target']['obsid'])))
            return
    obs_d = dict(
            sc_id=sc_id,
            obsid=obs['target']['obsid'],
            obs_idx=obs_idx,
            point_ra=obs['coords'].get('RA'),
            point_dec=obs['coords'].get('DEC'),
            point_roll=obs['coords'].get('ROLL'),
            target_id=obs['target'].get('target'),
            sci_instr=obs['target'].get('sci_instr'),
            sim_z_offset_steps=obs['target'].get('sim_z_offset_steps'),
            dither_state=obs['dither'].get('state'),
            dither_y_amp=obs['dither'].get('y_amp'),
            dither_y_period=obs['dither'].get('y_period'),
            dither_z_amp=obs['dither'].get('z_amp'),
            dither_z_period=obs['dither'].get('z_period'),
            mp_starcat_time=obs['times'].get('MP_STARCAT'),
            mp_starcat_vcdu_cnt=obs['times'].get('VCDU_cnt'),
            )
    if obs_d['sim_z_offset_steps'] is not None:
        obs_d['sim_z_offset_mm'] = (
            float(obs['target'].get('sim_z_offset_steps'))
            * 1000
            * 2.51432e-06)
    logger.debug("inserting obs %s from %s" % (obs_d['obsid'], st))
    db.insert(obs_d, 'starcheck_obs')
    for (manvr, instance) in izip(obs['manvr'], count(1)):
        logger.debug("inserting manvr at %s" % manvr['MP_TARGQUAT'])
        db.insert(dict(
                sc_id=sc_id,
                obsid=obs['target']['obsid'],
                obs_idx=obs_idx,
                instance=instance,
                duration_sec=manvr.get('duration_sec'),
                angle_deg=manvr.get('angle_deg'),
                slew_err_arcsec=manvr.get('slew_err_arcsec'),
                target_Q1=manvr['Q1'],
                target_Q2=manvr['Q2'],
                target_Q3=manvr['Q3'],
                target_Q4=manvr['Q4'],
                mp_targquat_time=manvr['MP_TARGQUAT'],
                mp_targquat_vcdu_cnt=manvr['VCDU_cnt']),
                  'starcheck_manvr')
    for star in obs['stars']:
        id_match = re.search('^(\D*)(\d+)', star['ID'])
        if id_match:
            s_id = id_match.group(2)
            s_idnote = id_match.group(1)
        else:
            if star['ID'] == '---':
                s_id = None
                s_idnote = star['ID']
            else:
                raise ValueError("id weird")
        star_d = dict(
                sc_id=sc_id,
                obsid=obs['target']['obsid'],
                obs_idx=obs_idx,
                mp_starcat_time=obs['times']['MP_STARCAT'],
                idx=star['IDX'],
                slot=star['SLOT'],
                id=s_id,
                idnote=s_idnote,
                type=star['TYPE'],
                sz=star['SZ'],
                minmag=star['MINMAG'],
                mag=star['MAG'],
                maxmag=star['MAXMAG'],
                yang=star['YANG'],
                zang=star['ZANG'],
                dim=star['DIM'],
                res=star['RES'],
                halfw=star['HALFW'],
                notes=star.get('NOTES'))
        star_d['pass'] = star.get('PASS')
        logger.debug("inserting %s idx of catalog" % star['IDX'])
        db.insert(star_d, 'starcheck_catalog')
    logger.debug("inserting %d warnings" % len(obs['warnings']))
    for warn in obs['warnings']:
        warn_d = dict(
            sc_id=sc_id,
            obsid=obs['target']['obsid'],
            obs_idx=obs_idx,
            warning=warn)
        idx_m = re.search('\[\s*(\d+)\]',
                          warn)
        if idx_m:
            warn_d['idx'] = idx_m.group(0)
        db.insert(warn_d, 'starcheck_warnings')
    db.commit()


config = dict(data_root='/data/aca/archive/starcheck',
              touch_file="starcheck_parser.touch",
              dbi='sqlite',
              server='starcheck.db3',
              mp_top_level='/data/mpcrit1/mplogs')


def get_options():
    parser = argparse.ArgumentParser(
        description="Update starcheck database")
    defaults = dict(config)
    parser.set_defaults(**defaults)
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--dbi",
                        help="dbi for starcheck database")
    parser.add_argument("--server",
                        help="server or file for database")
    parser.add_argument("--mp-top-level",
                        help="top level SOT MP dir")
    parser.add_argument("--touch-file")
    opt = parser.parse_args()
    return opt


def update(data_root, dbi, server, mp_top_level, touch_file):

    if dbi == 'sqlite':
        server = os.path.join(data_root, server)
    touch_file = os.path.join(data_root, touch_file)

    starchecks = Ska.Shell.bash(
        "find %s" % mp_top_level
        + "/[12]???/[A-Z][A-Z][A-Z][0-9][0-9][0-9][0-9]/ofls?/ "
        + " -maxdepth 1 -wholename '*/ofls?/starcheck.txt' "
        + "-cnewer %s" % touch_file)

    db = Ska.DBI.DBI(dbi=dbi, server=server, autocommit=False)
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
        for (obs, obs_idx) in izip(starcheck['obsdata'], count(0)):
            ingest_obs(obs, obs_idx, sc_id, st, db, existing=existing)
        logger.info("Done with %s; updating touch file" % st)
        Ska.Shell.bash("touch -r %s %s" % (st, touch_file))


def main():
    """
    Command line interface to update starcheck database from files
    in MP directories
    """
    opt = get_options()
    kwargs = vars(opt)
    update(**kwargs)

if __name__ == '__main__':
    main()
