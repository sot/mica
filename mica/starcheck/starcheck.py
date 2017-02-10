import os
import re
import logging
import argparse
from itertools import count
from operator import itemgetter
import numpy as np

from six.moves import zip
from Chandra.Time import DateTime
import Ska.Shell
import Ska.DBI
from kadi import events
from astropy.table import Table

from mica.starcheck.starcheck_parser import read_starcheck
from mica.common import MICA_ARCHIVE

logger = logging.getLogger('starcheck ingest')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

DEFAULT_CONFIG = dict(dbi='sqlite',
                      server=os.path.join(MICA_ARCHIVE, 'starcheck', 'starcheck.db3'),
                      mp_top_level='/data/mpcrit1/mplogs')
FILES = dict(data_root=os.path.join(MICA_ARCHIVE, 'starcheck'),
             touch_file=os.path.join(MICA_ARCHIVE, 'starcheck', "starcheck_parser.touch"),
             sql_def='starcheck.sql')

# this doesn't belong in starcheck
def get_timeline_at_date(date, timelines_db=None):
    date = DateTime(date).date
    if timelines_db is None:
        timelines_db = Ska.DBI.DBI(dbi='sqlite',
                         server='/proj/sot/ska/data/cmd_states/cmd_states.db3')
    return timelines_db.fetchone(
        "select * from timeline_loads where datestop >= '%s' "
        " and datestart <= '%s' and scs <= 130 order by datestart desc"
        % (date, date))


def get_catalog_at_date(date, config=None, timelines_db=None):
    """
    For a given date, return a dictionary describing the starcheck catalog that should apply.
    The content of that dictionary is from the database tables that parsed the starcheck report.
    A catalog is defined as applying, in this function, to any time from the end of the
    previous dwell through the end of the dwell in which the catalog was used.
    """

    if config is None:
        config = DEFAULT_CONFIG
    date = DateTime(date).date
    if timelines_db is None:
        timelines_db = Ska.DBI.DBI(dbi='sqlite',
                                   server='/proj/sot/ska/data/cmd_states/cmd_states.db3')
    last_tl = timelines_db.fetchone(
        "select max(datestop) as datestop from timelines")['datestop']
    first_tl = timelines_db.fetchone(
        "select min(datestart) as datestart from timelines")['datestart']
    db = Ska.DBI.DBI(dbi='sqlite', server=config['server'])
    # if we're outside of timelines, just try from the starcheck database
    if date > last_tl or date < first_tl:
        # Get one entry that is the last one before the specified time, in the most
        # recently ingested products directory
        starcheck = db.fetchone(
            """select * from starcheck_obs, starcheck_id
               where mp_starcat_time <= '{}' and mp_starcat_time > '{}'
               and starcheck_id.id = starcheck_obs.sc_id
               order by sc_id desc, mp_starcat_time desc """.format(
                date, (DateTime(date) - 1).date))
        if starcheck:
            cat_info = get_starcheck_catalog(starcheck['obsid'], mp_dir=starcheck['dir'])
            if date < first_tl:
                cat_info['status'] = 'ran before timelines'
            if date > last_tl:
                cat_info['status'] = 'planned'
            return cat_info
    # If timelines apply, try to use those to find a commanded catalog
    # Use kadi to get the first dwell that *ends* after the given time
    dwells = events.dwells.filter(stop__gte=date, subset=slice(None, 1))
    dwell = dwells[0]
    # Check the timelines that cover the range between the last dwell and the dwell
    # for this catalog, that should cover NMM when the catalog is commanded
    timelines = timelines_db.fetchall(
        """select * from timeline_loads where scs < 131
           and datestop > '{}' and datestart < '{}' order by datestart""".format(
           dwell.get_previous().stop, dwell.start))
    for timeline in timelines[::-1]:
        starchecks = db.fetchall(
            """select * from starcheck_obs, starcheck_id
               where dir = '{}'
               and mp_starcat_time >= '{}'
               and mp_starcat_time <= '{}' and mp_starcat_time <= '{}'
               and starcheck_id.id = starcheck_obs.sc_id
               order by mp_starcat_time """.format(
                timeline['mp_dir'],
                timeline['datestart'],
                timeline['datestop'], dwell.start))
        # The last one should be the one before beginning of the dwell
        if len(starchecks):
            cat_info = get_starcheck_catalog(starchecks[-1]['obsid'], mp_dir=starchecks[-1]['dir'])
            if cat_info['obs']['mp_starcat_time'] > DateTime().date:
                cat_info['status'] = 'approved'
            else:
                cat_info['status'] = 'ran'
            return cat_info
    return None


def get_mp_dir(obsid, config=None, timelines_db=None):
    """
    Get the mission planning directory for an obsid and some status information.  If the obsid catalog was
    used more than once (multi-obi or rescheduled after being maneuvered-to but not science observed during
    an SCS107 vehicle interval), return the directory and details of the last one used on the spacecraft.
    """
    if config is None:
        config = DEFAULT_CONFIG
    if timelines_db is None:
        timelines_db = Ska.DBI.DBI(dbi='sqlite',
                                   server='/proj/sot/ska/data/cmd_states/cmd_states.db3')
    last_tl = timelines_db.fetchone(
        "select max(datestop) as datestop from timelines")['datestop']
    with Ska.DBI.DBI(dbi='sqlite', server=config['server']) as db:
        starchecks = db.fetchall(
            """select * from starcheck_obs, starcheck_id
               where obsid = %d and starcheck_id.id = starcheck_obs.sc_id
               order by sc_id """ % obsid)
        # If this predates timelines data just return the last record that matches
        # and hope for the best
        if len(starchecks) and (starchecks[-1]['mp_starcat_time'] < '2001:001:00:00:00.000'):
            return (starchecks[-1]['dir'], 'ran_pretimelines', starchecks[-1]['mp_starcat_time'])
        # And if there's nothing do, a quick check to see if it is in other database
        # to throw an error if it is just missing in mica
        if not len(starchecks):
            tl_obsids = timelines_db.fetchall("select * from tl_obsids where obsid = {}".format(
                    obsid))
            if len(tl_obsids):
                raise ValueError(
                    "Obsid {} in parsed loads in tl_obsids but not in mica starcheck".format(obsid))
            return (None, None, None)
        # Go through the entries backwards (which are in ingest/date order)
        for sc in starchecks[::-1]:
            sc_date = sc['mp_starcat_time']
            # if this is in a schedule not-yet-in-timelines
            # We'll have to hope that is the most useful entry if there are multiples
            if sc_date > last_tl:
                return (sc['dir'], 'planned', sc_date)
            tl = get_timeline_at_date(sc_date, timelines_db=timelines_db)
            if not len(tl):
                raise ValueError("No timeline covering date {}".format(sc_date))
            # If the approved products in timelines are from a different directory, no-go
            if tl['mp_dir'] != sc['dir']:
                continue
            if sc_date < DateTime().date:
                return (sc['dir'], 'ran', sc_date)
            else:
                return (sc['dir'], 'approved', sc_date)
        raise ValueError("get_mp_dir should find something or nothing")


def obsid(obsid, mp_dir=None, config=None):
    if mp_dir is None:
        mp_dir, status, mp_date = get_mp_dir(obsid)
    if mp_dir is None:
        raise LookupError("No starcheck catalog found for {}".format(obsid))
    if config is None:
        config = DEFAULT_CONFIG
    db = Ska.DBI.DBI(dbi='sqlite', server=config['server'])
    sc = db.fetchone(
        "select id from starcheck_id where dir = '%s'" % mp_dir)
    if sc is None:
        raise LookupError("Catalog for {} missing.  Should already be in mica.starcheck db".format(obsid))
    sc_id = sc['id']
    sc = {}
    for d in ['manvr', 'catalog', 'obs', 'warnings']:
        sc[d] = db.fetchall(
            "select * from starcheck_%s where obsid = %d and sc_id = %d"
            % (d, obsid, sc_id))
    return sc


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



def get_starcheck_catalog(obsid, mp_dir=None,
                          config=None, tstart=None):
    if config is None:
        config = DEFAULT_CONFIG
    sc_dbi = config['dbi']
    sc_server = config['server']
    dbh = Ska.DBI.DBI(dbi=sc_dbi, server=sc_server)
    # tstart will override mp_dir
    if tstart is not None:
        # if this should be in timelines
        if DateTime(tstart).date > '2002:007':
            timeline = get_timeline_at_date(DateTime(tstart).date)
            if mp_dir is not None:
                if timeline['mp_dir'] != mp_dir:
                    raise ValueError(
                        "mp_dir was specified but inconsistent with tstart")
            mp_dir = timeline['mp_dir']
    if mp_dir is None:
        # in this mode, just get last NPNT interval
        # not clear about what obi we'll get, but..
        aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase',
                             user='aca_read')
        cmd_state = aca_db.fetchone(
            "select max(datestart) as datestart from cmd_states "
            "where obsid = %d and pcad_mode = 'NPNT'"
            % obsid)
        if cmd_state['datestart'] is None:
            # just the most recent starcheck if not in cmd_states
            obs = dbh.fetchall("""select * from starcheck_obs
                                  where obsid = %d""" % obsid)
            sc_id = np.max(obs['sc_id'])
            mp_dir = dbh.fetchone("select dir from starcheck_id "
                                  "where id = %d" % sc_id)['dir']
        else:
            timeline = aca_db.fetchone(
                "select * from timeline_loads where datestart in "
                "(select max(datestart) from timeline_loads "
                "where datestart <= '%s')" % cmd_state['datestart'])
            mp_dir = timeline['mp_dir']
    sc_id = dbh.fetchone("select id from starcheck_id "
                         "where dir = '%s'" % mp_dir)['id']
    obs = dbh.fetchone("select * from starcheck_obs "
                       "where sc_id = %d and obsid = %d" % (sc_id, obsid))
    cat = dbh.fetchall("""select * from starcheck_catalog
                      where sc_id = %d and obsid = %d
                      order by idx""" % (sc_id, obsid))
    manvr = dbh.fetchall("select * from starcheck_manvr "
                         "where sc_id = %d and obsid = %d"
                         % (sc_id, obsid))
    temp = dbh.fetchone("select pred_ccd_temp from starcheck_pred_temp "
                        "where sc_id = {} and obsid = {}".format(
            sc_id, obsid))
    pred_temp = None
    if temp:
        pred_temp = temp['pred_ccd_temp']
    return {'obs': obs,
            'cat': cat,
            'manvr': manvr,
            'mp_dir': mp_dir,
            'pred_temp': pred_temp}


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
    if config['dbi'] != 'sqlite':
        raise ValueError("Only sqlite DBI implemented")
    if (not os.path.exists(config['server'])
            or os.stat(config['server']).st_size == 0):
        if not os.path.exists(os.path.dirname(config['server'])):
            os.makedirs(os.path.dirname(config['server']))
        db_sql = os.path.join(os.environ['SKA_DATA'], 'mica', FILES['sql_def'])
        db_init_cmds = open(db_sql).read()
        db = Ska.DBI.DBI(dbi='sqlite', server=config['server'])
        db.execute(db_init_cmds)
        # make a touch file with a time before starchecks
        Ska.Shell.bash("touch -t 199801010000 %s" % (FILES['touch_file']),
                       env={'MAILCHECK': -1})
    else:
        db = Ska.DBI.DBI(dbi='sqlite', server=config['server'])

    starchecks = Ska.Shell.bash(
        "find %s" % config['mp_top_level']
        + "/[12]???/[A-Z][A-Z][A-Z][0-9][0-9][0-9][0-9]/ofls?/ "
        + " -maxdepth 1 -wholename '*/ofls?/starcheck.txt' "
        + "-cnewer %s" % FILES['touch_file'],
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
