# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import re
import numpy as np

from Chandra.Time import DateTime
from Chandra.cmd_states import get_cmd_states
import Ska.DBI
from kadi import events
from astropy.table import Table

from mica.utils import get_timeline_at_date
from mica.common import MICA_ARCHIVE

DEFAULT_CONFIG = dict(starcheck_db=dict(dbi='sqlite',
                                        server=os.path.join(MICA_ARCHIVE, 'starcheck', 'starcheck.db3')),
                      mp_top_level='/data/mpcrit1/mplogs',
                      timelines_db=dict(dbi='sqlite',
                                        server=os.path.join(os.environ['SKA'], 'data', 'cmd_states', 'cmd_states.db3')))
FILES = dict(data_root=os.path.join(MICA_ARCHIVE, 'starcheck'),
             touch_file=os.path.join(MICA_ARCHIVE, 'starcheck', "starcheck_parser.touch"),
             sql_def='starcheck.sql')


OBS_CACHE = {}


def get_starcat(obsid):
    """
    Get the planned star catalog for an obsid.  Uses mica.starcheck database.

    :param obsid: obsid
    :returns: astropy.Table of star catalog
    """
    if obsid not in OBS_CACHE:
        OBS_CACHE[obsid] = get_starcheck_catalog(obsid)
    return OBS_CACHE[obsid]['cat']


def get_dither(obsid):
    """
    Get the planned dither for an obsid.  Uses mica.starcheck database.  Note that this
    does not have dither values for ERs and for some early mission observations.

    Observations with dither dither disabled will have yaw and pitch amplitude set
    to 0 and period of -999.

    :param obsid: obsid
    :returns: dictionary of planned dither parameters: 'yaw_ampl', 'pitch_ampl',
    'yaw_period', and 'pitch_period'.  amplitudes in arcseconds, periods in seconds.
    """
    if obsid not in OBS_CACHE:
        OBS_CACHE[obsid] = get_starcheck_catalog(obsid)
    obs = OBS_CACHE[obsid]['obs']
    if obs['dither_state'] is None:
        raise ValueError("Dither parameters not in starcheck database for obsid {}".format(obsid))
    # For observations with no dither, period may be None, so just set to 1.0 to
    # avoid divide-by-None and divide-by-0 errors.  It is moot for 0 amplitude anyway.
    return {'yaw_ampl': obs['dither_y_amp'],
            'pitch_ampl': obs['dither_z_amp'],
            'yaw_period': obs['dither_y_period'] if obs['dither_y_amp'] != 0 else -999.0,
            'pitch_period': obs['dither_z_period'] if obs['dither_z_amp'] != 0 else -999.0}


def get_att(obsid):
    """
    Get the planned ACA attitude from the starcheck database.

    :param obsid: obsid
    :returns: list of RA, Dec, Roll
    """
    if obsid not in OBS_CACHE:
        OBS_CACHE[obsid] = get_starcheck_catalog(obsid)
    obs = OBS_CACHE[obsid]['obs']
    return [obs['point_ra'],
            obs['point_dec'],
            obs['point_roll']]


def get_monitor_windows(start=None, stop=None, min_obsid=40000, config=None):
    """
    Use the database of starcheck products to get a list of monitor windows
    This list is filtered by timelines content to only include catalogs that should have or will run.

    :param start: Optional start time for filtering windows as fetched from the database
    :param stop: Optional stop time for filtering windows as fetched from the database
    :param min_obsid: Minimum obsid value for filtering.  Default of 40000 is intended to fetch only ERs
    :param config: config dictionary. If supplied must include 'starcheck_db' and 'timelines_db' keys
                   with dictionaries of the required arguments to Ska.DBI to connect to those databases.
    :returns: astropy Table of monitor windows.  See get_starcheck_catalog_at_date for description of the values
              of the 'status' column.  The 'catalog' column contains the get_starcheck_catalog_at_date
              returned dictionary.
    """
    if config is None:
        config = DEFAULT_CONFIG
    start_date = DateTime(start or '1999:001').date
    stop_date = DateTime(stop).date
    with Ska.DBI.DBI(**config['starcheck_db']) as db:
        mons = db.fetchall("""select obsid, mp_starcat_time as mp_starcat_date, type, sz, yang, zang, dir
                          from starcheck_catalog, starcheck_id
                          where type = 'MON' and obsid > {}
                          and mp_starcat_date >= '{}'
                          and mp_starcat_date <= '{}'
                          and starcheck_id.id = starcheck_catalog.sc_id""".format(
                min_obsid, start_date, stop_date))
    mons = Table(mons)
    mons['catalog'] = None
    with Ska.DBI.DBI(**config['timelines_db']) as timelines_db:
        statuses = []
        # Check which ones actually ran or are likely to run
        for mon in mons:
            try:
                catalog = get_starcheck_catalog_at_date(mon['mp_starcat_date'],
                                                        timelines_db=timelines_db)
                if (catalog is not None and catalog['obs']['obsid'] == mon['obsid']
                        and catalog['mp_dir'] == mon['dir']):
                    mon['catalog'] = catalog
                    statuses.append(catalog['status'])
                else:
                    statuses.append('none')
            except LookupError:
                statuses.append('none')
        # Add that status info to the table and filter on it
        mons['status'] = statuses
        mons = mons[mons['status'] != 'none']
    return mons


def get_starcheck_catalog_at_date(date, starcheck_db=None, timelines_db=None):
    """
    For a given date, return a dictionary describing the starcheck catalog that should apply.
    The content of that dictionary is from the database tables that parsed the starcheck report.
    A catalog is defined as applying, in this function, to any time from the end of the
    previous dwell through the end of the dwell in which the catalog was used.

    Star catalog dictionary with keys:

    - cat: catalog rows as astropy.table
    - manvr: list of maneuvers to this attitude
    - pred_temp: predicted ACA CCD temperature
    - warnings: list of warnings below catalog in starcheck output
    - obs: dictionary of observation target and pointing information
    - mp_dir: directory with products that are the source of this catalog data
    - status: string describing status of that observation, described below.

    Status:

    - ran: observation was observed
    - planned: observation in a not-approved future schedule
    - approved: observation in an approved future schedule (ingested in timelines/cmd_states)
    - ran_pretimelines: ran, but before timelines database starts
    - timelines_gap: after timelines database start but missing data
    - no starcat: in the database but has no star catalog

    :param date: Chandra.Time compatible date
    :param starcheck_db: optional handle to already-open starcheck database
    :param timelines_db: optional handle to already-open timelines database
    :returns: dictionary with starcheck content described above


    """
    date = DateTime(date).date
    if starcheck_db is None:
        starcheck_db = Ska.DBI.DBI(**DEFAULT_CONFIG['starcheck_db'])
    db = starcheck_db
    if timelines_db is None:
        timelines_db = Ska.DBI.DBI(**DEFAULT_CONFIG['timelines_db'])
    last_tl = timelines_db.fetchone(
        "select max(datestop) as datestop from timelines")['datestop']
    first_tl = timelines_db.fetchone(
        "select min(datestart) as datestart from timelines")['datestart']
    # Check kadi to get the first dwell that *ends* after the given time
    dwells = events.dwells.filter(stop__gte=date, subset=slice(None, 1))
    # if we're outside of timelines or not yet in kadi, just try from the starcheck database
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
                cat_info['status'] = 'ran_pretimelines'
            if date > last_tl:
                cat_info['status'] = 'planned'
            return cat_info

    # We want to search for legitimate commanding that would cover the time when a star
    # catalog would have been commanded for this dwell.  This is generally the time range
    # between the end of the previous dwell and the beginning of this dwell.  However, if
    # there are multiple dwells from one maneuver, use the beginning of NMM from that one
    # maneuver else, use the end of the last dwell.  Don't use nman_start time by default
    # because that doesn't appear to work if the catalog was commanded in a nonstandard
    # nmm sequence like dark cal.

    # There is a tiny window of time in cmd_states but not yet in kadi, but this code tries to
    # grab the dwell and maneuver that would be related to a date in that range
    if date < last_tl and len(dwells) == 0:
        pcad_states = get_cmd_states.fetch_states(start=DateTime(date) - 2, vals=['pcad_mode'])
        dwell = pcad_states[(pcad_states['pcad_mode'] == 'NPNT') & (pcad_states['datestop'] >= date)][0]
        manvr = pcad_states[(pcad_states['pcad_mode'] == 'NMAN')
                            & (pcad_states['datestop'] <= dwell['datestart'])][-1]
        start_cat_search = manvr['datestart']
        dwell_start = dwell['datestart']
    else:
        # If we have a dwell from kadi, use it to search for commanding
        dwell = dwells[0]
        dwell_start = dwell.start
        # Try to use the beginning of the previous nman period to define when the catalog
        # should have been commanded.  If there are multiple dwells for attitude, try to
        # use nman_start if available.
        if dwell.manvr.n_dwell > 1 and dwell.manvr.nman_start is not None:
            start_cat_search = dwell.manvr.nman_start
        else:
            start_cat_search = dwell.get_previous().stop

    timelines = timelines_db.fetchall(
            """select * from timeline_loads where scs < 131
           and datestop > '{}' and datestart < '{}' order by datestart""".format(
            start_cat_search, dwell_start))
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
                timeline['datestop'], dwell_start))
        # The last one should be the one before beginning of the dwell
        if len(starchecks):
            # Use the obsid and the known products directory to use the more generic get_starcheck_catalog
            # to fetch the right one from the database
            cat_info = get_starcheck_catalog(starchecks[-1]['obsid'], mp_dir=starchecks[-1]['dir'])
            cat_info['status'] = 'ran' if date < DateTime().date else 'approved'
            return cat_info
    return None


def get_mp_dir(obsid, starcheck_db=None, timelines_db=None):
    """
    Get the mission planning directory for an obsid and some status information.  If the obsid catalog was
    used more than once (multi-obi or rescheduled after being used in a vehicle-only interval), return the
    directory and details of the last one used on the spacecraft.

    The returned directory describes the directory that was used for the products with this star catalog.
    The returned status has possible values:

    - ran: observation was observed
    - planned: observation in a not-approved future schedule
    - approved: observation in an approved future schedule (ingested in timelines/cmd_states)
    - ran_pretimelines: ran, but before timelines database starts
    - timelines_gap: after timelines database start but missing data
    - no starcat: in the database but has no star catalog

    The return 'date' is the date/time of the `MP_STARCAT` time.

    :param obsid: obsid
    :param starcheck_db: optional handle to already-open starcheck database
    :param timelines_db: optional handle to already-open timelines database
    :returns: directory, status, date .   Described above.

    """
    if starcheck_db is None:
        starcheck_db = Ska.DBI.DBI(**DEFAULT_CONFIG['starcheck_db'])
    db = starcheck_db
    if timelines_db is None:
        timelines_db = Ska.DBI.DBI(**DEFAULT_CONFIG['timelines_db'])
    last_tl = timelines_db.fetchone(
        "select max(datestop) as datestop from timelines")['datestop']
    starchecks = db.fetchall(
        """select * from starcheck_obs, starcheck_id
           where obsid = %d and starcheck_id.id = starcheck_obs.sc_id
           order by sc_id """ % obsid)
    # If this has a star catalog that predates timelines data just return the last record that matches
    # and hope for the best
    if (len(starchecks) and (starchecks[-1]['mp_starcat_time'] is not None) and
        (starchecks[-1]['mp_starcat_time'] < '2001:001:00:00:00.000')):
        return (starchecks[-1]['dir'], 'ran_pretimelines', starchecks[-1]['mp_starcat_time'])
    # And if there are no entries, just return None for dir, status, time
    if not len(starchecks):
        return (None, None, None)
    # Go through the entries backwards (which are in ingest/date order)
    for sc in starchecks[::-1]:
        sc_date = sc['mp_starcat_time']
        # If it has no star catalog, there's not much to do
        if sc_date is None:
            return (sc['dir'], 'no starcat', sc_date)
        # If this is in a schedule not-yet-in-timelines or has no star catalog,
        # we'll have to hope that is the most useful entry if there are multiples
        if sc_date > last_tl:
            return (sc['dir'], 'planned', sc_date)
        tl = get_timeline_at_date(sc_date, timelines_db=timelines_db)
        # If there is a gap in timelines at this date, just return the most recent starcheck entry
        if tl is None or not len(tl):
            return (starchecks[-1]['dir'], 'timelines_gap', starchecks[-1]['mp_starcat_time'])
        # If the approved products in timelines are from a different directory, no-go
        if tl['mp_dir'] != sc['dir']:
            continue
        if sc_date < DateTime().date:
            return (sc['dir'], 'ran', sc_date)
        else:
            return (sc['dir'], 'approved', sc_date)
    raise ValueError("get_mp_dir should find something or nothing")


def get_starcheck_catalog(obsid, mp_dir=None, starcheck_db=None, timelines_db=None):
    """
    For a given obsid, return a dictionary describing the starcheck catalog that should apply.
    The content of that dictionary is from the database tables of that parsed the starcheck report
    and has keys:

    - cat: catalog rows as astropy.table
    - manvr: list of maneuvers to this attitude
    - pred_temp: predicted ACA CCD temperature
    - warnings: list of warnings below catalog in starcheck output
    - obs: dictionary of observation target and pointing information
    - mp_dir: directory with products that are the source of this catalog data
    - status: string describing status of that observation, described below.

    Status:

    - ran: observation was observed
    - planned: observation in a not-approved future schedule
    - approved: observation in an approved future schedule (ingested in timelines/cmd_states)
    - ran_pretimelines: ran, but before timelines database starts
    - timelines_gap: after timelines database start but missing data

    :param obsid: obsid
    :param mp_dir: mission planning directory (in the form '/2017/FEB1317/oflsa/') to which to limit
                   searches for the obsid.  If 'None', get_mp_dir() will be used to select appropriate directory.
    :param starcheck_db: optional handle to already-open starcheck database
    :param timelines_db: optional handle to already-open timelines database
    :returns: dictionary with starcheck content described above
    """
    if starcheck_db is None:
        starcheck_db = Ska.DBI.DBI(**DEFAULT_CONFIG['starcheck_db'])
    if timelines_db is None:
        timelines_db = Ska.DBI.DBI(**DEFAULT_CONFIG['timelines_db'])
    status = None
    if mp_dir is None:
        mp_dir, status, obs_date = get_mp_dir(obsid, starcheck_db=starcheck_db, timelines_db=timelines_db)
    # if it is still none, there's nothing to try here
    if mp_dir is None:
        return None
    db = starcheck_db # shorthand for the rest of the routine
    sc_id = db.fetchone("select id from starcheck_id "
                                  "where dir = '%s'" % mp_dir)
    if sc_id is None:
        return None
    sc_id = sc_id['id']
    sc = {'mp_dir': mp_dir,
          'status': status}
    sc['obs'] = db.fetchone("select * from starcheck_obs where obsid = {} and sc_id = {}".format(
            obsid, sc_id))
    for d in ['manvr', 'catalog', 'warnings']:
        table_entries = db.fetchall(
            "select * from starcheck_%s where obsid = %d and sc_id = %d"
            % (d, obsid, sc_id))
        sc[d] = Table(table_entries) if len(table_entries) else []
    sc['cat'] = sc['catalog']
    del sc['catalog']
    pred_temp = db.fetchone("select pred_ccd_temp from starcheck_pred_temp "
                            "where sc_id = {} and obsid = {}".format(
            sc_id, obsid))
    if pred_temp is not None and 'pred_ccd_temp' in pred_temp:
        sc['pred_temp'] = pred_temp['pred_ccd_temp']
    return sc

