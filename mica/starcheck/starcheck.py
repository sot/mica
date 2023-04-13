# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import re

import numpy as np
import Ska.DBI
import astropy.units as u
from astropy.table import Table

from cxotime import CxoTime
from kadi.commands import get_observations

from mica.common import MICA_ARCHIVE
from mica.utils import load_name_to_mp_dir

DEFAULT_CONFIG = dict(
    starcheck_db=dict(
        dbi="sqlite", server=os.path.join(MICA_ARCHIVE, "starcheck", "starcheck.db3")
    ),
    mp_top_level="/data/mpcrit1/mplogs",
)
FILES = dict(
    data_root=os.path.join(MICA_ARCHIVE, "starcheck"),
    touch_file=os.path.join(MICA_ARCHIVE, "starcheck", "starcheck_parser.touch"),
    sql_def="starcheck.sql",
)


OBS_CACHE = {}


def get_starcat(obsid, mp_dir=None):
    """
    Get the planned star catalog for an obsid.  Uses mica.starcheck database.

    :param obsid: obsid
    :param mp_dir: optional load specifier like 'FEB1317A' or '/2017/FEB1317/oflsa/'
    :returns: astropy.Table of star catalog
    """
    try:
        return get_starcheck_catalog(obsid, mp_dir)["cat"]
    except Exception:
        return None


def get_dither(obsid, mp_dir=None):
    """
    Get the planned dither for an obsid.  Uses mica.starcheck database.  Note that this
    does not have dither values for ERs and for some early mission observations.

    Observations with dither dither disabled will have yaw and pitch amplitude set
    to 0 and period of -999.

    :param obsid: obsid
    :param mp_dir: optional load specifier like 'FEB1317A' or '/2017/FEB1317/oflsa/'
    :returns: dictionary of planned dither parameters: 'yaw_ampl', 'pitch_ampl',
    'yaw_period', and 'pitch_period'.  amplitudes in arcseconds, periods in seconds.
    """
    obs = get_starcheck_catalog(obsid, mp_dir)["obs"]
    if obs["dither_state"] is None:
        raise ValueError(
            "Dither parameters not in starcheck database for obsid {}".format(obsid)
        )
    # For observations with no dither, period may be None, so just set to 1.0 to
    # avoid divide-by-None and divide-by-0 errors.  It is moot for 0 amplitude anyway.
    return {
        "yaw_ampl": obs["dither_y_amp"],
        "pitch_ampl": obs["dither_z_amp"],
        "yaw_period": obs["dither_y_period"] if obs["dither_y_amp"] != 0 else -999.0,
        "pitch_period": obs["dither_z_period"] if obs["dither_z_amp"] != 0 else -999.0,
    }


def get_att(obsid, mp_dir=None):
    """
    Get the planned ACA attitude from the starcheck database.

    :param obsid: obsid
    :param mp_dir: optional load specifier like 'FEB1317A' or '/2017/FEB1317/oflsa/'
    :returns: list of RA, Dec, Roll
    """
    obs = get_starcheck_catalog(obsid, mp_dir)["obs"]
    return [obs["point_ra"], obs["point_dec"], obs["point_roll"]]


def get_monitor_windows(start=None, stop=None, min_obsid=40000, config=None):
    """
    Use the database of starcheck products to get a list of monitor windows.

    This includes only catalogs that ran or will run.

    NOTE: by default only monitor windows in ER's are included, set ``min_obsid=0``
    to change this.

    :param start: Optional start time for filtering windows as fetched from the
        database
    :param stop: Optional stop time for filtering windows as fetched from the
        database
    :param min_obsid: Minimum obsid value for filtering.  Default of 40000 is
        intended to fetch only ERs

    :param config: config dictionary. If supplied must include 'starcheck_db'
                   key with a dictionary of the required arguments to Ska.DBI to connect
                   to that database.

    :returns: astropy Table of monitor windows.  See
              get_starcheck_catalog_at_date for description of the values of the
              'status' column.  The 'catalog' column contains the
              get_starcheck_catalog_at_date returned dictionary.
    """
    if config is None:
        config = DEFAULT_CONFIG

    start_date = CxoTime(start or "1999:001:12:00:00").date
    stop_date = CxoTime(stop).date
    with Ska.DBI.DBI(**config["starcheck_db"]) as db:
        mons = db.fetchall(
            f"""select obsid, mp_starcat_time as mp_starcat_date, type, sz, yang, zang, dir
                    from starcheck_catalog, starcheck_id
                    where type = 'MON' and obsid > {min_obsid}
                    and mp_starcat_date >= '{start_date}'
                    and mp_starcat_date <= '{stop_date}'
                    and starcheck_id.id = starcheck_catalog.sc_id"""  # noqa: E501
        )
    mons = Table(mons)
    mons["catalog"] = None

    obss = get_observations(start=start_date, stop=stop_date)
    obss_keys = set()
    for obs in obss:
        if obs["source"] != "CMD_EVT":
            key = (
                obs["obsid"],
                obs.get("starcat_date"),
                load_name_to_mp_dir(obs["source"]),
            )
            obss_keys.add(key)

    statuses = []
    date_now = CxoTime.now().date
    for mon in mons:
        sc_date = mon["mp_starcat_date"]
        key = (mon["obsid"], sc_date, mon["dir"])
        if key not in obss_keys:
            statuses.append("none")
        else:
            statuses.append("run" if sc_date < date_now else "approved")

    mons = mons[np.array(statuses) != "none"]

    return mons


def get_starcheck_catalog_at_date(date, starcheck_db=None):
    """
    For a given date, return a dictionary describing the starcheck catalog that should
    apply. The content of that dictionary is from the database tables that parsed the
    starcheck report. A catalog is defined as applying, in this function, to any time
    from the end of the previous dwell through the end of the dwell in which the catalog
    was used.

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
    - approved: observation in an approved future schedule

    :param date: Chandra.Time compatible date
    :param starcheck_db: optional handle to already-open starcheck database
    :returns: dictionary with starcheck content described above


    """
    date = CxoTime(date)
    if starcheck_db is None:
        starcheck_db = Ska.DBI.DBI(**DEFAULT_CONFIG["starcheck_db"])

    # Get observations within +/- 7 days and find the one where ``date`` is between end
    # of the previous obs dwell through the end of the current obs dwell.
    obss = Table(get_observations(start=date - 7 * u.day, stop=date + 7 * u.day))
    obss.sort("obs_start")
    for obs_prev, obs in zip(obss[:-1], obss[1:]):
        if obs_prev["obs_stop"] < date.date <= obs["obs_stop"]:
            break
    else:
        raise ValueError(f"could not find observation at {date}")

    mp_dir = load_name_to_mp_dir(obs["source"])

    # Kadi observations are labeled by the commanded obsid.  Starcheck catalogs
    # are labeled by the planned obsid used in the starcheck output.  For observations
    # after SCS107 with SOSA, the kadi observation obsids and starcheck obsids don't match.
    # Use a helper function to line these up by starcat time.
    obsid = get_starcheck_db_obsid(obs["starcat_date"], mp_dir)

    # Use the obsid and the known products directory to use the more generic
    # get_starcheck_catalog to fetch the right one from the database
    cat_info = get_starcheck_catalog(obsid, mp_dir=mp_dir, starcheck_db=starcheck_db)
    cat_info["status"] = "ran" if date < CxoTime().date else "approved"
    return cat_info


def get_mp_dir_from_starcheck_db(obsid, starcheck_db=None):
    """Get a guess for MP_DIR using the starchecks database.

    :param obsid: obsid
    :param starcheck_db: optional handle to already-open starcheck database
    :returns: directory, status, date
    """
    if starcheck_db is None:
        starcheck_db = Ska.DBI.DBI(**DEFAULT_CONFIG["starcheck_db"])

    starchecks = starcheck_db.fetchall(
        """select * from starcheck_obs, starcheck_id
            where obsid = %d and starcheck_id.id = starcheck_obs.sc_id
            order by sc_id """
        % obsid
    )
    # If this has a star catalog that predates command observations then return
    # the last record that matches and hope for the best
    obss = get_observations(start="1999:001", stop="2003:001", scenario="flight")
    if (
        len(starchecks)
        and (starchecks[-1]["mp_starcat_time"] is not None)
        and (starchecks[-1]["mp_starcat_time"] < obss[0]["obs_start"])
    ):
        out = (
            starchecks[-1]["dir"],
            "ran_pre_commands",
            starchecks[-1]["mp_starcat_time"],
        )
    else:
        out = (None, None, None)

    return out


def get_mp_dir(obsid, starcheck_db=None):
    """
    Get the mission planning directory for an obsid and some status information.
    If the obsid catalog was used more than once (multi-obi or rescheduled after
    being used in a vehicle-only interval), return the directory and details of
    the last one used on the spacecraft.

    Only observations from approved schedules are considered.

    The returned ``directory`` is a string like "/2006/DEC2506/oflsc/" that
    describes the directory that was used for the products with this star
    catalog.

    The returned ``status`` has possible values:

    - "ran": observation date before current time
    - "approved": observation date after current time
    - "no starcat": observation exists but has no star catalog

    The return ``date`` is the date/time of the ``MP_STARCAT`` time.

    :param obsid: obsid
    :param starcheck_db: optional handle to already-open starcheck database
    :returns: directory, status, date

    """
    try:
        obs = get_observations(obsid=obsid)[-1]
    except ValueError as err:
        if "No matching observations" in str(err):
            mp_dir, status, sc_date = get_mp_dir_from_starcheck_db(obsid, starcheck_db)
        else:
            raise err
    else:
        sc_date = obs.get("starcat_date")
        if sc_date is None:
            status = "no starcat"
        else:
            if sc_date < CxoTime.now().date:
                status = "ran"
            else:
                status = "approved"
        mp_dir = load_name_to_mp_dir(obs["source"])

    return mp_dir, status, sc_date


def get_starcheck_db_obsid(mp_starcat_time, mp_dir, starcheck_db=None):
    """
    Get starcheck database obsid for a given MP_STARCAT time and products directory.

    :param mp_starcat_time: MP_STARCAT time
    :param mp_dir: products directory (e.g. /2006/DEC2506/oflsc/)
    :param starcheck_db: optional handle to already-open starcheck database

    :returns: obsid
    """
    if starcheck_db is None:
        starcheck_db = Ska.DBI.DBI(**DEFAULT_CONFIG["starcheck_db"])

    result = starcheck_db.fetchone(
        "SELECT starcheck_obs.obsid FROM starcheck_obs "
        "INNER JOIN starcheck_id "
        "ON starcheck_id.id = starcheck_obs.sc_id "
        f"AND starcheck_id.dir = '{mp_dir}'"
        f"AND starcheck_obs.mp_starcat_time = '{mp_starcat_time}'"
    )
    if result is None:
        raise ValueError(f"no starcheck entry for {mp_starcat_time=} and {mp_dir=}")

    obsid = result["obsid"]
    return obsid


def get_starcheck_catalog(obsid, mp_dir=None, starcheck_db=None):
    """
    For a given obsid, return a dictionary describing the starcheck catalog that should
    apply. The content of that dictionary is from the database tables of that parsed the
    starcheck report and has keys:

    - cat: catalog rows as astropy.table
    - manvr: list of maneuvers to this attitude
    - pred_temp: predicted ACA CCD temperature
    - warnings: list of warnings below catalog in starcheck output
    - obs: dictionary of observation target and pointing information
    - mp_dir: directory with products that are the source of this catalog data
    - status: string describing status of that observation, described below.

    Status:

    - ran: observation approved and has date before current time
    - approved: observation approved and has date after current time
    - ran_pre_commands: ran, but before commands archive starts

    :param obsid: obsid
    :param mp_dir: optional load specifier like 'FEB1317A' or '/2017/FEB1317/oflsa/'. By
                   default the as-run loads (via ``get_mp_dir()``) are used.
    :param starcheck_db: optional handle to already-open starcheck database
    :returns: dictionary with starcheck content described above
    """
    # Keep the original kwarg from user for the caching key later
    mp_dir_kwarg = mp_dir

    if (obsid, mp_dir) in OBS_CACHE:
        return OBS_CACHE[obsid, mp_dir]

    if starcheck_db is None:
        starcheck_db = Ska.DBI.DBI(**DEFAULT_CONFIG["starcheck_db"])
    status = None
    if mp_dir is None:
        mp_dir, status, obs_date = get_mp_dir(obsid)

    # if it is still none, there's nothing to try here
    if mp_dir is None:
        return None

    # mp_dir is in the standard short-form "MAY3018A".  Translate to
    # '/2018/MAY3018/oflsa/'
    if re.match(r"[A-Z]{3}\d{4}[A-Z]$", mp_dir):
        mp_dir = load_name_to_mp_dir(mp_dir)

    db = starcheck_db  # shorthand for the rest of the routine
    sc_id = db.fetchone("select id from starcheck_id where dir = '%s'" % mp_dir)
    if sc_id is None:
        return None
    sc_id = sc_id["id"]
    sc = {"mp_dir": mp_dir, "status": status}
    sc["obs"] = db.fetchone(
        "select * from starcheck_obs where obsid = {} and sc_id = {}".format(
            obsid, sc_id
        )
    )
    for d in ["manvr", "catalog", "warnings"]:
        table_entries = db.fetchall(
            "select * from starcheck_%s where obsid = %d and sc_id = %d"
            % (d, obsid, sc_id)
        )
        sc[d] = Table(table_entries) if len(table_entries) else []
    sc["cat"] = sc["catalog"]
    del sc["catalog"]
    pred_temp = db.fetchone(
        "select pred_ccd_temp from starcheck_pred_temp "
        "where sc_id = {} and obsid = {}".format(sc_id, obsid)
    )
    if pred_temp is not None and "pred_ccd_temp" in pred_temp:
        sc["pred_temp"] = pred_temp["pred_ccd_temp"]

    # Cache result for later re-use (use whatever user supplied as mp_dir kwarg)
    # since mp_dir gets munged in the meantime.
    OBS_CACHE[obsid, mp_dir_kwarg] = sc

    return sc
