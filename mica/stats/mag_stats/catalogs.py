import os
import numpy as np

from astropy.table import Table
from astropy import table
from Ska.DBI import DBI

from kadi import events
from cxotime import CxoTime
import tables
from kadi.commands import get_cmds
from Chandra.Time import date2secs

TIMELINES = None
STARCAT_CMDS = None
DWELLS_NP = None
DWELLS_MAP = None
STARS_OBS = None
STARS_OBS_INDICES = None
STARS_OBS_MAP = None
STARS_OBS_NP = None


def _load_startcat_commands():
    """
    Get all star catalog commands.  This is a definitive list of commands
    *actually run* on the spacecraft (to the best of our knowledge).

    Note this one case in the cmds archive of duplicate commands.  Not clear
    if this is real or an artifact of bad timelines.

    <Table length=7>
     idx            date            type     tlmsid   scs   step  timeline_id   vcdu
    uint16        bytes21         bytes12   bytes10  uint8 uint16    uint32    int32
    ------ --------------------- ---------- -------- ----- ------ ----------- --------
        4 2006:295:18:59:08.748 MP_STARCAT AOSTRCAT   128     25   277930627   181961
        4 2006:295:18:59:08.748 MP_STARCAT AOSTRCAT   130    706   277867085   181961
        6 2006:295:18:59:12.999 COMMAND_SW AOMANUVR   130    757   277867085   181977
        6 2006:295:18:59:12.999 COMMAND_SW AOMANUVR   128     76   277930627   181977

    :return:
    """
    global STARCAT_CMDS
    global DWELLS_NP
    global DWELLS_MAP

    STARCAT_CMDS = Table(get_cmds('2003:001', type='MP_STARCAT'))['date', 'timeline_id']
    STARCAT_CMDS.rename_column('date', 'mp_starcat_time')
    # print(len(starcat_cmds))
    STARCAT_CMDS = table.unique(STARCAT_CMDS, keys='mp_starcat_time')
    # print(len(starcat_cmds))

    # Get all maneuvers and use this to create a table of NPM dwells
    # with a star catalog.
    #
    # Note that obsid is not unique for maneuvers since
    # after a safing action the obsid in telemetry can stay fixed while maneuvers
    # continue happening.
    manvrs = events.manvrs.filter('2003:001')
    rows = [(manvr.obsid, manvr.start, manvr.kalman_start, manvr.npnt_stop) for manvr in manvrs]

    # Form dwells table. At this point it includes maneuvers that might not have
    # a star catalog with NPM dwell.  That gets fixed in the next bit.
    dwells = Table(rows=rows, names=('obsid', 'manvr_start', 'start', 'stop'))

    # Clip starcat_cmds to be within time range of dwells
    ok = (STARCAT_CMDS['mp_starcat_time'] < dwells['manvr_start'][-1])
    STARCAT_CMDS = STARCAT_CMDS[ok]

    # Now make a time-based association between every star catalog command and the
    # subsequent maneuver command.  By the way that commands are built it is req'd
    # that a star catalog command be followed by a tlmsid=MANUVR command which
    # corresponds to `manvr_start` in the kadi manvrs events.  Empirically this is
    # always less than 1000 sec delay (so allow up to 1200).
    idxs = np.searchsorted(dwells['manvr_start'], STARCAT_CMDS['mp_starcat_time'])

    # Check that the association is correct.
    dt = CxoTime(dwells['manvr_start'][idxs]) - CxoTime(STARCAT_CMDS['mp_starcat_time'])
    ok = (dt.sec > 1) & (dt.sec < 1200)

    print(f'WARNING: {np.sum(~ok)} entries where star-catalog/maneuver association is not right')
    #assert np.all((dt.sec > 1) & (dt.sec < 1200))

    STARCAT_CMDS = STARCAT_CMDS[ok]

    # Now re-select dwells to be only dwells with a star catalog and set the
    # corresponding command time for MP_STARCAT.  This is the unique key to
    # allow going from a star catalog to the NPM dwell that uses the catalog.
    dwells = dwells[idxs][ok]
    dwells['mp_starcat_time'] = STARCAT_CMDS['mp_starcat_time']
    dwells.add_index('mp_starcat_time')

    # Kick out a few cases with None values for start or stop
    ok = [row['start'] is not None and row['stop'] is not None for row in dwells]
    dwells = dwells[ok]

    # Now make a numpy array version for speed and a map (fast version of .loc)
    dwells['tstart'] = date2secs(dwells['start'].tolist())
    dwells['tstop'] = date2secs(dwells['stop'].tolist())
    DWELLS_NP = dwells['mp_starcat_time', 'tstart', 'tstop'].as_array()
    DWELLS_MAP = {DWELLS_NP['mp_starcat_time'][idx].decode('ascii'): idx
                  for idx in range(len(DWELLS_NP))}


def _load_observed_catalogs():
    """
    Get actually observed star catalog entries for the mission.

    - `timelines` has mission planning `dir` and only includes loads that were approved
       and at least partially run.

    - `starcheck_id` table maps `sc_id` (called `id` there) to mission planning `dir`

    - `starcheck_catalog` has (`obsid`, `slot`, `id`, `sc_id`, `type`, `mag` and more) for every
       catalog ever created (including those never run) so we can map from `sc_id` to `dir` to see
       if a star catalog entry *might* have been run.  But not all approved load dirs were
       completed so we need to also cross-correlate with obsids from telemetry.
    """
    global TIMELINES
    global STARS_OBS
    global STARS_OBS_INDICES
    global STARS_OBS_MAP
    global STARS_OBS_NP

    # Timelines defines the mission planning `dir` of loads that were approved and
    # at least partially run.

    cmd_states_db3 = os.path.join(os.environ['SKA'], 'data', 'cmd_states', 'cmd_states.db3')
    with DBI(server=cmd_states_db3, dbi='sqlite') as timelines_db:
        TIMELINES = Table(timelines_db.fetchall('select dir, id from timelines '
                                                'where datestart > "2003:007"'))
    TIMELINES.rename_column('id', 'timeline_id')
    TIMELINES.add_index('dir')
    TIMELINES.add_index('timeline_id')

    # Get table of every star catalog entry (a single star or fid light) that
    # was every planned.

    # The output `cat_entries` table is not useful by itself because it has
    # multiple repeats for each entry.  The input starcheck tables only provide
    # the source MP directory for a catalog, not the timeline ID.  Thus the output
    # has an entry for each possible timeline_id that might correspond to the
    # catalog.  In the next step of joining with definitive MP_STARCAT commands
    # this all gets resolved back down.

    server = os.path.join(os.environ['SKA'],
                          'data', 'mica', 'archive', 'starcheck', 'starcheck.db3')
    with DBI(server=server, dbi='sqlite') as db:
        cat_entries = db.fetchall(
            'select obsid, slot, starcheck_catalog.id, sc_id, type, mag, mp_starcat_time, dir '
            'from starcheck_catalog, starcheck_id '
            'where starcheck_catalog.id is not NULL '
            'AND starcheck_catalog.sc_id = starcheck_id.id')

    # For this application filter down to guide stars right away
    ok = np.in1d(cat_entries['type'], ['GUI', 'BOT'])
    cat_entries = Table(cat_entries[ok])
    cat_entries.rename_column('id', 'agasc_id')
    cat_entries = table.join(cat_entries, TIMELINES, keys='dir')
    # print(len(cat_entries))
    # cat_entries[-5:].pprint(max_width=-1)

    # Now only accept catalog entries that were actually run by
    # joining `cat_entries` (which has a lot of candidates that might have been run)
    # with `starcat_cmds` that were actually run in a timeline.

    ces = table.join(cat_entries, STARCAT_CMDS, keys=['mp_starcat_time', 'timeline_id'])
    # print(len(ces))
    # print(ces[-20:])
    # ces[ces['obsid'] == 19809]

    # Add mag_aca_err column

    filename = os.path.join(os.environ['SKA'], 'data', 'agasc', 'proseco_agasc_1p7.h5')
    with tables.open_file(filename) as h5:
        agasc_ids = h5.root.data.col('AGASC_ID')
        mag_errs = h5.root.data.col('MAG_ACA_ERR') * 0.01

    tt = Table([agasc_ids, mag_errs], names=['agasc_id', 'mag_aca_err'])
    ces = table.join(ces, tt, keys='agasc_id')
    # print(ces[-5:])

    # Make final table of observed guide star catalog entries `star_obs`
    # ------------------------------------------------------------------

    # Group table by agasc_id

    STARS_OBS = ces.group_by('agasc_id')
    STARS_OBS.add_index('agasc_id')
    STARS_OBS.add_index('mp_starcat_time')

    # Make a table that contains the indices into `star_obs` and gives the
    # repeat count `n_obs` for each star.

    indices = STARS_OBS.groups.indices
    rows = []
    for idx0, idx1 in zip(indices[:-1], indices[1:]):
        agasc_id = STARS_OBS['agasc_id'][idx0]
        rows.append((agasc_id, idx1 - idx0, idx0, idx1))
    STARS_OBS_INDICES = Table(rows=rows, names=['agasc_id', 'n_obs', 'idx0', 'idx1'])
    STARS_OBS_INDICES.sort('n_obs')
    STARS_OBS_INDICES.add_index('agasc_id')
    # print(stars_obs_indices[:5])
    # print(stars_obs_indices[-10:])

    STARS_OBS_MAP = {row['agasc_id']: (row['idx0'], row['idx1']) for row in STARS_OBS_INDICES}
    STARS_OBS_NP = STARS_OBS.as_array()


_load_startcat_commands()
_load_observed_catalogs()
