# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import os
import time
import sys
import argparse
import numpy as np
import tables
from astropy.table import Table
import logging
import warnings

from Chandra.Time import DateTime
import agasc
from kadi import events
from Ska.engarchive import fetch, fetch_sci
import mica.archive.obspar
from mica.starcheck.starcheck import get_starcheck_catalog_at_date
import Ska.astro
from Quaternion import Quat
from chandra_aca import dark_model
from chandra_aca.transform import radec_to_eci

# Ignore known numexpr.necompiler and table.conditions warning
warnings.filterwarnings(
    'ignore',
    message="using `oa_ndim == 0` when `op_axes` is NULL is deprecated.*",
    category=DeprecationWarning)


logger = logging.getLogger('star_stats')
logger.setLevel(logging.INFO)
if not len(logger.handlers):
    logger.addHandler(logging.StreamHandler())

STAT_VERSION = 0.6

GUIDE_COLS = {
    'obs': [
        ('obsid', 'int'),
        ('obi', 'int'),
        ('kalman_tstart', 'float'),
        ('npnt_tstop', 'float'),
        ('kalman_datestart', 'S21'),
        ('npnt_datestop', 'S21'),
        ('revision', 'S15')],
    'cat': [
        ('slot', 'int'),
        ('idx', 'int'),
        ('type', 'S5'),
        ('yang', 'float'),
        ('zang', 'float'),
        ('sz', 'S4'),
        ('mag', 'float')],
    'stat': [
        ('n_samples', 'int'),
        ('n_track', 'int'),
        ('f_track', 'float'),
        ('f_racq', 'float'),
        ('f_srch', 'float'),
        ('f_none', 'float'),
        ('n_kalman', 'int'),
        ('no_track', 'float'),
        ('f_within_0.3', 'float'),
        ('f_within_1', 'float'),
        ('f_within_3', 'float'),
        ('f_within_5', 'float'),
        ('f_outside_5', 'float'),
        ('f_obc_bad', 'float'),
        ('f_common_col', 'float'),
        ('f_quad_bound', 'float'),
        ('f_sat_pix', 'float'),
        ('f_def_pix', 'float'),
        ('f_ion_rad', 'float'),
        ('f_mult_star', 'float'),
        ('aoacmag_min', 'float'),
        ('aoacmag_mean', 'float'),
        ('aoacmag_max', 'float'),
        ('aoacmag_5th', 'float'),
        ('aoacmag_16th', 'float'),
        ('aoacmag_50th', 'float'),
        ('aoacmag_84th', 'float'),
        ('aoacmag_95th', 'float'),
        ('aoacmag_std', 'float'),
        ('aoacyan_mean', 'float'),
        ('aoaczan_mean', 'float'),
        ('dy_min', 'float'),
        ('dy_mean', 'float'),
        ('dy_std', 'float'),
        ('dy_max', 'float'),
        ('dz_min', 'float'),
        ('dz_mean', 'float'),
        ('dz_std', 'float'),
        ('dz_max', 'float'),
        ('dr_min', 'float'),
        ('dr_mean', 'float'),
        ('dr_std', 'float'),
        ('dr_5th', 'float'),
        ('dr_95th', 'float'),
        ('dr_max', 'float'),
        ('n_track_interv', 'int'),
        ('n_long_track_interv', 'int'),
        ('n_long_no_track_interv', 'int'),
        ('n_racq_interv', 'int'),
        ('n_srch_interv', 'int'),
        ],
    'agasc': [
        ('agasc_id', 'int'),
        ('color', 'float'),
        ('ra', 'float'),
        ('dec', 'float'),
        ('epoch', 'float'),
        ('pm_ra', 'int'),
        ('pm_dec', 'int'),
        ('var', 'int'),
        ('pos_err', 'int'),
        ('mag_aca', 'float'),
        ('mag_aca_err', 'int'),
        ('mag_band', 'int'),
        ('pos_catid', 'int'),
        ('aspq1', 'int'),
        ('aspq2', 'int'),
        ('aspq3', 'int'),
        ('acqq1', 'int'),
        ('acqq2', 'int'),
        ('acqq4', 'int')],
    'temp': [
        ('n100_warm_frac', 'float'),
        ('tccd_mean', 'float'),
        ('tccd_max', 'float')],
    'bad': [
        ('known_bad', 'bool'),
        ('bad_comment', 'S15')],

    }


def get_options():
    parser = argparse.ArgumentParser(
        description="Update guide stats table")
    parser.add_argument("--check-missing",
                        action='store_true',
                        help="check for missing observations in table and reprocess")
    parser.add_argument("--obsid",
                        help="specific obsid to process.  Not required in regular update mode")
    parser.add_argument("--start",
                        help="start time for processing")
    parser.add_argument("--stop",
                        help="stop time for processing")
    parser.add_argument("--datafile",
                        default="gs.h5")
    opt = parser.parse_args()
    return opt


def _deltas_vs_obc_quat(vals, times, catalog):
    # Misalign is the identity matrix because this is the OBC quaternion
    aca_misalign = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    q_att = Quat(q=np.array([vals['AOATTQT1'],
                             vals['AOATTQT2'],
                             vals['AOATTQT3'],
                             vals['AOATTQT4']]).transpose())
    Ts = q_att.transform
    acqs = catalog

    R2A = 206264.81

    dy = {}
    dz = {}
    yag = {}
    zag = {}
    star_info = {}
    for slot in range(0, 8):
        if slot not in acqs['slot']:
            continue
        agasc_id = acqs[acqs['slot'] == slot][0]['id']
        if agasc_id is None:
            logger.info("No agasc id for slot {}, skipping".format(slot))
            continue
        try:
            # This is not perfect for star catalogs for agasc 1.4 and 1.5
            star = agasc.get_star(agasc_id, date=times[0],
                                  agasc_file='miniagasc_*',
                                  use_supplement=False)
        except:
            logger.info("agasc error on slot {}:{}".format(
                    slot, sys.exc_info()[0]))
            continue
        ra = star['RA_PMCORR']
        dec = star['DEC_PMCORR']
        star_pos_eci = radec_to_eci(ra, dec)
        d_aca = np.dot(np.dot(aca_misalign, Ts.transpose(0, 2, 1)),
                       star_pos_eci).transpose()
        yag[slot] = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
        zag[slot] = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
        dy[slot] = vals['AOACYAN{}'.format(slot)] - yag[slot]
        dz[slot] = vals['AOACZAN{}'.format(slot)] - zag[slot]
        star_info[slot] = star

    return dy, dz, star_info, yag, zag


def get_data(start, stop, obsid=None, starcheck=None):
    # Get telemetry
    msids = ['AOACASEQ', 'AOACQSUC', 'AOFREACQ', 'AOFWAIT', 'AOREPEAT',
             'AOACSTAT', 'AOACHIBK', 'AOFSTAR', 'AOFATTMD', 'AOACPRGS',
             'AOATUPST', 'AONSTARS', 'AOPCADMD', 'AORFSTR1', 'AORFSTR2',
             'AOATTQT1', 'AOATTQT2', 'AOATTQT3', 'AOATTQT4']
    per_slot = ['AOACQID', 'AOACFCT', 'AOIMAGE',
                'AOACMAG', 'AOACYAN', 'AOACZAN',
                'AOACICC', 'AOACIDP', 'AOACIIR', 'AOACIMS',
                'AOACIQB', 'AOACISP']
    slot_msids = [field + '%s' % slot
                  for field in per_slot
                  for slot in range(0, 8)]

    start_time = DateTime(start).secs
    stop_time = DateTime(stop)

    dat = fetch.MSIDset(msids + slot_msids,
                        start_time,
                        stop_time)
    if len(dat['AOACASEQ']) == 0:
        raise ValueError("No telemetry for obsid {}".format(obsid))
    # Interpolate the MSIDset onto the original time grid (which shouldn't do much)
    # but also remove all rows where any one msid has a bad value
    dat.interpolate(times=dat['AOACASEQ'].times, bad_union=True)
    eng_data = Table([col.vals for col in dat.values()], names=dat.keys())
    eng_data['times'] = dat.times

    times = eng_data['times']
    if starcheck is None:
        return eng_data, times, None
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    # Filter the catalog to be just guide stars
    catalog = catalog[(catalog['type'] == 'GUI') | (catalog['type'] == 'BOT')]
    # Get the position deltas relative to onboard solution
    dy, dz, star_info, yag, zag = _deltas_vs_obc_quat(eng_data, times, catalog)
    # And add the deltas to the table
    for slot in range(0, 8):
        if slot not in dy:
            continue
        eng_data['dy{}'.format(slot)] = dy[slot].data
        eng_data['dz{}'.format(slot)] = dz[slot].data
        eng_data['cat_yag{}'.format(slot)] = yag[slot]
        eng_data['cat_zag{}'.format(slot)] = zag[slot]
        cat_entry = catalog[catalog['slot'] == slot][0]
        dmag = eng_data['AOACMAG{}'.format(slot)] - cat_entry['mag']
        eng_data['dmag'] = dmag.data
    eng_data['time'] = times
    return eng_data, star_info


def consecutive(data, stepsize=1):
        return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


def calc_gui_stats(data, star_info):
    logger.info("calculating statistics")
    gui_stats = {}
    for slot in range(0, 8):
        if 'dy{}'.format(slot) not in data.colnames:
            continue
        stats = {}
        aoacfct = data['AOACFCT{}'.format(slot)]
        stats['n_samples'] = len(aoacfct)
        if len(aoacfct) == 0:
            gui_stats[slot] = stats
            continue
        stats['n_track'] = np.count_nonzero(aoacfct == 'TRAK')
        stats['f_track'] = stats['n_track'] / stats['n_samples']
        stats['f_racq'] = np.count_nonzero(aoacfct == 'RACQ') / stats['n_samples']
        stats['f_srch'] = np.count_nonzero(aoacfct == 'SRCH') / stats['n_samples']
        stats['f_none'] = np.count_nonzero(aoacfct == 'NONE') / stats['n_samples']
        if np.all(aoacfct != 'TRAK'):
            gui_stats[slot] = stats
            continue

        trak = data[aoacfct == 'TRAK']
        ok_flags = ((trak['AOACIIR{}'.format(slot)] == 'OK ')
                    & (trak['AOACISP{}'.format(slot)] == 'OK '))
        stats['n_kalman'] = np.count_nonzero(ok_flags)
        stats['no_track'] = (stats['n_samples'] - stats['n_track']) / stats['n_samples']
        stats['f_obc_bad'] = (stats['n_track'] - stats['n_kalman']) / stats['n_track']
        stats['f_common_col'] = np.count_nonzero(trak['AOACICC{}'.format(slot)] == 'ERR') / stats['n_track']
        stats['f_sat_pix'] = np.count_nonzero(trak['AOACISP{}'.format(slot)] == 'ERR') / stats['n_track']
        stats['f_def_pix'] = np.count_nonzero(trak['AOACIDP{}'.format(slot)] == 'ERR') / stats['n_track']
        stats['f_ion_rad'] = np.count_nonzero(trak['AOACIIR{}'.format(slot)] == 'ERR') / stats['n_track']
        stats['f_mult_star'] = np.count_nonzero(trak['AOACIMS{}'.format(slot)] == 'ERR') / stats['n_track']
        stats['f_quad_bound'] = np.count_nonzero(trak['AOACIQB{}'.format(slot)] == 'ERR') / stats['n_track']

        track_interv = consecutive(np.flatnonzero(
                data['AOACFCT{}'.format(slot)] == 'TRAK'))
        stats['n_track_interv'] = len(track_interv)
        track_interv_durations = np.array([len(interv) for interv in track_interv])
        stats['n_long_track_interv'] = np.count_nonzero(track_interv_durations > 60)

        not_track_interv = consecutive(np.flatnonzero(
                data['AOACFCT{}'.format(slot)] != 'TRAK'))
        not_track_interv_durations = np.array([len(interv) for interv in not_track_interv])
        stats['n_long_no_track_interv'] = np.count_nonzero(not_track_interv_durations > 60)

        stats['n_racq_interv'] = len(consecutive(np.flatnonzero(
                    data['AOACFCT{}'.format(slot)] == 'RACQ')))
        stats['n_srch_interv'] = len(consecutive(np.flatnonzero(
                    data['AOACFCT{}'.format(slot)] == 'SRCH')))
        stats['n_track_interv'] = len(consecutive(np.flatnonzero(
                    data['AOACFCT{}'.format(slot)] == 'TRAK')))


        # reduce this to just the samples that don't have IR or SP set and are in Kalman on guide stars
        # and are after the first 60 seconds
        kal = trak[ok_flags & (trak['AOACASEQ'] == 'KALM') & (trak['AOPCADMD'] == 'NPNT')
                   & (trak['AOFSTAR'] == 'GUID') & (trak['time'] > (data['time'][0] + 60))]
        dy = kal['dy{}'.format(slot)]
        dz = kal['dz{}'.format(slot)]
        # cheating here and ignoring spherical trig
        dr = (dy ** 2 + dz ** 2) ** .5
        stats['star_tracked'] = np.any(dr < 5.0)
        stats['spoiler_tracked'] = np.any(dr > 5.0)
        deltas = {'dy': dy, 'dz': dz, 'dr': dr}
        stats['dr_5th'] = np.percentile(deltas['dr'], 5)
        stats['dr_95th'] = np.percentile(deltas['dr'], 95)
        for ax in deltas:
            stats['{}_mean'.format(ax)] = np.mean(deltas[ax])
            stats['{}_std'.format(ax)] = np.std(deltas[ax])
            stats['{}_max'.format(ax)] = np.max(deltas[ax])
            stats['{}_min'.format(ax)] = np.min(deltas[ax])
        mag = kal['AOACMAG{}'.format(slot)]
        stats['aoacmag_min'] = np.min(mag)
        stats['aoacmag_mean'] = np.mean(mag)
        stats['aoacmag_max'] = np.max(mag)
        stats['aoacmag_std'] = np.std(mag)
        for perc in [5, 16, 50, 84, 95]:
            stats[f'aoacmag_{perc}th'] = np.percentile(mag, perc)
        stats['aoacyan_mean'] = np.mean(kal['AOACYAN{}'.format(slot)])
        stats['aoaczan_mean'] = np.mean(kal['AOACZAN{}'.format(slot)])

        for dist in ['0.3', '1', '3', '5']:
            stats['f_within_{}'.format(dist)] = np.count_nonzero(dr < float(dist)) / len(kal)
        stats['f_outside_5'] = np.count_nonzero(dr > 5) / len(kal)

        gui_stats[slot] = stats

    return gui_stats


def _get_obsids_to_update(check_missing=False, table_file=None, start=None, stop=None):
    if check_missing:
        last_tstart = start if start is not None else '2007:271:12:00:00'
        kadi_obsids = events.obsids.filter(start=last_tstart)
        try:
            h5 = tables.open_file(table_file, 'r')
            tbl = h5.root.data[:]
            h5.close()
        except:
            raise ValueError
        # get all obsids that aren't already in tbl
        obsids = [o.obsid for o in kadi_obsids if o.obsid not in tbl['obsid']]
    else:
        try:
            h5 = tables.open_file(table_file, 'r')
            tbl = h5.get_node('/', 'data')
            last_tstart = tbl.cols.kalman_tstart[tbl.colindexes['kalman_tstart'][-1]]
            h5.close()
        except:
            last_tstart = start if start is not None else  '2002:012:12:00:00'
        kadi_obsids = events.obsids.filter(start=last_tstart, stop=stop)
        # Skip the first obsid (as we already have it in the table)
        obsids = [o.obsid for o in kadi_obsids][1:]
    return obsids


def calc_stats(obsid):
    obspar = mica.archive.obspar.get_obspar(obsid)
    if not obspar:
        raise ValueError("No obspar for {}".format(obsid))
    manvr = None
    dwell = None
    try:
        manvrs = events.manvrs.filter(obsid=obsid, n_dwell__gt=0)
        dwells = events.dwells.filter(obsid=obsid)
        if dwells.count() == 1 and manvrs.count() == 0:
            # If there is more than one dwell for the manvr but they have
            # different obsids (unusual) so don't throw an overlapping interval kadi error
            # just get the maneuver to the attitude with this dwell
            dwell = dwells[0]
            manvr = dwell.manvr
        elif dwells.count() == 0:
            # If there's just nothing, that doesn't need an error here
            # and gets caught outside the try/except
            pass
        else:
            # Else just take the first matches from each
            manvr = manvrs[0]
            dwell = dwells[0]
    except ValueError:
        multi_manvr = events.manvrs.filter(start=obspar['tstart'] - 100000,
                                           stop=obspar['tstart'] + 100000)
        multi = multi_manvr.select_overlapping(events.obsids(obsid=obsid))
        deltas = [np.abs(m.tstart - obspar['tstart']) for m in multi]
        manvr = multi[np.argmin(deltas)]
        dwell = manvr.dwell_set.first()

    if not manvr or not dwell:
        raise ValueError("No manvr or dwell for {}".format(obsid))
    if not manvr.get_next():
        raise ValueError("No *next* manvr so can't calculate dwell")
    if not manvr.guide_start:
        raise ValueError("No guide transition for {}".format(obsid))
    if not manvr.kalman_start:
        raise ValueError("No Kalman transition for {}".format(obsid))


    logger.info("Found obsid manvr at {}".format(manvr.start))
    logger.info("Found dwell at {}".format(dwell.start))
    starcheck = get_starcheck_catalog_at_date(manvr.guide_start)
    if starcheck is None or 'cat' not in starcheck or not len(starcheck['cat']):
        raise ValueError('No starcheck catalog found for {}'.format(
                manvr.get_obsid()))
    starcat_time = DateTime(starcheck['cat']['mp_starcat_time'][0]).secs
    starcat_dtime = starcat_time - DateTime(manvr.start).secs
    # If it looks like the wrong starcheck by time, give up
    if abs(starcat_dtime) > 300:
        raise ValueError("Starcheck cat time delta is {}".format(starcat_dtime))
    if abs(starcat_dtime) > 30:
        logger.warning("Starcheck cat time delta of {} is > 30 sec".format(abs(starcat_dtime)))
    # The NPNT dwell should end when the next maneuver starts, but explicitly confirm via pcadmd
    pcadmd = fetch.Msid('AOPCADMD', manvr.kalman_start, manvr.get_next().tstart + 20)
    next_nman_start = pcadmd.times[pcadmd.vals != 'NPNT'][0]
    vals, star_info = get_data(start=manvr.kalman_start, stop=next_nman_start,
                               obsid=obsid, starcheck=starcheck)
    gui_stats = calc_gui_stats(vals, star_info)
    obsid_info = {'obsid': obsid,
                  'obi': obspar['obi_num'],
                  'kalman_datestart': manvr.kalman_start,
                  'kalman_tstart': DateTime(manvr.kalman_start).secs,
                  'npnt_tstop': DateTime(next_nman_start).secs,
                  'npnt_datestop': DateTime(next_nman_start).date,
                  'revision': STAT_VERSION}
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    guide_catalog = catalog[(catalog['type'] == 'GUI') | (catalog['type'] == 'BOT')]
    aacccdpt = fetch_sci.MSID('AACCCDPT', manvr.kalman_start, manvr.get_next().start)
    warm_threshold = 100.0
    tccd_mean = np.mean(aacccdpt.vals)
    tccd_max = np.max(aacccdpt.vals)
    warm_frac = dark_model.get_warm_fracs(warm_threshold, manvr.start, tccd_mean)
    temps = {'tccd_mean': tccd_mean, 'n100_warm_frac': warm_frac,
             'tccd_max': tccd_max}
    return obsid_info, gui_stats, star_info, guide_catalog, temps


def table_gui_stats(obsid_info, gui_stats, star_info, catalog, temp):
    logger.info("arranging stats into tabular data")
    cols = (GUIDE_COLS['obs'] + GUIDE_COLS['cat'] + GUIDE_COLS['stat']
            + GUIDE_COLS['agasc'] + GUIDE_COLS['temp'] + GUIDE_COLS['bad'])

    # Initialize all values to zero
    table = Table(np.zeros((1, 8), dtype=cols).flatten())
    # Set all columns with mag info to 99.0 initial value instead of zero
    for col in table.dtype.names:
        if 'aoacmag' in col:
            table[col] = 99.0

    for col in np.dtype(GUIDE_COLS['obs']).names:
        if col in obsid_info:
            table[col][:] = obsid_info[col]
    # Make a mask to identify 'missing' slots
    missing_slots = np.zeros(8, dtype=bool)
    for slot in range(0, 8):
        row = table[slot]
        if slot not in catalog['slot']:
            missing_slots[slot] = True
            continue
        for col in np.dtype(GUIDE_COLS['cat']).names:
            row[col] = catalog[catalog['slot'] == slot][0][col]
        for col in np.dtype(GUIDE_COLS['stat']).names:
            if col in gui_stats[slot]:
                row[col] = gui_stats[slot][col]
        if slot not in star_info:
            continue
        row['color'] = star_info[slot]['COLOR1']
        row['ra'] = star_info[slot]['RA_PMCORR']
        row['dec'] = star_info[slot]['DEC_PMCORR']
        for col in np.dtype(GUIDE_COLS['agasc']).names:
            if col in ['color', 'ra', 'dec']:
                continue
            row[col] = star_info[slot][col.upper()]
        row['tccd_mean'] = temp['tccd_mean']
        row['tccd_max'] = temp['tccd_max']
        row['n100_warm_frac'] = temp['n100_warm_frac']
        row['known_bad'] = False
        row['bad_comment'] = ''
    # Exclude any rows that are missing
    table = table[~missing_slots]
    return table


def _save_gui_stats(t, table_file=None):
    if table_file is None:
        return
    if not os.path.exists(table_file):
        cols = (GUIDE_COLS['obs'] + GUIDE_COLS['cat'] + GUIDE_COLS['stat']
                + GUIDE_COLS['agasc'] + GUIDE_COLS['temp'] + GUIDE_COLS['bad'])
        desc, byteorder = tables.descr_from_dtype(np.dtype(cols))
        filters = tables.Filters(complevel=5, complib='zlib')
        h5 = tables.open_file(table_file, 'a')
        tbl = h5.create_table('/', 'data', desc, filters=filters,
                             expectedrows=1e6)
        tbl.cols.obsid.create_index()
        tbl.cols.kalman_tstart.create_csindex()
        tbl.cols.agasc_id.create_index()
        h5.close()
        del h5
    h5 = tables.open_file(table_file, 'a')
    tbl = h5.get_node('/', 'data')
    have_obsid_coord = tbl.get_where_list(
        '(obsid == {}) & (obi == {})'.format(
            t[0]['obsid'], t[0]['obi']), sort=True)
    if len(have_obsid_coord):
        obsid_rec = tbl.read_coordinates(have_obsid_coord)
        if len(obsid_rec) != len(t):
            raise ValueError(
                "Could not update {}; different number of slots".format(
                    t[0]['obsid']))
        # preserve any 'known_bad' status
        for row in obsid_rec:
            slot = row['slot']
            t['known_bad'][t['slot'] == slot] = row['known_bad']
            t['bad_comment'][t['slot'] == slot] = row['bad_comment']
        tbl.modify_coordinates(have_obsid_coord, t.as_array())
    else:
        tbl.append(t.as_array())
    logger.info("saving stats to h5 table")
    tbl.flush()
    h5.flush()
    h5.close()


def update(opt):
    if opt.obsid:
        obsids = [int(opt.obsid)]
    else:
        obsids = _get_obsids_to_update(table_file=opt.datafile, check_missing=opt.check_missing,
                                        start=opt.start, stop=opt.stop)
    for obsid in obsids:

        logger.info("Processing obsid {}".format(obsid))
        try:
            obsid_info, gui_stats, star_info, guide_catalog, temp = calc_stats(obsid)
        except Exception as e:
            open(os.path.splitext(opt.datafile)[0] + '_skipped.dat', 'a').write(
                "{}: {}\n".format(obsid, e))
            logger.info("Skipping obsid {}: {}".format(obsid, e))
            continue
        if not len(gui_stats):
            open(os.path.splitext(opt.datafile)[0] + '_skipped.dat', 'a').write(
                "{}: No stats\n".format(obsid))
            logger.info("Skipping obsid {}, no stats determined".format(obsid))
            continue
        t = table_gui_stats(obsid_info, gui_stats, star_info, guide_catalog, temp)
        _save_gui_stats(t, opt.datafile)


def main():
    opt = get_options()
    update(opt)


if __name__ == '__main__':
    main()
