import os
import sys
import agasc
from kadi import events
from Ska.engarchive import fetch
from astropy.table import Table, Column
from Chandra.Time import DateTime
import mica.archive.obspar
import mica.starcheck
import tables
import numpy as np
import Ska.quatutil
import Ska.astro
from mica.quaternion import Quat
import logging

logger = logging.getLogger('star_stats')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

ID_DIST_LIMIT = 1.5

ACQ_COLS = {
    'obs': [
        ('obsid', 'int'),
        ('obi', 'int'),
        ('acq_start', 'S21'),
        ('guide_start', 'S21'),
        ('guide_tstart', 'float'),
        ('one_shot_length', 'float'),
        ('revision', 'S15')],
    'cat': [
        ('slot', 'int'),
        ('idx', 'int'),
        ('type', 'S5'),
        ('yang', 'float'),
        ('zang', 'float'),
        ('halfw', 'int'),
        ('mag', 'float')],
    'stat': [
        ('acqid', 'bool'),
        ('star_tracked', 'bool'),
        ('spoiler_tracked', 'bool'),
        ('img_func', 'S7'),
        ('n_trak_interv', 'int'),
        ('max_trak_cdy',  'float'),
        ('min_trak_cdy',  'float'),
        ('mean_trak_cdy', 'float'),
        ('max_trak_cdz',  'float'),
        ('min_trak_cdz',  'float'),
        ('mean_trak_cdz', 'float'),
        ('max_trak_mag',  'float'),
        ('min_trak_mag',  'float'),
        ('mean_trak_mag', 'float'),
        ('cdy', 'float'),
        ('cdz', 'float'),
        ('dy', 'float'),
        ('dz', 'float'),
        ('ion_rad', 'bool'),
        ('def_pix', 'bool'),
        ('mult_star', 'bool'),
        ('sat_pix', 'bool'),
        ('mag_obs', 'float'),
        ('yang_obs', 'float'),
        ('zang_obs', 'float')],
    'agasc': [
        ('agasc_id', 'int'),
        ('color1', 'float'),
        ('ra', 'float'),
        ('dec', 'float'),
        ('epoch', 'float'),
        ('pm_ra', 'int'),
        ('pm_dec', 'int'),
        ('var', 'int'),
        ('pos_err', 'int'),
        ('mag_aca', 'float'),
        ('mag_err', 'int'),
        ('mag_band', 'int'),
        ('pos_catid', 'int'),
        ('aspq1', 'int'),
        ('aspq2', 'int'),
        ('aspq3', 'int'),
        ('acqq1', 'int'),
        ('acqq2', 'int'),
        ('acqq4', 'int')],
    'bad': [
        ('known_bad', 'bool')]
    }

SKA = os.environ['SKA']
table_file = os.path.join(SKA, 'data', 'acq_stats', 'acq_stats.h5')


def deltas_vs_obc_quat(vals, times, catalog):
    # Ignore misalign
    aca_misalign = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    q_att = Quat(np.array([vals['AOATTQT1'],
                           vals['AOATTQT2'],
                           vals['AOATTQT3'],
                           vals['AOATTQT4']]).transpose())
    Ts = q_att.transform
    acqs = catalog

    R2A = 206264.81

    dy = {}
    dz = {}
    star_info = {}
    for slot in range(0, 8):
        if slot not in acqs['slot']:
            continue
        agasc_id = acqs[acqs['slot'] == slot][0]['id']
        if agasc_id is None:
            logger.info("No agasc id for slot {}, skipping".format(slot))
            continue
        try:
            star = agasc.get_star(agasc_id)
        except:
            logger.info("agasc error on slot {}:{}".format(
                    slot, sys.exc_info()[0]))
            continue
        ra = star['RA']
        dec = star['DEC']
        if (star['PM_RA'] != -9999 or star['PM_DEC'] != -9999):
            # Compute the multiplicative factor to convert from
            # the AGASC proper motion field to degrees.  The AGASC PM
            # is specified in milliarcsecs / year, so this is
            # dyear * (degrees / milliarcsec)
            dyear = ((DateTime(times[0]).secs
                     - DateTime("{}:001".format(int(star['EPOCH']))).secs)
                     / (86400 * 365.25))
            pm_to_degrees = dyear / (3600. * 1000.)
            if star['PM_RA'] != -9999:
                ra_scale = np.cos(np.radians(dec))
                ra = star['RA'] + star['PM_RA'] * pm_to_degrees / ra_scale
            if star['PM_DEC'] != -9999:
                dec = star['DEC'] + star['PM_DEC'] * pm_to_degrees
        star_pos_eci = Ska.quatutil.radec2eci(ra, dec)
        d_aca = np.dot(np.dot(aca_misalign, Ts.transpose(0, 2, 1)),
                       star_pos_eci).transpose()
        yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
        zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
        dy[slot] = vals['AOACYAN{}'.format(slot)] - yag
        dz[slot] = vals['AOACZAN{}'.format(slot)] - zag
        star_info[slot] = star

    return dy, dz, star_info


def get_delta_quat(eng_data, times, manvr):
    guide = eng_data[times >= DateTime(manvr.guide_start).secs - 1]
    first_guide = guide[0]
    # skip the duplicate and go to index 2
    second_guide = guide[2]
    q1 = Quat([first_guide['AOATTQT1'],
               first_guide['AOATTQT2'],
               first_guide['AOATTQT3'],
               first_guide['AOATTQT4']])
    q2 = Quat([second_guide['AOATTQT1'],
               second_guide['AOATTQT2'],
               second_guide['AOATTQT3'],
               second_guide['AOATTQT4']])
    dq = q2 / q1
    dot_q = np.sum(q1.q * q2.q)
    if dot_q > 1:
        dot_q = 1.0
    return dq, dot_q


def q_mult(q1, q2):
    mult = None
    if q1.ndim == 1:
        if q2.ndim != 1:
            mult = np.zeros((len(q2), 4))
        q1 = q1[np.newaxis]
    if q2.ndim == 1:
        q2 = q2[np.newaxis]
    if mult is None:
        mult = np.zeros((len(q1), 4))
    mult[:, 0] = (q1[:, 3] * q2[:, 0] - q1[:, 2] * q2[:, 1]
                  + q1[:, 1] * q2[:, 2] + q1[:, 0] * q2[:, 3])
    mult[:, 1] = (q1[:, 2] * q2[:, 0] + q1[:, 3] * q2[:, 1]
                  - q1[:, 0] * q2[:, 2] + q1[:, 1] * q2[:, 3])
    mult[:, 2] = (-q1[:, 1] * q2[:, 0] + q1[:, 0] * q2[:, 1]
                  + q1[:, 3] * q2[:, 2] + q1[:, 2] * q2[:, 3])
    mult[:, 3] = (-q1[:, 0] * q2[:, 0] - q1[:, 1] * q2[:, 1]
                  - q1[:, 2] * q2[:, 2] + q1[:, 3] * q2[:, 3])
    return Quat(mult)


def search_agasc(yang, zang, field_agasc, q_aca):
    """
    Search the retrieved agasc region for a star at the specified
    yang, zang, and return the star if there is a match.

    :param yang:
    :param zang:
    :param field_agasc: the retrieved agasc star field
    :param q_aca: pointing quaternion for obsid
    :rtype: recarray of the matching star or None
    """

    for agasc_star in field_agasc:
        ra, dec = Ska.quatutil.yagzag2radec(
            yang * 1.0 / 3600,
            zang * 1.0 / 3600,
            q_aca)
        # 3600*(sph_dist in degrees) for arcseconds
        dist = 3600 * Ska.astro.sph_dist(agasc_star['RA_PMCORR'],
                                         agasc_star['DEC_PMCORR'],
                                         ra, dec)
        if dist <= ID_DIST_LIMIT:
            return agasc_star

    return None


def get_modern_data(manvr, dwell, starcheck):
    if 'cat' not in starcheck:
        raise ValueError('No starcheck catalog found for {}'.format(manvr.get_obsid()))
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    # Filter the catalog to be just acquisition stars
    catalog = catalog[(catalog['type'] == 'ACQ') | (catalog['type'] == 'BOT')]
    slot_for_pos = [cat_row['slot'] for cat_row in catalog]
    pos_for_slot = dict([(slot, idx) for idx, slot in enumerate(slot_for_pos)])
    # Also, save out the starcheck index for each slot for later
    index_for_slot = dict([(cat_row['slot'], cat_row['idx'])
                           for cat_row in catalog])

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

    start_time = DateTime(manvr.acq_start).secs
    stop_time = DateTime(dwell.start).secs + 100
    raw_eng_data = fetch.MSIDset(msids + slot_msids,
                                 start_time,
                                 stop_time,
                                 filter_bad=True)
    eng_data = Table([raw_eng_data[col].vals for col in msids],
                     names=msids)
    for field in slot_msids:
        eng_data.add_column(
            Column(
                name=field, data=raw_eng_data[field].vals))
        times = Table([raw_eng_data['AOACASEQ'].times],
                      names=['time'])
    if not len(eng_data['AOACASEQ']):
        raise ValueError("No telemetry for obsid {}".format(manvr.get_obsid()))

    # Estimate the offsets from the expected catalog positions
    dy, dz, star_info = deltas_vs_obc_quat(eng_data, times['time'], catalog)
    # And add the deltas to the table
    for slot in range(0, 8):
        if slot not in dy:
            continue
        eng_data.add_column(Column(name='dy{}'.format(slot),
                                   data=dy[slot].data))
        eng_data.add_column(Column(name='dz{}'.format(slot),
                                   data=dz[slot].data))
        cat_entry = catalog[catalog['slot'] == slot][0]
        dmag = eng_data['AOACMAG{}'.format(slot)] - cat_entry['mag']
        eng_data.add_column(Column(name='dmag{}'.format(slot),
                                   data=dmag.data))

    # Get the one-shot delta quaternion and the dot product of the deltas
    delta_quat, dot_q = get_delta_quat(eng_data, times['time'], manvr)
    one_shot_length = np.degrees(2 * np.arccos(dot_q))
    one_shot_length = np.min([one_shot_length, 360 - one_shot_length])
    one_shot_length = one_shot_length * 3600

    # Update a copy of the telemetry structure with quaternions
    # corrected by the one-shot delta
    corr_eng_data = eng_data.copy()
    uncorr_times = (times['time'] < DateTime(manvr.guide_start).secs + 1.0)
    q_orig = Quat(np.array([eng_data[uncorr_times]['AOATTQT1'],
                            eng_data[uncorr_times]['AOATTQT2'],
                            eng_data[uncorr_times]['AOATTQT3'],
                            eng_data[uncorr_times]['AOATTQT4']]).transpose())
    q_corr = q_mult(delta_quat.q, q_orig.q)
    corr_eng_data['AOATTQT1'][uncorr_times] = q_corr.q.transpose()[0]
    corr_eng_data['AOATTQT2'][uncorr_times] = q_corr.q.transpose()[1]
    corr_eng_data['AOATTQT3'][uncorr_times] = q_corr.q.transpose()[2]
    corr_eng_data['AOATTQT4'][uncorr_times] = q_corr.q.transpose()[3]
    corr_dy, corr_dz, si = deltas_vs_obc_quat(corr_eng_data, times['time'], catalog)
    # delete the now-extra copy of the data
    del corr_eng_data
    # And add the corrected deltas to the table
    for slot in range(0, 8):
        if slot not in corr_dy:
            continue
        eng_data.add_column(Column(name='corr_dy{}'.format(slot),
                                   data=corr_dy[slot].data))
        eng_data.add_column(Column(name='corr_dz{}'.format(slot),
                                   data=corr_dz[slot].data))

    # Also add the acquisition id in a useful way
    for slot in range(0, 8):
        if slot not in pos_for_slot:
            continue
        eng_data.add_column(
            Column(
                name='POS_ACQID{}'.format(slot),
                data=eng_data['AOACQID{}'.format(pos_for_slot[slot])]))

    return eng_data, times['time'], one_shot_length, star_info


def calc_acq_stats(manvr, vals, times):
    logger.info("calculating statistics")
    acq_stats = {}
    for slot in range(0, 8):
        if 'dy{}'.format(slot) not in vals.colnames:
            continue
        stats = {}
        guide_times = (times >= DateTime(manvr.guide_start).secs - 1)
        stats['acqid'] = vals['POS_ACQID{}'.format(slot)][guide_times][0] == 'ID  '
        acq_times = ((times > DateTime(manvr.acq_start).secs)
                     & (times < DateTime(manvr.guide_start).secs + 1))
        acq_data = vals[acq_times]
        aoacfct = acq_data['AOACFCT{}'.format(slot)]
        # Does it look like the desired star was tracked?
        stats['star_tracked'] = False
        stats['spoiler_tracked'] = False
        if np.any(aoacfct == 'TRAK'):
            trak = acq_data[aoacfct == 'TRAK']
            corr_dy = trak['corr_dy{}'.format(slot)]
            corr_dz = trak['corr_dz{}'.format(slot)]
            # cheating here and ignoring spherical trig
            corr_dr = (corr_dy ** 2 + corr_dz ** 2) ** .5
            if np.min(corr_dr) < 5.0:
                stats['star_tracked'] = True
            if np.max(corr_dr) >= 5.0:
                stats['spoiler_tracked'] = True
            stats['mean_trak_cdy'] = np.mean(corr_dy)
            stats['min_trak_cdy'] = np.min(corr_dy)
            stats['max_trak_cdy'] = np.max(corr_dy)
            stats['mean_trak_cdz'] = np.mean(corr_dz)
            stats['min_trak_cdz'] = np.min(corr_dz)
            stats['max_trak_cdz'] = np.max(corr_dz)
            stats['mean_trak_mag'] = np.mean(trak['AOACMAG{}'.format(slot)])
            stats['min_trak_mag'] = np.min(trak['AOACMAG{}'.format(slot)])
            stats['max_trak_mag'] = np.max(trak['AOACMAG{}'.format(slot)])
            # Did we try to find the star more than once?
            stats['n_trak_interv'] = 0
            for i, stat in enumerate(aoacfct):
                if (i + 1) < len(aoacfct):
                    if stat != 'TRAK' and aoacfct[i + 1] == 'TRAK':
                        stats['n_trak_interv'] += 1

        stats['img_func'] = None
        tracked_at_trans = aoacfct[-2] == 'TRAK'
        # I think it makes slightly more sense to use the
        # "-2" values of these because the quaternions were updated
        # but the star positions weren't at -1
        if tracked_at_trans:
            stats['dy'] = acq_data[-2]['dy{}'.format(slot)]
            stats['dz'] = acq_data[-2]['dz{}'.format(slot)]
            stats['cdy'] = acq_data[-2]['corr_dy{}'.format(slot)]
            stats['cdz'] = acq_data[-2]['corr_dz{}'.format(slot)]
            stats['mag_obs'] = acq_data[-2]['AOACMAG{}'.format(slot)]
            stats['yang_obs'] = acq_data[-2]['AOACYAN{}'.format(slot)]
            stats['zang_obs'] = acq_data[-2]['AOACZAN{}'.format(slot)]
            if ((stats['cdy'] ** 2 + stats['cdz'] ** 2) ** .5) < 5.0:
                stats['img_func'] = 'star'
            else:
                stats['img_func'] = 'spoiler'
        else:
            stats['img_func'] = aoacfct[-2]

        guide_data = vals[guide_times]
        flag_map = {'AOACIDP': 'def_pix',
                    'AOACIIR': 'ion_rad',
                    'AOACIMS': 'mult_star',
                    'AOACISP': 'sat_pix'}
        for flag in ['AOACIDP', 'AOACISP', 'AOACIMS', 'AOACIIR']:
            if guide_data['{}{}'.format(flag, slot)][0] == 'ERR':
                stats[flag_map[flag]] = True
        acq_stats[slot] = stats
    return acq_stats


def get_obsids_to_update():
    try:
        h5 = tables.openFile(table_file, 'r')
        tbl = h5.getNode('/', 'data')
        last_tstart = tbl.cols.guide_tstart[tbl.colindexes['guide_tstart'][-1]]
        h5.close()
    except:
        last_tstart = '1998:001'
    kadi_obsids = events.obsids.filter(start=last_tstart)
    obsids = [o.obsid for o in kadi_obsids]
    # Skip the first obsid (as we already have it in the table)
    return obsids[1:]


def get_stats(obsid):
    obspar = mica.archive.obspar.get_obspar(obsid)
    if not obspar:
        raise ValueError("No obspar for {}".format(obsid))
    manvr = None
    dwell = None
    try:
        manvrs = events.manvrs.filter(obsid=obsid)
        dwells = events.dwells.filter(obsid=obsid)
        if manvrs.count() and dwells.count() == 1:
            manvr = manvrs[0]
            dwell = dwells[0]
    except ValueError:
        multi_manvr = events.manvrs.filter(start=obspar['tstart'] - 10000,
                                           stop=obspar['tstart'] + 10000)
        multi = multi_manvr.select_overlapping(events.obsids(obsid=obsid))
        deltas = [np.abs(m.tstart - obspar['tstart']) for m in multi]
        manvr = multi[np.argmin(deltas)]
        dwell = manvr.dwell_set.first()

    if not manvr or not dwell:
        raise ValueError("No manvr or dwell for {}".format(obsid))

    logger.info("Found obsid manvr at {}".format(manvr.start))
    logger.info("Found dwell at {}".format(dwell.start))
    acq_start = manvr.acq_start
    guide_start = manvr.guide_start
    try:
        starcheck = mica.starcheck.get_starcheck_catalog(int(obsid))
    except:
        raise ValueError("Problem looking up starcheck for {}".format(obsid))
    vals, times, one_shot, star_info = get_modern_data(manvr, dwell, starcheck)

    acq_stats = calc_acq_stats(manvr, vals, times)
    obsid_info = {'obsid': obsid,
                  'obi': obspar['obi_num'],
                  'acq_start': acq_start,
                  'guide_start':  guide_start,
                  'guide_tstart': DateTime(guide_start).secs,
                  'one_shot_length': one_shot,
                  'revision': '1.0'}
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    # Filter the catalog to be just acquisition stars
    catalog = catalog[(catalog['type'] == 'ACQ') | (catalog['type'] == 'BOT')]
    return obsid_info, acq_stats, star_info, catalog


def table_acq_stats(obsid_info, acq_stats, star_info, catalog):
    logger.info("arranging stats into tabular data")
    cols = (ACQ_COLS['obs'] + ACQ_COLS['cat'] + ACQ_COLS['stat']
            + ACQ_COLS['agasc'] + ACQ_COLS['bad'])
    table = Table(np.zeros((1, 8), dtype=cols).flatten())
    for col in np.dtype(ACQ_COLS['obs']).names:
        if col in obsid_info:
            table[col][:] = obsid_info[col]
    # Make a mask to identify 'missing' slots
    missing_slots = np.zeros(8, dtype=bool)
    for slot in range(0, 8):
        row = table[slot]
        if slot not in catalog['slot'] or slot not in acq_stats:
            missing_slots[slot] = True
            continue
        for col in np.dtype(ACQ_COLS['cat']).names:
            row[col] = catalog[catalog['slot'] == slot][0][col]
        for col in np.dtype(ACQ_COLS['stat']).names:
            if col in acq_stats[slot]:
                row[col] = acq_stats[slot][col]
        if slot not in star_info:
            continue
        for col in np.dtype(ACQ_COLS['agasc']).names:
            row[col] = star_info[slot][col.upper()]
        row['known_bad'] = False
    # Exclude any rows that are missing
    table = table[~missing_slots]
    return table


def save_acq_stats(t):
    if not os.path.exists(table_file):
        cols = (ACQ_COLS['obs'] + ACQ_COLS['cat'] + ACQ_COLS['stat']
                + ACQ_COLS['agasc'] + ACQ_COLS['bad'])
        desc, byteorder = tables.descr_from_dtype(np.dtype(cols))
        filters = tables.Filters(complevel=5, complib='zlib')
        h5 = tables.openFile(table_file, 'a')
        tbl = h5.createTable('/', 'data', desc, filters=filters,
                             expectedrows=1e6)
        tbl.cols.obsid.createIndex()
        tbl.cols.guide_tstart.createCSIndex()
        h5.close()
        del h5
    h5 = tables.openFile(table_file, 'a')
    tbl = h5.getNode('/', 'data')
    have_obsid_coord = tbl.getWhereList(
        '(obsid == {}) & (obi == {})'.format(
            t[0]['obsid'], t[0]['obi']), sort=True)
    if len(have_obsid_coord):
        obsid_rec = tbl.readCoordinates(have_obsid_coord)
        if len(obsid_rec) != len(t):
            raise ValueError(
                "Could not update {}; different number of slots".format(
                    t[0]['obsid']))
        # preserve any 'known_bad' status
        for row in obsid_rec:
            slot = row['slot']
            t['known_bad'][t['slot'] == slot] = row['known_bad']
        tbl.modifyCoordinates(have_obsid_coord, t._data)
    else:
        tbl.append(t._data)
    logger.info("saving stats to h5 table")
    tbl.flush()
    h5.flush()
    h5.close()


def update():
    obsids = get_obsids_to_update()
    for obsid in obsids:
        logger.info("Processing obsid {}".format(obsid))
        try:
            obsid_info, acq_stats, star_info, catalog = get_stats(obsid)
        except ValueError:
            logger.info("Skipping obsid {}".format(obsid))
            continue
        if not len(acq_stats):
            logger.info("Skipping obsid {}".format(obsid))
            continue
        t = table_acq_stats(obsid_info, acq_stats, star_info, catalog)
        save_acq_stats(t)

def main():
    update()


if __name__ == '__main__':
    main()