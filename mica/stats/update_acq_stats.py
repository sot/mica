# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import sys
import argparse
import agasc
from kadi import events
from Ska.engarchive import fetch, fetch_sci
from astropy.table import Table, Column
from Chandra.Time import DateTime
import mica.archive.obspar
from mica.starcheck import get_starcheck_catalog, get_starcheck_catalog_at_date
import tables
import numpy as np
import Ska.astro
from Quaternion import Quat
from chandra_aca import dark_model
from chandra_aca.transform import radec_to_eci, yagzag_to_radec
import logging
import smtplib
from email.mime.text import MIMEText

import warnings

# Ignore known numexpr.necompiler and table.conditions warning
warnings.filterwarnings(
    'ignore',
    message="using `oa_ndim == 0` when `op_axes` is NULL is deprecated.*",
    category=DeprecationWarning,
)


logger = logging.getLogger('star_stats')
logger.setLevel(logging.INFO)
if not len(logger.handlers):
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
        ('revision', 'S15'),
    ],
    'cat': [
        ('slot', 'int'),
        ('idx', 'int'),
        ('type', 'S5'),
        ('yang', 'float'),
        ('zang', 'float'),
        ('halfw', 'int'),
        ('mag', 'float'),
    ],
    'stat': [
        ('acqid', 'bool'),
        ('star_tracked', 'bool'),
        ('spoiler_tracked', 'bool'),
        ('img_func', 'S7'),
        ('n_trak_interv', 'int'),
        ('max_trak_cdy', 'float'),
        ('min_trak_cdy', 'float'),
        ('mean_trak_cdy', 'float'),
        ('max_trak_cdz', 'float'),
        ('min_trak_cdz', 'float'),
        ('mean_trak_cdz', 'float'),
        ('max_trak_mag', 'float'),
        ('min_trak_mag', 'float'),
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
        ('zang_obs', 'float'),
    ],
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
        ('acqq4', 'int'),
    ],
    'temp': [('n100_warm_frac', 'float'), ('ccd_temp', 'float')],
    'bad': [('known_bad', 'bool'), ('bad_comment', 'S15')],
}

SKA = os.environ['SKA']
table_file = os.path.join(SKA, 'data', 'acq_stats', 'acq_stats.h5')
DEFAULT_EMAIL = ['aca@cfa.harvard.edu']


def get_options():
    parser = argparse.ArgumentParser(description="Update acq stats table")
    parser.add_argument(
        "--check-missing",
        action='store_true',
        help="check for missing observations in table and reprocess",
    )
    parser.add_argument(
        "--obsid",
        help="specific obsid to process.  Not required in regular update mode",
    )
    parser.add_argument(
        '--email',
        action="append",
        help="email warning recipient, specify multiple times "
        + "for multiple recipients",
    )
    opt = parser.parse_args()
    return opt


def _deltas_vs_obc_quat(vals, times, catalog):
    # Ignore misalign
    aca_misalign = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    q_att = Quat(
        q=np.array(
            [vals['AOATTQT1'], vals['AOATTQT2'], vals['AOATTQT3'], vals['AOATTQT4']]
        ).transpose()
    )
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
            star = agasc.get_star(
                agasc_id, date=times[0], agasc_file='agasc*', use_supplement=False
            )
        except Exception as e:
            logger.info(f"agasc error on slot {slot} id {agasc_id} err {e} {type(e)}")
            continue
        ra = star['RA_PMCORR']
        dec = star['DEC_PMCORR']
        star_pos_eci = radec_to_eci(ra, dec)
        d_aca = np.dot(
            np.dot(aca_misalign, Ts.transpose(0, 2, 1)), star_pos_eci
        ).transpose()
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
    q1 = Quat(
        [
            first_guide['AOATTQT1'],
            first_guide['AOATTQT2'],
            first_guide['AOATTQT3'],
            first_guide['AOATTQT4'],
        ]
    )
    q2 = Quat(
        [
            second_guide['AOATTQT1'],
            second_guide['AOATTQT2'],
            second_guide['AOATTQT3'],
            second_guide['AOATTQT4'],
        ]
    )
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
    mult[:, 0] = (
        q1[:, 3] * q2[:, 0]
        - q1[:, 2] * q2[:, 1]
        + q1[:, 1] * q2[:, 2]
        + q1[:, 0] * q2[:, 3]
    )
    mult[:, 1] = (
        q1[:, 2] * q2[:, 0]
        + q1[:, 3] * q2[:, 1]
        - q1[:, 0] * q2[:, 2]
        + q1[:, 1] * q2[:, 3]
    )
    mult[:, 2] = (
        -q1[:, 1] * q2[:, 0]
        + q1[:, 0] * q2[:, 1]
        + q1[:, 3] * q2[:, 2]
        + q1[:, 2] * q2[:, 3]
    )
    mult[:, 3] = (
        -q1[:, 0] * q2[:, 0]
        - q1[:, 1] * q2[:, 1]
        - q1[:, 2] * q2[:, 2]
        + q1[:, 3] * q2[:, 3]
    )
    return Quat(q=mult)


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
        ra, dec = yagzag_to_radec(yang * 1.0 / 3600, zang * 1.0 / 3600, q_aca)
        # 3600*(sph_dist in degrees) for arcseconds
        dist = 3600 * Ska.astro.sph_dist(
            agasc_star['RA_PMCORR'], agasc_star['DEC_PMCORR'], ra, dec
        )
        if dist <= ID_DIST_LIMIT:
            return agasc_star

    return None


def get_modern_data(manvr, dwell, starcheck):
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    # Filter the catalog to be just acquisition stars
    catalog = catalog[(catalog['type'] == 'ACQ') | (catalog['type'] == 'BOT')]
    slot_for_pos = [cat_row['slot'] for cat_row in catalog]
    pos_for_slot = dict([(slot, idx) for idx, slot in enumerate(slot_for_pos)])
    # Also, save out the starcheck index for each slot for later
    index_for_slot = dict([(cat_row['slot'], cat_row['idx']) for cat_row in catalog])

    # Get telemetry
    msids = [
        'AOACASEQ',
        'AOACQSUC',
        'AOFREACQ',
        'AOFWAIT',
        'AOREPEAT',
        'AOACSTAT',
        'AOACHIBK',
        'AOFSTAR',
        'AOFATTMD',
        'AOACPRGS',
        'AOATUPST',
        'AONSTARS',
        'AOPCADMD',
        'AORFSTR1',
        'AORFSTR2',
        'AOATTQT1',
        'AOATTQT2',
        'AOATTQT3',
        'AOATTQT4',
    ]
    per_slot = [
        'AOACQID',
        'AOACFCT',
        'AOIMAGE',
        'AOACMAG',
        'AOACYAN',
        'AOACZAN',
        'AOACICC',
        'AOACIDP',
        'AOACIIR',
        'AOACIMS',
        'AOACIQB',
        'AOACISP',
    ]
    slot_msids = [field + '%s' % slot for field in per_slot for slot in range(0, 8)]

    start_time = DateTime(manvr.acq_start).secs
    stop_time = DateTime(dwell.start).secs + 100
    raw_eng_data = fetch.MSIDset(
        msids + slot_msids, start_time, stop_time, filter_bad=True
    )
    eng_data = Table([raw_eng_data[col].vals for col in msids], names=msids)
    for field in slot_msids:
        eng_data.add_column(Column(name=field, data=raw_eng_data[field].vals))
        times = Table([raw_eng_data['AOACASEQ'].times], names=['time'])
    if not len(eng_data['AOACASEQ']):
        raise ValueError("No telemetry for obsid {}".format(manvr.get_obsid()))

    # Estimate the offsets from the expected catalog positions
    dy, dz, star_info = _deltas_vs_obc_quat(eng_data, times['time'], catalog)
    # And add the deltas to the table
    for slot in range(0, 8):
        if slot not in dy:
            continue
        eng_data.add_column(Column(name='dy{}'.format(slot), data=dy[slot].data))
        eng_data.add_column(Column(name='dz{}'.format(slot), data=dz[slot].data))
        cat_entry = catalog[catalog['slot'] == slot][0]
        dmag = eng_data['AOACMAG{}'.format(slot)] - cat_entry['mag']
        eng_data.add_column(Column(name='dmag{}'.format(slot), data=dmag.data))

    # Get the one-shot delta quaternion and the dot product of the deltas
    delta_quat, dot_q = get_delta_quat(eng_data, times['time'], manvr)
    one_shot_length = np.degrees(2 * np.arccos(dot_q))
    one_shot_length = np.min([one_shot_length, 360 - one_shot_length])
    one_shot_length = one_shot_length * 3600

    # Update a copy of the telemetry structure with quaternions
    # corrected by the one-shot delta
    corr_eng_data = eng_data.copy()
    uncorr_times = times['time'] < DateTime(manvr.guide_start).secs + 1.0
    q_orig = Quat(
        q=np.array(
            [
                eng_data[uncorr_times]['AOATTQT1'],
                eng_data[uncorr_times]['AOATTQT2'],
                eng_data[uncorr_times]['AOATTQT3'],
                eng_data[uncorr_times]['AOATTQT4'],
            ]
        ).transpose()
    )
    q_corr = q_mult(delta_quat.q, q_orig.q)
    corr_eng_data['AOATTQT1'][uncorr_times] = q_corr.q.transpose()[0]
    corr_eng_data['AOATTQT2'][uncorr_times] = q_corr.q.transpose()[1]
    corr_eng_data['AOATTQT3'][uncorr_times] = q_corr.q.transpose()[2]
    corr_eng_data['AOATTQT4'][uncorr_times] = q_corr.q.transpose()[3]
    corr_dy, corr_dz, si = _deltas_vs_obc_quat(corr_eng_data, times['time'], catalog)
    # delete the now-extra copy of the data
    del corr_eng_data
    # And add the corrected deltas to the table
    for slot in range(0, 8):
        if slot not in corr_dy:
            continue
        eng_data.add_column(
            Column(name='corr_dy{}'.format(slot), data=corr_dy[slot].data)
        )
        eng_data.add_column(
            Column(name='corr_dz{}'.format(slot), data=corr_dz[slot].data)
        )

    # Also add the acquisition id in a useful way
    for slot in range(0, 8):
        if slot not in pos_for_slot:
            continue
        eng_data.add_column(
            Column(
                name='POS_ACQID{}'.format(slot),
                data=eng_data['AOACQID{}'.format(pos_for_slot[slot])],
            )
        )

    return eng_data, times['time'], one_shot_length, star_info


def calc_acq_stats(manvr, vals, times):
    logger.info("calculating statistics")
    acq_stats = {}
    guide_times = times >= DateTime(manvr.guide_start).secs - 1
    acq_times = (times > DateTime(manvr.acq_start).secs) & (
        times < DateTime(manvr.guide_start).secs + 1
    )
    acq_data = vals[acq_times]
    for slot in range(0, 8):
        if 'dy{}'.format(slot) not in vals.colnames:
            continue
        stats = {}
        stats['acqid'] = vals['POS_ACQID{}'.format(slot)][guide_times][0] == 'ID  '
        aoacfct = acq_data['AOACFCT{}'.format(slot)]
        # Does it look like the desired star was tracked?
        stats['star_tracked'] = False
        stats['spoiler_tracked'] = False
        if np.any(aoacfct == 'TRAK'):
            trak = acq_data[aoacfct == 'TRAK']
            corr_dy = trak['corr_dy{}'.format(slot)]
            corr_dz = trak['corr_dz{}'.format(slot)]
            # cheating here and ignoring spherical trig
            corr_dr = (corr_dy**2 + corr_dz**2) ** 0.5
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
            if ((stats['cdy'] ** 2 + stats['cdz'] ** 2) ** 0.5) < 5.0:
                stats['img_func'] = 'star'
            else:
                stats['img_func'] = 'spoiler'
        else:
            stats['img_func'] = aoacfct[-2]

        guide_data = vals[guide_times]
        flag_map = {
            'AOACIDP': 'def_pix',
            'AOACIIR': 'ion_rad',
            'AOACIMS': 'mult_star',
            'AOACISP': 'sat_pix',
        }
        for flag in ['AOACIDP', 'AOACISP', 'AOACIMS', 'AOACIIR']:
            if guide_data['{}{}'.format(flag, slot)][0] == 'ERR':
                stats[flag_map[flag]] = True
        acq_stats[slot] = stats
    return acq_stats


def _get_obsids_to_update(check_missing=False):
    if check_missing:
        last_tstart = '2007:271:12:00:00'
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
            last_tstart = tbl.cols.guide_tstart[tbl.colindexes['guide_tstart'][-1]]
            h5.close()
        except:
            last_tstart = '2002:012:12:00:00'
        kadi_obsids = events.obsids.filter(start=last_tstart)
        # Skip the first obsid (as we already have it in the table)
        obsids = [o.obsid for o in kadi_obsids][1:]
    return obsids


def warn_on_acq_anom(acqs, emails):
    """
    Log and email warning about any acquisition stars with observed positions outside the
    expected search box (plus a pad).  For acquisition stars with tracked positions in the wrong
    search box, note the box (classic acquisition anomaly).

    :param acqs: astropy table of acquisition stats, includes expected and observed position,
                 obsid, slot, and halfw used in this method.
    :param emails: list of addresses to receive email warning if any are generated
    """
    # Find tracked objects in the acq stats table outside the intended search box plus padding

    # Note that dy/dz are observed yag/zag (t_guide) - predicted yag/zag (t_guide) using AOATTQT
    # (estimated attitude). Observed yag/zag are from AOAC{Y,Z}AN, and t_guide is the time of the
    # first sample with AOACASEQ = 'GUID'. t_guide is the same as manvrs.guide_start in kadi.
    # The one-shot attitude update occurs in telemetry on the sample after the GUID transition,
    # so the estimated attitude for dy/dz gives a reasonable approximation of the OBC estimated
    # attitude at the time of commanding the search boxes. (It would be more accurate to use the
    # time of transition to acquisition, but this approximation is at least closer than using
    # catalog yag/zag.)
    box_pad = 16  # arcsecs
    anom_match = (
        (acqs['img_func'] != 'NONE')
        & (acqs['img_func'] != 'SRCH')
        & (
            (np.abs(acqs['dy']) >= (acqs['halfw'] + box_pad))
            | (np.abs(acqs['dz']) >= (acqs['halfw'] + box_pad))
        )
    )
    for anom in acqs[anom_match]:
        # Check to see if the star is actually found in another box.
        other_box_match = (
            np.abs(anom['yang_obs'] - acqs['yang']) <= (acqs['halfw'] + box_pad)
        ) & (np.abs(anom['zang_obs'] - acqs['zang']) <= (acqs['halfw'] + box_pad))
        if np.any(other_box_match):
            text = (
                "Acquisition Anomaly.  Star for slot {} actually in box {} \n".format(
                    anom['slot'], acqs[other_box_match][0]['slot']
                )
            )
        else:
            text = "Does not appear to be classic star-in-wrong-box anomaly\n"
        # Make a dictionary of the anom record for use in string formatting
        output_dict = {col: anom[col] for col in anom.colnames}
        output_dict['dy'] = anom['dy']
        output_dict['dz'] = anom['dz']
        text += """Large Deviation from Expected ACQ Star Position in {obsid}
      Slot {slot} Expected (Y-Pos, Z-Pos) = ({yang:.1f}, {zang:.1f})
      Slot {slot} Observed (Y-Pos, Z-Pos) = ({yang_obs:.1f}, {zang_obs:.1f})
      Halfwidth {halfw:03d}        (dy, dz) = ({dy:.1f}, {dz:.1f})

      Expected here is catalog Y-Pos/Z-pos.  dy, dz calculation corrects these for estimated attitude.
""".format(**output_dict)
        # Log and Send message for slot.  Obsid can have more than one email
        logger.warning(text)
        msg = MIMEText(text)
        msg['From'] = 'aca@head.cfa.harvard.edu'
        msg['Subject'] = "Acq Anomaly: Obsid {} (mica processing)".format(anom['obsid'])
        msg['To'] = ",".join(emails)
        s = smtplib.SMTP('head.cfa.harvard.edu')
        s.sendmail('aca@head.cfa.harvard.edu', emails, msg.as_string())
        s.quit()


def calc_stats(obsid):
    obspar = mica.archive.obspar.get_obspar(obsid, version='last')
    if not obspar:
        raise ValueError("No obspar for {}".format(obsid))
    manvr = None
    dwell = None
    try:
        manvrs = events.manvrs.filter(obsid=obsid, n_dwell__gt=0)
        if len(manvrs) == 0:
            raise ValueError
        dwells = events.dwells.filter(obsid=obsid)
        # Use the last manvr and the first dwell
        manvr = manvrs[manvrs.count() - 1]
        dwell = dwells[0]
    except ValueError:
        multi_manvr = events.manvrs.filter(
            start=obspar['tstart'] - 10000, stop=obspar['tstart'] + 10000
        )
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
        starcheck = get_starcheck_catalog_at_date(manvr.acq_start)
    except Exception:
        # No matching observations for some known observations with problems. Use
        # hard-coded input for get_starcheck_catalog.
        if obsid in [1966]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2002/JAN1202/oflsa')
        elif obsid in [3105, 2741, 61334, 61333, 61332, 61331, 3358, 3357]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2002/JAN2802/oflsd/')
        elif obsid in [61261]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2002/MAR1902/oflsa/')
        elif obsid in [3471, 3086, 61250, 61249, 3094, 3066, 3115, 2833, 3464, 3175]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2002/MAR2502/oflsb/')
        elif obsid in [3663, 61185, 61184, 3392, 61183]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2002/MAY2402/oflsa/')
        elif obsid in [60983]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2002/OCT2102/oflsc/')
        elif obsid in [60640, 60639, 60638, 60637, 60636, 60635, 60634, 60633]:
            raise ValueError("Starcheck not available for PCHECK_JUL2003")
        elif obsid in [60616, 60615]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2003/JUL2103/oflsc/')
        elif obsid in [3911]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2003/JUL2803/oflsc/')
        elif obsid in [4162]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2003/SEP2903/oflsa/')
        elif obsid in [60401]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2004/JAN1904/oflsb/')
        elif obsid in [59921, 5035]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2004/DEC1404/oflsc/')
        elif obsid in [58548, 58547, 58546, 7753]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2007/JAN2907/oflsb/')
        elif obsid in [7936, 7463]:
            starcheck = get_starcheck_catalog(obsid, mp_dir='/2007/MAY2807/oflsb/')
        else:
            raise ValueError("Problem looking up starcheck for {}".format(obsid))
    if starcheck is None or 'cat' not in starcheck or not len(starcheck['cat']):
        raise ValueError('No starcheck catalog found for {}'.format(manvr.get_obsid()))
    starcat_time = DateTime(starcheck['cat']['mp_starcat_time'][0]).secs
    starcat_dtime = starcat_time - DateTime(manvr.start).secs
    # If it looks like the wrong starcheck by time, give up
    if abs(starcat_dtime) > 300:
        raise ValueError("Starcheck cat time delta is {}".format(starcat_dtime))
    if abs(starcat_dtime) > 30:
        logger.warning(
            "Starcheck cat time delta of {} is > 30 sec".format(abs(starcat_dtime))
        )
    vals, times, one_shot, star_info = get_modern_data(manvr, dwell, starcheck)

    acq_stats = calc_acq_stats(manvr, vals, times)
    obsid_info = {
        'obsid': obsid,
        'obi': obspar['obi_num'],
        'acq_start': acq_start,
        'guide_start': guide_start,
        'guide_tstart': DateTime(guide_start).secs,
        'one_shot_length': one_shot,
        'revision': '1.0',
    }
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    # Filter the catalog to be just acquisition stars
    catalog = catalog[(catalog['type'] == 'ACQ') | (catalog['type'] == 'BOT')]
    time = DateTime(guide_start).secs
    ccd_temp = np.mean(fetch_sci.MSID('AACCCDPT', time - 250, time + 250).vals)
    warm_threshold = 100.0
    warm_frac = dark_model.get_warm_fracs(warm_threshold, time, ccd_temp)
    temps = {'ccd_temp': ccd_temp, 'n100_warm_frac': warm_frac}
    return obsid_info, acq_stats, star_info, catalog, temps


def table_acq_stats(obsid_info, acq_stats, star_info, catalog, temp):
    logger.info("arranging stats into tabular data")
    cols = (
        ACQ_COLS['obs']
        + ACQ_COLS['cat']
        + ACQ_COLS['stat']
        + ACQ_COLS['agasc']
        + ACQ_COLS['temp']
        + ACQ_COLS['bad']
    )
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
        row['ccd_temp'] = temp['ccd_temp']
        row['n100_warm_frac'] = temp['n100_warm_frac']
        row['known_bad'] = False
        row['bad_comment'] = ''
    # Exclude any rows that are missing
    table = table[~missing_slots]
    return table


def _save_acq_stats(t):
    if not os.path.exists(table_file):
        cols = (
            ACQ_COLS['obs']
            + ACQ_COLS['cat']
            + ACQ_COLS['stat']
            + ACQ_COLS['agasc']
            + ACQ_COLS['temp']
            + ACQ_COLS['bad']
        )
        desc, byteorder = tables.descr_from_dtype(np.dtype(cols))
        filters = tables.Filters(complevel=5, complib='zlib')
        h5 = tables.open_file(table_file, 'a')
        tbl = h5.create_table('/', 'data', desc, filters=filters, expectedrows=1e6)
        tbl.cols.obsid.create_index()
        tbl.cols.guide_tstart.create_csindex()
        tbl.cols.agasc_id.create_index()
        h5.close()
        del h5
    h5 = tables.open_file(table_file, 'a')
    tbl = h5.get_node('/', 'data')
    have_obsid_coord = tbl.get_where_list(
        '(obsid == {}) & (obi == {})'.format(t[0]['obsid'], t[0]['obi']), sort=True
    )
    if len(have_obsid_coord):
        obsid_rec = tbl.read_coordinates(have_obsid_coord)
        if len(obsid_rec) != len(t):
            raise ValueError(
                "Could not update {}; different number of slots".format(t[0]['obsid'])
            )
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
        obsids = _get_obsids_to_update(check_missing=opt.check_missing)
    for obsid in obsids:
        import time

        t = time.localtime()
        # Don't run during kadi update
        if t.tm_hour == 7:
            logger.info("Sleeping")
            time.sleep(3720)
        logger.info("Processing obsid {}".format(obsid))
        try:
            obsid_info, acq_stats, star_info, catalog, temp = calc_stats(obsid)
        except Exception as e:
            logger.info("Skipping obsid {}: {}".format(obsid, e))
            continue
        if not len(acq_stats):
            logger.info("Skipping obsid {}, no stats determined".format(obsid))
            continue

        t = table_acq_stats(obsid_info, acq_stats, star_info, catalog, temp)
        # Check for acquisition anomalies
        warn_on_acq_anom(t, opt.email)
        _save_acq_stats(t)


def main():
    opt = get_options()
    if opt.email is None:
        opt.email = DEFAULT_EMAIL
    update(opt)


if __name__ == '__main__':
    main()
