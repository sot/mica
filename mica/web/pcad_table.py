from kadi import events
import numpy as np
from Chandra.Time import DateTime
from Ska.engarchive import fetch
from astropy.table import Table, Column
from mica.quaternion import Quat
import Ska.quatutil
import mica.starcheck
import agasc

R2A = 206264.81


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

def deltas_vs_obc_quat(vals, times, catalog):
    # Ignore misalign
    aca_misalign = np.array([[1.0,0,0], [0,1,0],[0,0,1]])
    q_att = Quat(np.array([vals['AOATTQT1'],
                           vals['AOATTQT2'],
                           vals['AOATTQT3'],
                           vals['AOATTQT4']]).transpose())
    Ts = q_att.transform

    dy = {}
    dz = {}
    for slot in range(0, 8):
        if np.any(catalog['slot'] == slot):
            agasc_id = catalog[catalog['slot'] == slot][0]['id']
        else:
            continue
        if agasc_id < 20:
            continue
        star = agasc.get_star(agasc_id, date=times[0])
        ra = star['RA_PMCORR']
        dec = star['DEC_PMCORR']
        star_pos_eci = Ska.quatutil.radec2eci(ra, dec)
        d_aca = np.dot(np.dot(aca_misalign, Ts.transpose(0, 2, 1)), star_pos_eci).transpose()
        yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
        zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
        dy[slot] = vals['AOACYAN{}'.format(slot)] - yag
        dz[slot] = vals['AOACZAN{}'.format(slot)] - zag

    return dy, dz

def get_acq_table(obsid):

    obsid = int(obsid)

    acq_start = None

    try:
        manvrs = events.manvrs.filter(obsid=obsid)
        acq_start = manvrs[0].acq_start
    except:
        try:
            sc = mica.starcheck.get_starcheck_catalog(obsid)
            acq_start = sc['manvr'][-1]['end_date']
            fetch.data_source.set('maude')
        except:
            raise

    start_time = DateTime(acq_start).secs
    stop_time = start_time + (60 * 8)
    print(start_time)
    print(stop_time)
    acq_data = fetch.MSIDset(msids + slot_msids, start_time, stop_time)
    acq_data.interpolate(1.025)
    vals = Table([acq_data[col].vals for col in msids + slot_msids],
                 names=msids + slot_msids)
    times = {'time': acq_data['AOACASEQ'].times}

    def compress_data(data, dtime):
        non_repeat = (data['AOREPEAT'] == '0 ') | (data['AOREPEAT'] == '0')
        return data[non_repeat], dtime[non_repeat]

    vals, times['time'] = compress_data(vals, times['time'])

    # Get the catalog for the stars
    # This is used both to map ACQID to the right slot and
    # to get the star positions to estimate deltas later
    timeline_at_acq = mica.starcheck.starcheck.get_timeline_at_date(acq_start)
    mp_dir = None if (timeline_at_acq is None) else timeline_at_acq['mp_dir']
    starcheck = mica.starcheck.get_starcheck_catalog(int(obsid), mp_dir=mp_dir)
    if 'cat' not in starcheck:
        raise ValueError('No starcheck catalog found for {}'.format(obsid))
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    # Filter the catalog to be just acquisition stars
    acq_catalog = catalog[(catalog['type'] == 'ACQ') | (catalog['type'] == 'BOT')]
    slot_for_pos = [cat_row['slot'] for cat_row in acq_catalog]
    pos_for_slot = dict([(slot, idx) for idx, slot in enumerate(slot_for_pos)])
    # Also, save out the starcheck index for each slot for later
    index_for_slot = dict([(cat_row['slot'], cat_row['idx']) for cat_row in acq_catalog])

    # Estimate the offsets from the expected catalog positions
    dy, dz = deltas_vs_obc_quat(vals, times['time'], acq_catalog)
    for slot in range(0, 8):
        if slot in dy:
            vals['acq_dy{}'.format(slot)] = dy[slot].data
            vals['acq_dz{}'.format(slot)] = dz[slot].data
            cat_entry = catalog[catalog['slot'] == slot][0]
            dmag = vals['AOACMAG{}'.format(slot)] - cat_entry['mag']
            vals['acq_dmag{}'.format(slot)] = dmag.data

    # Filter the catalog to be just tracked things
    gui_catalog = catalog[(catalog['type'] != 'ACQ')]
    # Estimate the offsets from the expected guide catalog positions
    gui_dy, gui_dz = deltas_vs_obc_quat(vals, times['time'], gui_catalog)
    for slot in range(0, 8):
        if slot in gui_dy:
            vals['gui_dy{}'.format(slot)] = gui_dy[slot].data
            vals['gui_dz{}'.format(slot)] = gui_dz[slot].data
            cat_entry = catalog[catalog['slot'] == slot][0]
            dmag = vals['AOACMAG{}'.format(slot)] - cat_entry['mag']
            vals['gui_dmag{}'.format(slot)] = dmag.data

    # make a list of dicts of the table
    simple_data = []
    kalm_start = None
    found_acq = None
    for drow, trow in zip(vals, times['time']):
#        if (drow['AOACASEQ'] == 'AQXN'):
#            found_acq = True
#        if (found_acq is not None) and (kalm_start is None) and (drow['AOACASEQ'] == 'KALM'):
#            kalm_start = trow
#        if (kalm_start is not None) and (trow > kalm_start + 5):
#            continue
        slot_data = {'slots': [],
                     'time': trow,
                     'aorfstr1_slot': slot_for_pos[int(drow['AORFSTR1'])],
                     'aorfstr2_slot': slot_for_pos[int(drow['AORFSTR2'])],
                     }
        for m in msids:
            slot_data[m] = drow[m]
        for slot in range(0, 8):
            row_dict = {'slot': slot,
                        'catpos': pos_for_slot[slot],
                        'index': index_for_slot[slot]}
            for col in per_slot + ['acq_dy', 'acq_dz', 'acq_dmag', 'gui_dy', 'gui_dz', 'gui_dmag']:
                if col not in ['AOACQID']:
                    col_key = '{}{}'.format(col, slot)
                    if col_key in drow.colnames:
                        row_dict[col] = drow[col_key]
            prefix = 'acq_' if drow['AOACASEQ'] == 'AQXN' else 'gui_'
            for col in ['dy', 'dz', 'dmag']:
                res_key = '{}{}{}'.format(prefix, col, slot)
                if res_key in drow.colnames:
                    row_dict[col] = drow[res_key]
            row_dict['POS_ACQID'] = drow['AOACQID{}'.format(pos_for_slot[slot])]
            slot_data['slots'].append(row_dict)
        simple_data.append(slot_data)

    return simple_data





