from kadi import events
import numpy as np
from Chandra.Time import DateTime
from Ska.engarchive import fetch
from astropy.table import Table, Column
from mica.quaternion import Quat
import Ska.quatutil
import mica.starcheck
import agasc

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
    acqs = catalog[(catalog['type'] == 'BOT') | (catalog['type'] == 'ACQ')]

    # Compute the multiplicative factor to convert from the AGASC proper motion
    # field to degrees.  The AGASC PM is specified in milliarcsecs / year, so this
    # is dyear * (degrees / milliarcsec)
    agasc_equinox = DateTime('2000:001:00:00:00.000')
    dyear = (DateTime(times[0]) - agasc_equinox) / 365.25
    pm_to_degrees = dyear / (3600. * 1000.)
    R2A = 206264.81

    dy = {}
    dz = {}
    for slot in range(0, 8):
        agasc_id = acqs[acqs['slot'] == slot][0]['id']
        star = agasc.get_star(agasc_id)
        ra = star['RA']
        dec = star['DEC']
        if star['PM_RA'] != -9999:
            ra = star['RA'] + star['PM_RA'] * pm_to_degrees
        if star['PM_DEC'] != -9999:
            dec = star['DEC'] + star['PM_DEC'] * pm_to_degrees
        star_pos_eci = Ska.quatutil.radec2eci(ra, dec)
        d_aca = np.dot(np.dot(aca_misalign, Ts.transpose(0, 2, 1)),
                       star_pos_eci).transpose()
        yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
        zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
        dy[slot] = vals['AOACYAN{}'.format(slot)] - yag
        dz[slot] = vals['AOACZAN{}'.format(slot)] - zag

    return dy, dz

def get_acq_table(obsid):

    manvrs = events.manvrs.filter(obsid=obsid)
    if not len(manvrs):
        return None
    manvr = manvrs[0]

    start_time = DateTime(manvr.acq_start).secs
    stop_time = DateTime(manvr.guide_start).secs + 3
    acq_data = fetch.MSIDset(msids + slot_msids, start_time, stop_time)

    vals = Table([acq_data[col].vals for col in msids],
                 names=msids)
    for field in slot_msids:
        vals.add_column(Column(name=field, data=acq_data[field].vals))
        times = Table([acq_data['AOACASEQ'].times], names=['time'])

    def compress_data(data, dtime):
        return data[data['AOREPEAT'] == '0 '], dtime[data['AOREPEAT'] == '0 ']

    vals, times = compress_data(vals, times)
    d_times = times['time'] - DateTime(manvr.guide_start).secs

    # Get the catalog for the stars
    # This is used both to map ACQID to the right slot and
    # to get the star positions to estimate deltas later
    starcheck = mica.starcheck.get_starcheck_catalog(int(obsid))
    if 'cat' not in starcheck:
        raise ValueError('No starcheck catalog found for {}'.format(obsid))
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    # Filter the catalog to be just acquisition stars
    catalog = catalog[(catalog['type'] == 'ACQ') | (catalog['type'] == 'BOT')]
    slot_for_pos = [cat_row['slot'] for cat_row in catalog]
    pos_for_slot = dict([(slot, idx) for idx, slot in enumerate(slot_for_pos)])
    # Also, save out the starcheck index for each slot for later
    index_for_slot = dict([(cat_row['slot'], cat_row['idx'])
                           for cat_row in catalog])

    # Estimate the offsets from the expected catalog positions
    dy, dz = deltas_vs_obc_quat(vals, times['time'], catalog)
    for slot in range(0, 8):
        vals.add_column(Column(name='dy{}'.format(slot),
                               data=dy[slot].data))
        vals.add_column(Column(name='dz{}'.format(slot),
                               data=dz[slot].data))
        cat_entry = catalog[catalog['slot'] == slot][0]
        dmag = vals['AOACMAG{}'.format(slot)] - cat_entry['mag']
        vals.add_column(Column(name='dmag{}'.format(slot),
                               data=dmag.data))

    # make a list of dicts of the table
    simple_data = []
    for drow, trow in zip(vals, times):
        slot_data = {'slots': [],
                     'time': trow['time'],
                     'aorfstr1_slot': slot_for_pos[int(drow['AORFSTR1'])],
                     'aorfstr2_slot': slot_for_pos[int(drow['AORFSTR2'])],
                     }
        for m in msids:
            slot_data[m] = drow[m]
        for slot in range(0, 8):
            row_dict = {'slot': slot,
                        'catpos': pos_for_slot[slot],
                        'index': index_for_slot[slot]}
            for col in per_slot:
                if col not in ['AOACQID']:
                    row_dict[col] = drow['{}{}'.format(col, slot)]
            for col in ['dy', 'dz', 'dmag']:
                row_dict[col] = drow['{}{}'.format(col, slot)]
            row_dict['POS_ACQID'] = drow['AOACQID{}'.format(pos_for_slot[slot])]
            slot_data['slots'].append(row_dict)
        simple_data.append(slot_data)

    return simple_data





