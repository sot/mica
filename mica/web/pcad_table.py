from kadi import events
#import numpy as np
from Chandra.Time import DateTime
from Ska.engarchive import fetch
from astropy.table import Table, Column
#import mica.archive.aca_l0
import mica.starcheck

msids = ['AOACASEQ', 'AOACQSUC', 'AOFREACQ', 'AOFWAIT', 'AOREPEAT',
         'AOACSTAT', 'AOACHIBK', 'AOFSTAR', 'AOFATTMD', 'AOACPRGS',
         'AOATUPST', 'AONSTARS', 'AOPCADMD', 'AORFSTR1', 'AORFSTR2']
per_slot = ['AOACQID', 'AOACFCT', 'AOIMAGE',
            'AOACMAG', 'AOACYAN', 'AOACZAN',
            'AOACICC', 'AOACIDP', 'AOACIIR', 'AOACIMS',
            'AOACIQB', 'AOACISP']

slot_msids = [field + '%s' % slot
              for field in per_slot
              for slot in range(0, 8)]


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

    starcheck = mica.starcheck.get_starcheck_catalog(int(obsid))
    if 'cat' not in starcheck:
        raise ValueError('No starcheck catalog found for {}'.format(obsid))
    catalog = Table(starcheck['cat'])
    catalog.sort('idx')
    slot_for_pos = [cat_row['slot'] for cat_row in catalog if
                    (cat_row['type'] == 'ACQ') or (cat_row['type'] == 'BOT')]
    pos_for_slot = dict([(slot, idx) for idx, slot in enumerate(slot_for_pos)])


    #l0_data = {}
    #for slot in range(0, 8):
    #    l0_data[slot] = mica.archive.aca_l0.get_slot_data(start_time - 5,
    #                                                      stop_time,
    #                                                      slot)

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
        #for slot in [1]:
            row_dict = {'slot': slot, 'catpos': pos_for_slot[slot]}
            for col in per_slot:
                if col not in ['AOACQID']:
                    row_dict[col] = drow['{}{}'.format(col, slot)]
            row_dict['POS_ACQID'] = drow['AOACQID{}'.format(pos_for_slot[slot])]
            #idx = np.searchsorted(l0_data[slot]['TIME'], trow['time'])
            #for col in ['GLBSTAT', 'IMGSIZE', 'IMGSTAT', 'BGDAVG', 'BGDSTAT', 'IMGFUNC1']:
            #    row_dict[col] = l0_data[slot][idx - 1][col]
            slot_data['slots'].append(row_dict)
        simple_data.append(slot_data)

    return simple_data

#from django.template import Template, Context
#acq_template = open('templates/acq.html').read()
#template = Template(acq_template)
#c = Context({'vals': simple_data})
#f = open('out.html', 'w')
#f.write(template.render(c))
#f.close()




