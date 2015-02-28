from kadi import events
import numpy as np
from Chandra.Time import DateTime
from Ska.engarchive import fetch
from astropy.table import Table, Column

msids = [
    'COBSRQID', 'CORADMEN', '3TSCPOS', '3TSCMOVE',
    'AOACASEQ', 'AOACQSUC', 'AOFREACQ', 'AOFWAIT', 'AOREPEAT',
    'AOACSTAT', 'AOACHIBK', 'AOFSTAR', 'AOFATTMD', 'AOACPRGS',
    'AOATUPST', 'AONSTARS', 'AOPCADMD']

per_slot = ['AOACQID', 'AOACFCT', 'AOIMAGE',
            'AOACMAG', 'AOACYAN', 'AOACZAN',
            'AOACICC', 'AOACIDP', 'AOACIIR', 'AOACIMS',
            'AOACIQB', 'AOACISP']

slot_msids = [field + '%s' % slot
              for field in per_slot
              for slot in range(0, 8)]


def get_telem_table(start_time, stop_time):

    pcad_data = fetch.MSIDset(msids + slot_msids,
                              DateTime(start_time).secs - 328,
                              DateTime(stop_time).secs + 328)
    range_times = pcad_data['AOACASEQ'].times[
        (pcad_data['AOACASEQ'].times >= DateTime(start_time).secs)
        & (pcad_data['AOACASEQ'].times <= DateTime(stop_time).secs)]
    pcad_data.interpolate(times=range_times)
    vals = Table([pcad_data[col].vals for col in msids],
                 names=msids)
    for field in slot_msids:
        vals.add_column(Column(name=field, data=pcad_data[field].vals))
        times = Table([pcad_data['AOACASEQ'].times], names=['time'])

    def compress_data(data, dtime):
        return data[data['AOREPEAT'] == '0 '], dtime[data['AOREPEAT'] == '0 ']

    vals, times = compress_data(vals, times)

    # make a list of dicts of the table
    simple_data = []
    for drow, trow in zip(vals, times):
        slot_data = {'slots': [],
                     'time': trow['time'],
                     'date': DateTime(trow['time']).date}
        for m in msids:
            if m == 'COBSRQID':
                slot_data[m] = int(drow[m])
            else:
                slot_data[m] = drow[m]
        for slot in range(0, 8):
            row_dict = {'slot': slot}
            for col in per_slot:
                row_dict[col] = drow['{}{}'.format(col, slot)]
            slot_data['slots'].append(row_dict)
        simple_data.append(slot_data)

    return simple_data

