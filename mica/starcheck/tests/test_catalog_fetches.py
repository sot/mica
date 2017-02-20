import numpy as np
from Quaternion import normalize, Quat
from Chandra.Time import DateTime
from kadi import events
from Ska.engarchive import fetch
from Ska.quatutil import radec2yagzag, yagzag2radec

from .. import starcheck

def get_cmd_quat(date):
    date = DateTime(date)
    cmd_quats = fetch.MSIDset(['AOCMDQT{}'.format(i) for i in [1, 2, 3]], date.secs, date.secs + 120)
    cmd_q4 = np.sqrt(np.abs(1 - cmd_quats['AOCMDQT1'].vals[0]**2
                            - cmd_quats['AOCMDQT2'].vals[0]**2
                            - cmd_quats['AOCMDQT3'].vals[0]**2))
    return Quat(normalize([cmd_quats['AOCMDQT1'].vals[0],
                                      cmd_quats['AOCMDQT2'].vals[0],
                                      cmd_quats['AOCMDQT3'].vals[0],
                                      cmd_q4]))


def get_trak_cat_from_telem(date, cmd_quat):
    date = DateTime(date)
    msids = ["{}{}".format(m, i) for m in ['AOACYAN', 'AOACZAN', 'AOACFID', 'AOIMAGE', 'AOACFCT']
             for i in range(0, 8)]
    telem = fetch.MSIDset(['AOACASEQ', 'CORADMEN', 'AOPCADMD']
                           + msids, date.secs, date.secs + 300)
    att = fetch.MSIDset(['AOATTQT{}'.format(i) for i in [1, 2, 3, 4]], date.secs, date.secs + 300)
    cat = {}
    for slot in range(0, 8):
        track = telem['AOACFCT{}'.format(slot)].vals == 'TRAK'
        fid = telem['AOACFID{}'.format(slot)].vals == 'FID '
        n = 30
        if np.count_nonzero(track) < n:
            continue
        if np.any(fid & track):
            cat[slot] = {'type': 'FID',
                         'yag': telem['AOACYAN{}'.format(slot)].vals[fid & track][0],
                         'zag': telem['AOACZAN{}'.format(slot)].vals[fid & track][0],}
        else:
            yags = []
            zags = []
            for sample in range(0, n):
                qref = Quat(normalize([att['AOATTQT{}'.format(i)].vals[track][sample]
                                       for i in [1, 2, 3, 4]]))
                ra, dec = yagzag2radec(telem['AOACYAN{}'.format(slot)].vals[track][sample] / 3600.,
                                       telem['AOACZAN{}'.format(slot)].vals[track][sample] / 3600.,
                                       qref)
                yag, zag = radec2yagzag(ra, dec, cmd_quat)
                yags.append(yag)
                zags.append(zag)
            # This doesn't detect MON just yet
            cat[slot] = {'type': 'STAR',
                         'yag': np.median(yags) * 3600.,
                         'zag': np.median(zags) * 3600.,}
    return cat, telem


def test_validate_catalogs_over_range():
    start = '2017:001'
    stop = '2017:005'
    dwells = events.dwells.filter(start, stop)
    for dwell in dwells:
        telem_quat = get_cmd_quat(dwell.start)
        cat, telem = get_trak_cat_from_telem(dwell.start, telem_quat)
        sc = starcheck.get_starcheck_catalog_at_date(dwell.start)
        sc_quat = Quat([sc['manvr'][-1]["target_Q{}".format(i)] for i in [1, 2, 3, 4]])
        dq = sc_quat.dq(telem_quat)
        if ((np.abs(dq.pitch) * 3600 > 1) or (np.abs(dq.yaw) * 3600 > 1)):
            if dwell.manvr.template != 'unknown':
                raise ValueError("Unexpected offset in pointing")
            else:
                print(dwell.start, "pointing offset but has unknown manvr template")
                continue
        trak_sc = sc['cat'][sc['cat']['type'] != 'ACQ']
        for slot in range(0, 8):
            if slot in cat:
                slot_match = trak_sc[trak_sc['slot'] == slot]
                if not len(slot_match):
                    raise ValueError("missing slot in catalog")
                trak_sc_slot = slot_match[0]
                if trak_sc_slot['type'] == 'FID':
                    assert cat[slot]['type'] == 'FID'
                assert np.abs(cat[slot]['yag'] - trak_sc_slot['yang']) < 20.0
                assert np.abs(cat[slot]['zag'] - trak_sc_slot['zang']) < 20.0

