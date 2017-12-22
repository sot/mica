# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import numpy as np
from Quaternion import normalize, Quat
from Chandra.Time import DateTime
from kadi import events
from Ska.engarchive import fetch
from Ska.quatutil import radec2yagzag, yagzag2radec
import pytest

from .. import starcheck

HAS_SC_ARCHIVE = os.path.exists(starcheck.FILES['data_root'])

def get_cmd_quat(date):
    date = DateTime(date)
    cmd_quats = fetch.MSIDset(['AOCMDQT{}'.format(i) for i in [1, 2, 3]],
                              date.secs, date.secs + 120)
    cmd_q4 = np.sqrt(np.abs(1 - cmd_quats['AOCMDQT1'].vals[0]**2
                            - cmd_quats['AOCMDQT2'].vals[0]**2
                            - cmd_quats['AOCMDQT3'].vals[0]**2))
    return Quat(normalize([cmd_quats['AOCMDQT1'].vals[0],
                           cmd_quats['AOCMDQT2'].vals[0],
                           cmd_quats['AOCMDQT3'].vals[0],
                           cmd_q4]))


def get_trak_cat_from_telem(start, stop, cmd_quat):
    start = DateTime(start)
    stop = DateTime(stop)
    msids = ["{}{}".format(m, i) for m in ['AOACYAN', 'AOACZAN', 'AOACFID', 'AOIMAGE', 'AOACFCT']
             for i in range(0, 8)]
    telem = fetch.MSIDset(['AOACASEQ', 'CORADMEN', 'AOPCADMD', 'AONSTARS', 'AOKALSTR']
                          + msids, start, stop)
    att = fetch.MSIDset(['AOATTQT{}'.format(i) for i in [1, 2, 3, 4]], start, stop)
    cat = {}
    for slot in range(0, 8):
        track = telem['AOACFCT{}'.format(slot)].vals == 'TRAK'
        fid = telem['AOACFID{}'.format(slot)].vals == 'FID '
        star = telem['AOIMAGE{}'.format(slot)].vals == 'STAR'
        n = 30
        if np.count_nonzero(track) < n:
            continue
        if np.any(fid & track):
            cat[slot] = {'type': 'FID',
                         'yag': telem['AOACYAN{}'.format(slot)].vals[fid & track][0],
                         'zag': telem['AOACZAN{}'.format(slot)].vals[fid & track][0]}
        else:
            n_samples = np.count_nonzero(track & star)
            if n_samples < (n + 4):
                continue
            # If there is tracked data with a star, let's try to get our n samples from about
            # the middle of the range
            mid_point = int(n_samples / 2.)
            yags = []
            zags = []
            for sample in range(mid_point - int(n / 2.), mid_point + int(n / 2.)):
                qref = Quat(normalize([att['AOATTQT{}'.format(i)].vals[track & star][sample]
                                       for i in [1, 2, 3, 4]]))
                ra, dec = yagzag2radec(
                    telem['AOACYAN{}'.format(slot)].vals[track & star][sample] / 3600.,
                    telem['AOACZAN{}'.format(slot)].vals[track & star][sample] / 3600.,
                    qref)
                yag, zag = radec2yagzag(ra, dec, cmd_quat)
                yags.append(yag)
                zags.append(zag)
            # This doesn't detect MON just yet
            cat[slot] = {'type': 'STAR',
                         'yag': np.median(yags) * 3600.,
                         'zag': np.median(zags) * 3600.}
    return cat, telem


@pytest.mark.skipif('not HAS_SC_ARCHIVE', reason='Test requires starcheck archive')
def test_validate_catalogs_over_range():
    start = '2017:001'
    stop = '2017:002'
    dwells = events.dwells.filter(start, stop)
    for dwell in dwells:
        print(dwell)
        telem_quat = get_cmd_quat(dwell.start)
        # try to get the tracked telemetry for 1ks at the beginning of the dwell,
        # or if the dwell is shorter than that, just get the dwell
        cat, telem = get_trak_cat_from_telem(dwell.start,
                                             np.min([dwell.tstart + 100, dwell.tstop]),
                                             telem_quat)
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
                # the if statement may see false positives if fewer than a total of 8 BOT/GUI/MON
                # slots were commanded
                if not len(slot_match):
                    raise ValueError("missing slot in catalog")
                trak_sc_slot = slot_match[0]
                if trak_sc_slot['type'] == 'FID':
                    assert cat[slot]['type'] == 'FID'
                    assert np.abs(cat[slot]['yag'] - trak_sc_slot['yang']) < 20.0
                    assert np.abs(cat[slot]['zag'] - trak_sc_slot['zang']) < 20.0
                else:
                    assert np.abs(cat[slot]['yag'] - trak_sc_slot['yang']) < 5.0
                    assert np.abs(cat[slot]['zag'] - trak_sc_slot['zang']) < 5.0


@pytest.mark.skipif('not HAS_SC_ARCHIVE', reason='Test requires starcheck archive')
def test_obsid_catalog_fetch():
    tests = [{'obsid': 19990,
              'mp_dir': '/2017/FEB2017/oflsa/',
              'n_cat_entries': 11},
             {'obsid': 17210,
              'mp_dir': '/2016/JAN2516/oflsa/',
              'n_cat_entries': 11},
             {'obsid': 62668}]
    # 62668 should be an undercover with no catalog
    for t in tests:
        sc = starcheck.get_starcheck_catalog(t['obsid'])
        if 'mp_dir' in t:
            assert t['mp_dir'] == sc['mp_dir']
        if 'n_cat_entries' in t:
            assert len(sc['cat']) == t['n_cat_entries']
        if t['obsid'] == 62668:
            assert sc is None
    # review a dark current replica obsid
    dcdir, dcstatus, dcdate = starcheck.get_mp_dir(49961)
    assert dcdir == '/2017/JUL0317/oflsb/'
    assert dcstatus == 'no starcat'
    assert dcdate is None

@pytest.mark.skipif('not HAS_SC_ARCHIVE', reason='Test requires starcheck archive')
def test_monitor_fetch():
    mons = starcheck.get_monitor_windows(start='2009:002', stop='2009:007')
    assert len(mons) == 10
