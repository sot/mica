# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import numpy as np
import pytest
from Chandra.Time import DateTime
from kadi import events
from kadi.commands import get_starcats as kadi_get_starcats
from Quaternion import Quat, normalize
from Ska.engarchive import fetch
from Ska.quatutil import radec2yagzag, yagzag2radec

from .. import starcheck
from mica.utils import load_name_to_mp_dir

HAS_SC_ARCHIVE = os.path.exists(starcheck.FILES["data_root"])


def get_cmd_quat(date):
    date = DateTime(date)
    cmd_quats = fetch.MSIDset(
        ["AOCMDQT{}".format(i) for i in [1, 2, 3]], date.secs, date.secs + 120
    )
    cmd_q4 = np.sqrt(
        np.abs(
            1
            - cmd_quats["AOCMDQT1"].vals[0] ** 2
            - cmd_quats["AOCMDQT2"].vals[0] ** 2
            - cmd_quats["AOCMDQT3"].vals[0] ** 2
        )
    )
    return Quat(
        normalize(
            [
                cmd_quats["AOCMDQT1"].vals[0],
                cmd_quats["AOCMDQT2"].vals[0],
                cmd_quats["AOCMDQT3"].vals[0],
                cmd_q4,
            ]
        )
    )


def get_trak_cat_from_telem(start, stop, cmd_quat):
    start = DateTime(start)
    stop = DateTime(stop)
    msids = [
        "{}{}".format(m, i)
        for m in ["AOACYAN", "AOACZAN", "AOACFID", "AOIMAGE", "AOACFCT"]
        for i in range(0, 8)
    ]
    telem = fetch.MSIDset(
        ["AOACASEQ", "CORADMEN", "AOPCADMD", "AONSTARS", "AOKALSTR"] + msids,
        start,
        stop,
    )
    att = fetch.MSIDset(["AOATTQT{}".format(i) for i in [1, 2, 3, 4]], start, stop)
    cat = {}
    for slot in range(0, 8):
        track = telem["AOACFCT{}".format(slot)].vals == "TRAK"
        fid = telem["AOACFID{}".format(slot)].vals == "FID "
        star = telem["AOIMAGE{}".format(slot)].vals == "STAR"
        n = 30
        if np.count_nonzero(track) < n:
            continue
        if np.any(fid & track):
            cat[slot] = {
                "type": "FID",
                "yag": telem["AOACYAN{}".format(slot)].vals[fid & track][0],
                "zag": telem["AOACZAN{}".format(slot)].vals[fid & track][0],
            }
        else:
            n_samples = np.count_nonzero(track & star)
            if n_samples < (n + 4):
                continue
            # If there is tracked data with a star, let's try to get our n samples from about
            # the middle of the range
            mid_point = int(n_samples / 2.0)
            yags = []
            zags = []
            for sample in range(mid_point - int(n / 2.0), mid_point + int(n / 2.0)):
                qref = Quat(
                    normalize(
                        [
                            att["AOATTQT{}".format(i)].vals[track & star][sample]
                            for i in [1, 2, 3, 4]
                        ]
                    )
                )
                ra, dec = yagzag2radec(
                    telem["AOACYAN{}".format(slot)].vals[track & star][sample] / 3600.0,
                    telem["AOACZAN{}".format(slot)].vals[track & star][sample] / 3600.0,
                    qref,
                )
                yag, zag = radec2yagzag(ra, dec, cmd_quat)
                yags.append(yag)
                zags.append(zag)
            # This doesn't detect MON just yet
            cat[slot] = {
                "type": "STAR",
                "yag": np.median(yags) * 3600.0,
                "zag": np.median(zags) * 3600.0,
            }
    return cat, telem


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_validate_catalogs_over_range():
    start = "2017:001:12:00:00"
    stop = "2017:002:12:00:00"
    dwells = events.dwells.filter(start, stop)
    for dwell in dwells:
        print(dwell)
        try:
            telem_quat = get_cmd_quat(dwell.start)
            # try to get the tracked telemetry for 1ks at the beginning of the dwell,
            # or if the dwell is shorter than that, just get the dwell
            cat, telem = get_trak_cat_from_telem(
                dwell.start, np.min([dwell.tstart + 100, dwell.tstop]), telem_quat
            )
        except ValueError as err:
            if "MSID" in str(err):
                pytest.skip("Eng archive MSID missing {}".format(err))

        sc = starcheck.get_starcheck_catalog_at_date(dwell.start)
        sc_quat = Quat([sc["manvr"][-1]["target_Q{}".format(i)] for i in [1, 2, 3, 4]])
        dq = sc_quat.dq(telem_quat)
        if (np.abs(dq.pitch) * 3600 > 1) or (np.abs(dq.yaw) * 3600 > 1):
            if dwell.manvr.template != "unknown":
                raise ValueError("Unexpected offset in pointing")
            else:
                print(dwell.start, "pointing offset but has unknown manvr template")
                continue
        trak_sc = sc["cat"][sc["cat"]["type"] != "ACQ"]
        for slot in range(0, 8):
            if slot in cat:
                slot_match = trak_sc[trak_sc["slot"] == slot]
                # the if statement may see false positives if fewer than a total of 8 BOT/GUI/MON
                # slots were commanded
                if not len(slot_match):
                    raise ValueError("missing slot in catalog")
                trak_sc_slot = slot_match[0]
                if trak_sc_slot["type"] == "FID":
                    assert cat[slot]["type"] == "FID"
                    assert np.abs(cat[slot]["yag"] - trak_sc_slot["yang"]) < 20.0
                    assert np.abs(cat[slot]["zag"] - trak_sc_slot["zang"]) < 20.0
                else:
                    assert np.abs(cat[slot]["yag"] - trak_sc_slot["yang"]) < 5.0
                    assert np.abs(cat[slot]["zag"] - trak_sc_slot["zang"]) < 5.0


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
@pytest.mark.parametrize(
    "test_case",
    [
        {"obsid": 19990, "mp_dir": "/2017/FEB2017/oflsa/", "n_cat_entries": 11},
        {"obsid": 17210, "mp_dir": "/2016/JAN2516/oflsa/", "n_cat_entries": 11},
        # 45312 is a gyro hold with no star catalog
        {"obsid": 45312, "mp_dir": "/2022/AUG2322/oflsa/", "n_cat_entries": 0},
    ],
)
def test_obsid_catalog_fetch(test_case):
    sc = starcheck.get_starcheck_catalog(test_case["obsid"])
    assert test_case["mp_dir"] == sc["mp_dir"]
    assert len(sc["cat"]) == test_case["n_cat_entries"]


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_obsid_catalog_fetch_dark_cal():
    # review a dark current replica obsid
    dcdir, dcstatus, dcdate = starcheck.get_mp_dir(49961)
    assert dcdir == "/2017/JUL0317/oflsb/"
    assert dcstatus == "no starcat"
    assert dcdate is None


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_monitor_fetch():
    mons = starcheck.get_monitor_windows(
        start="2009:002:12:00:00", stop="2009:007:12:00:00"
    )
    assert len(mons) == 10


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_methods():
    """
    Check that the get_starcat, get_dither, and get_att Spacecraft methods
    return reasonable values.
    """
    # Clear the cache for testing
    starcheck.OBS_CACHE.clear()

    # Get the catalog for any obsid to add something to the cache
    starcheck.get_starcat(2121)
    assert len(starcheck.OBS_CACHE) == 1

    # Get an check the values for the utility methods for obsid 19372
    obsid = 19372
    cat = starcheck.get_starcat(obsid)
    assert len(starcheck.OBS_CACHE) == 2

    # IDX and ID of the entries of this obsid catalog
    regress = {
        1: 1,
        2: 5,
        3: 6,
        4: 454035528,
        5: 454035856,
        6: 454428648,
        7: 454432416,
        8: 454433448,
        9: 454036528,
        10: 454039632,
        11: 454036472,
        12: 454038056,
    }
    assert len(regress) == len(cat)
    for row in cat:
        assert regress[row["idx"]] == row["id"]

    dither = starcheck.get_dither(obsid)
    assert dither == {
        "pitch_ampl": 8.0,
        "pitch_period": 707.1,
        "yaw_ampl": 8.0,
        "yaw_period": 1000.0,
    }

    att = starcheck.get_att(obsid)
    assert att == [209.04218, 47.227524, 357.020117]

    # Check a dark cal obs from SEP3019A and an old (pre-kadi) obsid
    for obsid in (47846, 2000):
        dither = starcheck.get_dither(obsid)
        assert dither == {
            "pitch_ampl": 0.0,
            "pitch_period": -999.0,
            "yaw_ampl": 0.0,
            "yaw_period": -999.0,
        }


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_vs_kadi():
    # Make sure this works for an obsid in kadi commands
    cat_mica = starcheck.get_starcat(8008)
    cat_kadi = kadi_get_starcats(obsid=8008, scenario="flight")[0]
    assert np.all(cat_mica["id"] == cat_kadi["id"])


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_db_obsid():
    obsid = starcheck.get_starcheck_db_obsid(
        "2022:017:05:15:06.000", mp_dir="/2022/JAN1722/oflsa/"
    )
    assert obsid == 45774

    # Science observation not run due to SCS-107. The as-run obsid is 45091 (obsid at
    # the time of the interrupt) but the planned obsid in the starcheck database is
    # 26269. The case above may be similar but didn't get documented.
    obsid = starcheck.get_starcheck_db_obsid(
        "2022:311:02:39:03.454", "/2022/NOV0722/oflsa/"
    )
    assert obsid == 25316


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_db_obsid_fail():
    with pytest.raises(ValueError, match="no starcheck entry"):
        starcheck.get_starcheck_db_obsid(
            "2022:017:05:15:06.000", mp_dir="/2022/JAN1722/oflsZZZ/"
        )


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_scs107():
    cat = starcheck.get_starcheck_catalog_at_date("2022:017:06:00:00.000")
    assert cat["obs"]["obsid"] == 45774
    cat_mica = cat["cat"]
    cat_kadi = kadi_get_starcats(
        obsid=26269, starcat_date="2022:017:05:15:06.000", scenario="flight"
    )[0]
    assert np.all(cat_mica["id"] == cat_kadi["id"])


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_with_mp_dir():
    """Obsid 21082 was scheduled in APR2318A and B, but not C, which was actually
    run. Use this to test the functionality of explicitly specifying load name."""

    starcheck.OBS_CACHE.clear()

    sc = starcheck.get_starcheck_catalog(21082, "/2018/APR2318/oflsa/")
    assert len(sc["cat"]) == 12
    assert (21082, "/2018/APR2318/oflsa/") in starcheck.OBS_CACHE

    sc = starcheck.get_starcheck_catalog(21082, "APR2318A")
    assert len(sc["cat"]) == 12
    assert (21082, "APR2318A") in starcheck.OBS_CACHE


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_for_no_star_catalog():
    """Test getting starcheck by date for a time that has no star catalog (gyro hold)"""
    cat = starcheck.get_starcheck_catalog_at_date("2023:099:04:21:40.719")
    exp = exp = {
        "mp_dir": "/2023/APR0323/oflsa/",
        "status": "ran",
        "obs": {
            "sc_id": 2477,
            "obsid": 44653,
            "obs_idx": 45,
            "point_ra": 304.0,
            "point_dec": -52.0,
            "point_roll": 107.709085,
            "target_id": None,
            "sci_instr": None,
            "sim_z_offset_steps": None,
            "sim_z_offset_mm": None,
            "grating": None,
            "dither_state": None,
            "dither_y_amp": None,
            "dither_y_period": None,
            "dither_z_amp": None,
            "dither_z_period": None,
            "mp_starcat_time": None,
            "mp_starcat_vcdu_cnt": None,
            "obsid_as_run": 44653,
        },
        "manvr": [],
        "warnings": [],
        "cat": [],
    }
    assert cat == exp


@pytest.mark.skipif("not HAS_SC_ARCHIVE", reason="Test requires starcheck archive")
def test_get_starcheck_for_no_starcheck_entry():
    """Test getting starcheck by date after a while in safe mode with no loads"""
    date = "2023:050:00:30:00"
    cat = starcheck.get_starcheck_catalog_at_date(date)
    assert cat is None


def test_load_name_to_mp_dir():
    mp_dir = load_name_to_mp_dir("DEC2506C")
    assert mp_dir == "/2006/DEC2506/oflsc/"
