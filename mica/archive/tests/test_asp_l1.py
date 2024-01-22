# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pathlib import Path
import tempfile

import numpy as np
import pytest
import ska_dbi
from astropy.table import Table
from kadi import events
from Quaternion import Quat, normalize
from Ska.engarchive import fetch
from testr import test_helper

from .. import asp_l1, obsid_archive

HAS_L1_ARCHIVE = Path(asp_l1.CONFIG["data_root"]).exists()


def compare_obc_and_asol(atts, times, recs, ptol=2, ytol=2, rtol=65):
    """
    Check that obc solution and ground solution have relatively close quats.
    Note that because the ground aspect solution is not continuous (only run over science obsids
    in aspect intervals) that this method only checks those intervals.

    :param atts: attitudes from get_atts
    :param times: times from get_atts
    :param recs: list of dicts of header records from get_atts
    :param ptol: dq quaternion pitch tolerance in arcsecs
    :param ytol: dq quaternion yaw tolerance in arcsecs
    :param rtol: dq quaternion roll tolerance in arcsecs
    """
    for ai in recs:
        telem = fetch.Msidset(["aoattqt*"], ai["TSTART"], ai["TSTOP"])
        obc_atts = np.vstack(
            [telem["aoattqt{}".format(idx)].vals for idx in [1, 2, 3, 4]]
        ).transpose()
        obc_times = telem["aoattqt1"].times
        # Test that, at the obc times, the onboard solution and the ground solution
        # are reasonably close
        idxs = np.searchsorted(times[:-1], obc_times)
        dps = []
        dys = []
        drs = []
        for att, obc_att in zip(atts[idxs], obc_atts):
            dq = Quat(normalize(att)).dq(Quat(normalize(obc_att)))
            dps.append(dq.pitch * 3600)
            dys.append(dq.yaw * 3600)
            drs.append(dq.roll0 * 3600)
        dps = np.array(dps)
        dys = np.array(dys)
        drs = np.array(drs)
        assert np.all(np.abs(dys) < ytol)
        assert np.all(np.abs(dps) < ptol)
        assert np.all(np.abs(drs) < rtol)


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
@pytest.mark.parametrize("obsid", [14333, 15175, 5438, 2121])
def test_get_atts_obsid(obsid):
    atts, times, recs = asp_l1.get_atts(obsid=obsid)
    compare_obc_and_asol(atts, times, recs)


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
def test_get_atts_time():
    start = "2014:001:00:00:00.000"
    stop = "2014:005:00:00:00.000"
    atts, times, recs = asp_l1.get_atts(start=start, stop=stop)
    assert len(atts) == len(times)
    compare_obc_and_asol(atts, times, recs)
    dwells = events.dwells.filter(start, stop)
    for dwell in dwells:
        if dwell.get_obsid() > 38000:
            continue
        # check that more than 90% of the kalman interval is in times fetched from get_atts
        ok = (times < dwell.tstop) & (times > dwell.tstart)
        assert (times[ok][-1] - times[ok][0]) > dwell.dur * 0.90
        # also assert that the number of ~.25sec samples works out
        assert (len(times[ok]) * .25625) > dwell.dur * .90


@pytest.mark.skipif("not test_helper.on_head_network()", reason="Not on HEAD network")
def test_update_l1_archive(tmp_path):

    config = asp_l1.CONFIG.copy()
    config["data_root"] = tmp_path / "asp1"
    config["temp_root"] = tmp_path / "temp"
    config["bad_obsids"] = tmp_path / "asp1" / "asp_l1_bad_obsids.dat"
    config["firstrun"] = True
    config["rebuild"] = True
    config["obsid"] = 1
    config["version"] = 4
    archive = obsid_archive.ObsArchive(config)
    obsids = archive.update()

    with ska_dbi.DBI(
        dbi="sqlite", server=config["data_root"] / "archfiles.db3"
    ) as db:
        dat = Table(db.fetchall("select * from archfiles"))
        dat.sort("filename")
        assert dat[
            "filename",
            "filetime",
            "ascdsver",
            "caldbver",
            "content",
            "revision",
            "obsid",
        ].pformat_all() == (
            [
                "             filename             filetime ascdsver caldbver  content   revision obsid",
                "--------------------------------- -------- -------- -------- ---------- -------- -----",
                " pcadf059904356N004_acal1.fits.gz 59904356  8.3.2.1    4.3.0     ACACAL        4     1",
                " pcadf059904356N004_acen1.fits.gz 59904356  8.3.2.1    4.3.0    ACACENT        4     1",
                "pcadf059904356N004_aqual1.fits.gz 59904356  8.3.2.1    4.3.0    ASPQUAL        4     1",
                " pcadf059904356N004_asol1.fits.gz 59904356  8.3.2.1    4.3.0     ASPSOL        4     1",
                " pcadf059904356N004_bpix1.fits.gz 59904356  8.3.2.1    4.3.0 ACA_BADPIX        4     1",
                "pcadf059904356N004_fidpr1.fits.gz 59904356  8.3.2.1    4.3.0   FIDPROPS        4     1",
                " pcadf059904356N004_gcal1.fits.gz 59904356  8.3.2.1    4.3.0    GYROCAL        4     1",
                " pcadf059904356N004_gspr1.fits.gz 59904356  8.3.2.1    4.3.0    GSPROPS        4     1",
            ]
        )


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
def test_get_atts_filter():
    # Obsid 19039 has a momentum dump that shows up in asp_sol_status
    atts, times, recs = asp_l1.get_atts(obsid=19039)
    uf_atts, uf_times, uf_recs = asp_l1.get_atts(obsid=19039, filter=False)
    # Confirm that approximately 212 seconds are filtered
    assert np.abs((len(uf_atts) - len(atts)) * 0.25625 - 212.2) < 5
