# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import pytest
import numpy as np
from .. import vv
from .. import process
from ... import common

HAS_L1_ARCHIVE = os.path.exists(process.asp_l1_arch.CONFIG["data_root"])
HAS_VV_ARCHIVE = (
    os.path.exists(vv.FILES["data_root"])
    & os.path.exists(vv.FILES["asp1_proc_table"])
    & os.path.exists(vv.FILES["h5_file"])
)
HAS_VV_TABLE = os.path.exists(vv.FILES["h5_file"])


@pytest.mark.skipif("not HAS_VV_ARCHIVE", reason="Test requires vv archive")
def test_get_vv_dir():
    obsdir = vv.get_vv_dir(16504, version=1)
    assert obsdir == os.path.abspath(
        os.path.join(common.MICA_ARCHIVE, "vv/16/16504_v01")
    )


@pytest.mark.skipif("not HAS_VV_ARCHIVE", reason="Test requires vv archive")
def test_get_vv_files():
    obsfiles = vv.get_vv_files(16504, version=1)
    assert sorted(obsfiles)[-1] == os.path.abspath(
        os.path.join(common.MICA_ARCHIVE, "vv/16/16504_v01/vv_report.pkl")
    )


@pytest.mark.skipif("not HAS_VV_TABLE", reason="Test requires vv h5 table")
def test_get_rms_data():
    data = vv.get_rms_data()
    dz_rms = data[
        (data["obsid"] == 16505) & (data["slot"] == 4) & (data["revision"] == 3)
    ]["dz_rms"][0]
    assert np.allclose(dz_rms, 0.047886185719034906)


@pytest.mark.skipif("not HAS_VV_ARCHIVE", reason="Test requires vv archive")
def test_get_vv():
    obs = vv.get_vv(16504, version=1)
    assert np.allclose(obs["slots"]["7"]["dz_rms"], 0.11610256063309182)


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
def test_run_vv():
    obi = process.get_arch_vv(2121, version=3)
    assert np.allclose(obi.info()["sim"]["max_d_dy"], 0.002197265625)


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
def test_run_vv_omitted_slot():
    # This test run on obsid with omitted slot is just testing for unhandled exceptions
    process.get_arch_vv(19991, version="last")


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
def test_run_vv_multi_interval():
    # This test run on obsid with multiple intervals is just testing for unhandled exceptions
    process.get_arch_vv(18980, version="last")


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
def test_run_vv_omitted_fid():
    process.get_arch_vv(18978, version="last")


@pytest.mark.skipif("not HAS_L1_ARCHIVE", reason="Test requires L1 archive")
def test_run_vv_7_track_slots():
    # Run on an obsid with only 7 slots *commanded* during Kalman
    process.get_arch_vv(19847, version="last")
