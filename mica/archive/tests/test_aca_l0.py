# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division

import os

import numpy as np
import pytest
from astropy.table import Table
from chandra_aca import maude_decom
from cxotime import CxoTime
from Ska.Numpy import interpolate

from mica.archive import aca_l0, asp_l1

has_l0_2012_archive = os.path.exists(os.path.join(aca_l0.CONFIG["data_root"], "2012"))


@pytest.mark.skipif("not has_l0_2012_archive", reason="Test requires 2012 L0 archive")
def test_l0_images_meta():
    """
    Confirm meta values match reference/regress values
    """
    imgs = aca_l0.get_l0_images(467055635, 467055639, slot=7)
    assert imgs[0].meta == {
        "BGDAVG": 253,
        "IMGCOL0": 7,
        "IMGFUNC1": 2,
        "IMGROW0": -12,
        "IMGSIZE": 8,
        "IMGSTAT": 0,
        "IMGSCALE": 1025,
        "INTEG": np.float32(1.696),
        "TIME": np.float64(467055637.49031752),
    }


@pytest.mark.skipif(not has_l0_2012_archive, reason="Test requires 2012 L0 archive")
def test_get_aca_images():
    """
    Confirm mica and maude images sources agree for a small time range.
    """
    tstart = CxoTime("2012:180").secs
    imgs_table_maude = maude_decom.get_aca_images(tstart, tstart + 30)
    imgs_table_maude.sort(keys=["TIME", "IMGNUM"])
    imgs_table_mica = aca_l0.get_aca_images(tstart, tstart + 30)
    assert len(imgs_table_maude) == len(imgs_table_mica)
    for col in ["IMG", "IMGROW0_8X8", "IMGCOL0_8X8", "TIME"]:
        assert np.allclose(
            imgs_table_maude[col], imgs_table_mica[col], atol=0.1, rtol=0
        )
    for col in ["IMGNUM"]:
        assert np.all(imgs_table_maude[col] == imgs_table_maude[col])


@pytest.mark.skipif(not has_l0_2012_archive, reason="Test requires 2012 L0 archive")
def test_get_aca_images_empty():
    """
    Confirm that get_aca_images returns an empty table when no images are found
    """
    tstart = CxoTime("2012:180").secs
    # select on a non-sensical imgsize to ensure no images are found
    imgs_table_mica = aca_l0.get_aca_images(tstart, tstart + 30, imgsize=[7])
    assert len(imgs_table_mica) == 0


has_l0_2007_archive = os.path.exists(os.path.join(aca_l0.CONFIG["data_root"], "2007"))
has_asp_l1 = os.path.exists(os.path.join(asp_l1.CONFIG["data_root"]))


@pytest.mark.skipif(
    "not has_l0_2007_archive or not has_asp_l1",
    reason="Test requires 2007 L0 archive and L1 aspect archive",
)
def test_get_l0_images():
    """
    Do a validation test of get_l0_images:
    - Get 20 mins of image data for slot 6 of obsid 8008 (very nice clean stars)
    - Do first moment centroids in row and col
    - Compare to aspect pipeline FM centroids for same slot data

    This is a deep test that all the signs are right.  If not then everything
    breaks badly because the star image doesn't move in sync with row0, col0.
    """
    start = "2007:002:06:00:00"
    stop = "2007:002:06:20:00"

    imgs = aca_l0.get_l0_images(start, stop, slot=6)

    files = asp_l1.get_files(8008, content=["ACACENT"])
    acen = Table.read(files[0])
    # Pick FM centroids for slot 6
    ok = (acen["alg"] == 1) & (acen["slot"] == 6)
    acen = acen[ok]

    # Row and col centroids
    rcs = []
    ccs = []
    times = [img.TIME for img in imgs]

    # Easy way to do FM centroids with mgrid
    rw, cw = np.mgrid[0:6, 0:6]
    # rw = [[0, 0, 0, 0, 0, 0],
    #       [1, 1, 1, 1, 1, 1],
    #       [2, 2, 2, 2, 2, 2],
    #       [3, 3, 3, 3, 3, 3],
    #       [4, 4, 4, 4, 4, 4],
    #       [5, 5, 5, 5, 5, 5]]

    for img in imgs:
        norm = np.sum(img)
        rcs.append(np.sum(img * rw) / norm + img.row0)
        ccs.append(np.sum(img * cw) / norm + img.col0)

    rcen = interpolate(acen["cent_i"], acen["time"], times)
    ccen = interpolate(acen["cent_j"], acen["time"], times)

    assert np.all(np.abs(rcen - rcs) < 0.05)
    assert np.all(np.abs(ccen - ccs) < 0.05)


@pytest.mark.skipif(
    "not has_l0_2007_archive or not has_asp_l1",
    reason="Test requires 2007 L0 archive and L1 aspect archive",
)
def test_get_slot_data_8x8():
    """
    Do a validation test of get_l0_images:
    - Get 20 mins of image data for slot 6 of obsid 8008 (very nice clean stars)
    - Do first moment centroids in row and col
    - Compare to aspect pipeline FM centroids for same slot data

    This is a deep test that all the signs are right.  If not then everything
    breaks badly because the star image doesn't move in sync with row0, col0.
    """
    start = "2007:002:06:00:00"
    stop = "2007:002:06:20:00"

    slot_data = aca_l0.get_slot_data(start, stop, slot=6, centered_8x8=True)

    files = asp_l1.get_files(8008, content=["ACACENT"])
    acen = Table.read(files[0])
    # Pick FM centroids for slot 6
    ok = (acen["alg"] == 1) & (acen["slot"] == 6)
    acen = acen[ok]

    # Row and col centroids
    times = slot_data["TIME"]

    # Easy way to do FM centroids with mgrid
    rw, cw = np.mgrid[0:8, 0:8]

    img_raw = slot_data["IMGRAW"]  # np.round(slot_data['IMGRAW']).astype(int)
    norm = np.sum(img_raw, axis=(1, 2))
    rcs = np.sum(img_raw * rw, axis=(1, 2)) / norm + slot_data["IMGROW0"] - 1
    ccs = np.sum(img_raw * cw, axis=(1, 2)) / norm + slot_data["IMGCOL0"] - 1

    rcen = interpolate(acen["cent_i"], acen["time"], times)
    ccen = interpolate(acen["cent_j"], acen["time"], times)

    assert np.all(np.abs(rcen - rcs) < 0.05)
    assert np.all(np.abs(ccen - ccs) < 0.05)
