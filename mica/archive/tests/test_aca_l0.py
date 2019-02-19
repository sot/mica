# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import os
import pytest
from astropy.table import Table
import numpy as np
from mica.archive import aca_l0, asp_l1
from Ska.Numpy import interpolate

has_l0_2012_archive = os.path.exists(os.path.join(aca_l0.CONFIG['data_root'], '2012'))


@pytest.mark.skipif('not has_l0_2012_archive', reason='Test requires 2012 L0 archive')
def test_l0_images_meta():
    """
    Confirm meta values match reference/regress values
    """
    imgs = aca_l0.get_l0_images(467055635, 467055639, slot=7)
    assert imgs[0].meta == {'BGDAVG': 253,
                            'IMGCOL0': 7,
                            'IMGFUNC1': 2,
                            'IMGROW0': -12,
                            'IMGSIZE': 8,
                            'IMGSTAT': 0,
                            'IMGSCALE': 1025,
                            'INTEG': np.float32(1.696),
                            'TIME': np.float64(467055637.49031752)}

has_l0_2007_archive = os.path.exists(os.path.join(aca_l0.CONFIG['data_root'], '2007'))


@pytest.mark.skipif('not has_l0_2007_archive', reason='Test requires 2007 L0 archive')
def test_get_l0_images():
    """
    Do a validation test of get_l0_images:
    - Get 20 mins of image data for slot 6 of obsid 8008 (very nice clean stars)
    - Do first moment centroids in row and col
    - Compare to aspect pipeline FM centroids for same slot data

    This is a deep test that all the signs are right.  If not then everything
    breaks badly because the star image doesn't move in sync with row0, col0.
    """
    start = '2007:002:06:00:00'
    stop = '2007:002:06:20:00'

    imgs = aca_l0.get_l0_images(start, stop, slot=6)

    files = asp_l1.get_files(8008, content=['ACACENT'])
    acen = Table.read(files[0])
    # Pick FM centroids for slot 6
    ok = (acen['alg'] == 1) & (acen['slot'] == 6)
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

    rcen = interpolate(acen['cent_i'], acen['time'], times)
    ccen = interpolate(acen['cent_j'], acen['time'], times)

    assert np.all(np.abs(rcen - rcs) < 0.05)
    assert np.all(np.abs(ccen - ccs) < 0.05)
