# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Basic functionality and regression tests for ACA dark cal module.
"""

import numpy as np
import pytest

from ..aca_dark import dark_cal
from chandra_aca.aca_image import ACAImage


def test_date_to_dark_id():
    assert dark_cal.date_to_dark_id('2011-01-15T12:00:00') == '2011015'


def test_dark_id_to_date():
    assert dark_cal.dark_id_to_date('2011015') == '2011:015'


def test_dark_temp_scale():
    scale = dark_cal.dark_temp_scale(-10., -14)
    assert np.allclose(scale, 0.70)

    scale = dark_cal.dark_temp_scale(-10., -14, scale_4c=2.0)
    assert scale == 0.5  # Should be an exact match


def test_get_dark_cal_id():
    assert dark_cal.get_dark_cal_id('2007:008', 'nearest') == '2007006'
    assert dark_cal.get_dark_cal_id('2007:008', 'before') == '2007006'
    assert dark_cal.get_dark_cal_id('2007:008', 'after') == '2007069'


@pytest.mark.parametrize('aca_image', [True, False])
def test_get_dark_cal_image(aca_image):
    image = dark_cal.get_dark_cal_image('2007:008', aca_image=aca_image)
    assert image.shape == (1024, 1024)
    if aca_image:
        assert image.row0 == -512
        assert image.col0 == -512
        assert image.aca[-511, -511] == image[1, 1]
        assert image.aca[511, 0] == image[1023, 512]
    else:
        assert type(image) is np.ndarray


@pytest.mark.parametrize('aca_image', [True, False])
def test_get_dark_cal_props(aca_image):
    props = dark_cal.get_dark_cal_props('2007:008')
    assert len(props['replicas']) == 5
    assert props['start'] == '2007:006:01:56:46.817'

    props = dark_cal.get_dark_cal_props('2007:008', include_image=True, aca_image=aca_image)
    assert len(props['replicas']) == 5
    assert props['start'] == '2007:006:01:56:46.817'
    assert props['image'].shape == (1024, 1024)
    assert type(props['image']) is (ACAImage if aca_image else np.ndarray)


def test_get_dark_cal_props_table():
    props = dark_cal.get_dark_cal_props_table('2007:001', '2008:001')
    assert np.allclose(props['eb'], [24.6, 25.89, 51.13, 1.9])
