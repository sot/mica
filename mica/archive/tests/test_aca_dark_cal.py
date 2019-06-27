# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Basic functionality and regression tests for ACA dark cal module.
"""
import os
import numpy as np
import pytest
import six

from ..aca_dark import dark_cal
from chandra_aca.aca_image import ACAImage

HAS_DARK_ARCHIVE = os.path.exists(dark_cal.MICA_FILES['dark_cals_dir'].abs)
IS_PY3 = six.PY3


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
def test_date_to_dark_id():
    assert dark_cal.date_to_dark_id('2011-01-15T12:00:00') == '2011015'


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
def test_dark_id_to_date():
    assert dark_cal.dark_id_to_date('2011015') == '2011:015'


def test_dark_temp_scale():
    scale = dark_cal.dark_temp_scale(-10., -14)
    assert np.allclose(scale, 0.70)

    scale = dark_cal.dark_temp_scale(-10., -14, scale_4c=2.0)
    assert scale == 0.5  # Should be an exact match


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
def test_get_dark_cal_id():
    assert dark_cal.get_dark_cal_id('2007:008', 'nearest') == '2007006'
    assert dark_cal.get_dark_cal_id('2007:008', 'before') == '2007006'
    assert dark_cal.get_dark_cal_id('2007:008', 'after') == '2007069'


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
@pytest.mark.parametrize('allow_negative', [True, False])
@pytest.mark.parametrize('aca_image', [True, False])
def test_get_dark_cal_image(aca_image, allow_negative):
    image = dark_cal.get_dark_cal_image('2007:008', aca_image=aca_image,
                                        allow_negative=allow_negative)
    assert image.shape == (1024, 1024)
    if aca_image:
        assert image.row0 == -512
        assert image.col0 == -512
        assert image.aca[-511, -511] == image[1, 1]
        assert image.aca[511, 0] == image[1023, 512]
    else:
        assert type(image) is np.ndarray

    # Raw dark cal images always have negative values unless clipped
    assert np.any(image < 0) == allow_negative


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
@pytest.mark.parametrize('allow_negative', [True, False])
@pytest.mark.parametrize('aca_image', [True, False])
def test_get_dark_cal_props(aca_image, allow_negative):
    props = dark_cal.get_dark_cal_props('2007:008')
    assert len(props['replicas']) == 5
    assert props['start'] == '2007:006:01:56:46.817'

    props = dark_cal.get_dark_cal_props('2007:008', include_image=True,
                                        aca_image=aca_image,
                                        allow_negative=allow_negative)
    assert len(props['replicas']) == 5
    assert props['start'] == '2007:006:01:56:46.817'
    assert props['image'].shape == (1024, 1024)
    assert np.any(props['image'] < 0) == allow_negative
    assert type(props['image']) is (ACAImage if aca_image else np.ndarray)


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
@pytest.mark.skipif('not IS_PY3', reason='Test built and tested only PY3')
def test_get_dark_cal_props_table():
    props = dark_cal.get_dark_cal_props_table('2007:001', '2008:001')
    assert np.allclose(props['eb'], [24.6, 25.89, 51.13, 1.9])
    assert props.colnames == ['ccd_temp', 'date', 'dec', 'dur', 'eb', 'el', 'id', 'l_l0',
                              'ra', 'start', 'stop', 'sun_el', 'zodib']


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
@pytest.mark.skipif('not IS_PY3', reason='Test built and tested only PY3')
def test_get_dark_cal_props_table_acdc():
    """Just acdc dark cals, giving a non-masked result
    """
    props = dark_cal.get_dark_cal_props_table('2019:150', '2019:160')
    assert not hasattr(props['t_ccd'], 'mask')
    assert props.colnames == ['date', 't_ccd', 'n_ccd_img', 'ccd_temp', 'datestart',
                              'datestop', 'filename']


@pytest.mark.skipif('not HAS_DARK_ARCHIVE', reason='Test requires dark archive')
@pytest.mark.skipif('not IS_PY3', reason='Test built and tested only PY3')
def test_get_dark_cal_props_table_mixed():
    """Mix of "classic" dark cals and acdc dark cals, giving a masked table
    """
    props = dark_cal.get_dark_cal_props_table('2019:001', '2019:160')
    assert np.allclose(props['eb'], [1.82, -49.53, -100])
    assert np.all(props['eb'].mask == [False, False, True])
    assert np.allclose(props['t_ccd'], [-100, -100, -11.086])
    assert np.all(props['t_ccd'].mask == [True, True, False])
    assert props.colnames == ['ccd_temp', 'date', 'dec', 'dur', 'eb', 'el', 'id', 'l_l0',
                              'ra', 'start', 'stop', 'sun_el', 'zodib', 't_ccd',
                              'n_ccd_img', 'datestart', 'datestop', 'filename']
