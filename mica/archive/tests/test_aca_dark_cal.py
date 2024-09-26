# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Basic functionality and regression tests for ACA dark cal module.
"""

import os

import cxotime
import numpy as np
import pytest
from chandra_aca.aca_image import ACAImage

from mica.archive.aca_dark import dark_cal
from mica.common import MissingDataError

HAS_DARK_ARCHIVE = os.path.exists(dark_cal.MICA_FILES["dark_cals_dir"].abs)


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_date_to_dark_id():
    assert dark_cal.date_to_dark_id("2011-01-15T12:00:00") == "2011015"


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_dark_id_to_date():
    assert dark_cal.dark_id_to_date("2011015") == "2011:015"


def test_dark_temp_scale():
    scale = dark_cal.dark_temp_scale(-10.0, -14)
    assert np.allclose(scale, 0.70)

    scale = dark_cal.dark_temp_scale(-10.0, -14, scale_4c=2.0)
    assert scale == 0.5  # Should be an exact match


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_get_dark_cal_id():
    assert dark_cal.get_dark_cal_id("2007:008:12:00:00", "nearest") == "2007006"
    assert dark_cal.get_dark_cal_id("2007:008:12:00:00", "before") == "2007006"
    assert dark_cal.get_dark_cal_id("2007:008:12:00:00", "after") == "2007069"

    dark_cal_ids = list(dark_cal.get_dark_cal_ids().values())
    # removing these two to make sure it is not like the default case
    dark_cal_ids.remove("2007006")
    dark_cal_ids.remove("2007069")

    assert (
        dark_cal.get_dark_cal_id(
            "2007:008:12:00:00", "nearest", dark_cal_ids=dark_cal_ids
        )
        == "2006329"
    )
    assert (
        dark_cal.get_dark_cal_id(
            "2007:008:12:00:00", "before", dark_cal_ids=dark_cal_ids
        )
        == "2006329"
    )
    assert (
        dark_cal.get_dark_cal_id(
            "2007:008:12:00:00", "after", dark_cal_ids=dark_cal_ids
        )
        == "2007251"
    )


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
@pytest.mark.parametrize("allow_negative", [True, False])
@pytest.mark.parametrize("aca_image", [True, False])
def test_get_dark_cal_image(aca_image, allow_negative):
    image = dark_cal.get_dark_cal_image(
        "2007:008:12:00:00", aca_image=aca_image, allow_negative=allow_negative
    )
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


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
@pytest.mark.parametrize("allow_negative", [True, False])
@pytest.mark.parametrize("aca_image", [True, False])
def test_get_dark_cal_props(aca_image, allow_negative):
    props = dark_cal.get_dark_cal_props("2007:008:12:00:00")
    assert len(props["replicas"]) == 5
    assert props["start"] == "2007:006:01:56:46.817"

    props = dark_cal.get_dark_cal_props(
        "2007:008:12:00:00",
        include_image=True,
        aca_image=aca_image,
        allow_negative=allow_negative,
    )
    assert len(props["replicas"]) == 5
    assert props["start"] == "2007:006:01:56:46.817"
    assert props["image"].shape == (1024, 1024)
    assert np.any(props["image"] < 0) == allow_negative
    assert type(props["image"]) is (ACAImage if aca_image else np.ndarray)


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_get_dark_cal_props_table():
    props = dark_cal.get_dark_cal_props_table("2007:001:12:00:00", "2008:001:12:00:00")
    assert np.allclose(props["eb"], [24.6, 25.89, 51.13, 1.9])
    assert props.colnames == [
        "ccd_temp",
        "date",
        "dec",
        "dur",
        "eb",
        "el",
        "id",
        "l_l0",
        "ra",
        "start",
        "stop",
        "sun_el",
        "zodib",
    ]


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_get_dark_cal_props_table_acdc():
    """Just acdc dark cals, giving a non-masked result"""
    props = dark_cal.get_dark_cal_props_table("2019:150:12:00:00", "2019:160:12:00:00")
    assert not hasattr(props["t_ccd"], "mask")
    assert props.colnames == [
        "date",
        "t_ccd",
        "n_ccd_img",
        "ccd_temp",
        "datestart",
        "datestop",
        "filename",
    ]


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_get_dark_cal_props_table_mixed():
    """Mix of "classic" dark cals and acdc dark cals, giving a masked table"""
    props = dark_cal.get_dark_cal_props_table("2019:001:12:00:00", "2019:160:12:00:00")
    assert np.allclose(props["eb"], [1.82, -49.53, -100])
    assert np.all(props["eb"].mask == [False, False, True])
    assert np.allclose(props["t_ccd"], [-100, -100, -11.086])
    assert np.all(props["t_ccd"].mask == [True, True, False])
    assert props.colnames == [
        "ccd_temp",
        "date",
        "dec",
        "dur",
        "eb",
        "el",
        "id",
        "l_l0",
        "ra",
        "start",
        "stop",
        "sun_el",
        "zodib",
        "t_ccd",
        "n_ccd_img",
        "datestart",
        "datestop",
        "filename",
    ]


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_scalar():
    dark_cal_id = dark_cal.get_dark_cal_id("2022:100")
    assert dark_cal_id[:4] == "2022" and dark_cal_id[4:] == "069"


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_vectorized():
    dark_id_ref = ["2022100", "2022105", "2022127"]
    date_ref = ["2022:100", "2022:105", "2022:127"]
    assert np.all(dark_id_ref == dark_cal.date_to_dark_id(date_ref))
    assert np.all(date_ref == dark_cal.dark_id_to_date(dark_id_ref))

    dark_cal_id_ref = ["2022069", "2022104", "2022104"]
    dark_cal_id = dark_cal.get_dark_cal_id(
        ["2022:100", "2022:105", "2022:127"], "before"
    )
    assert np.all(dark_cal_id_ref == dark_cal_id)

    dark_cal_id_ref = ["2022104", "2022133", "2022133"]
    dark_cal_id = dark_cal.get_dark_cal_id(
        ["2022:100", "2022:105", "2022:127"], "after"
    )
    assert np.all(dark_cal_id_ref == dark_cal_id)

    dark_cal_id_ref = ["2022104", "2022104", "2022133"]
    dark_cal_id = dark_cal.get_dark_cal_id(
        ["2022:100", "2022:105", "2022:127"], "nearest"
    )
    assert np.all(dark_cal_id_ref == dark_cal_id)


@pytest.mark.skipif("not HAS_DARK_ARCHIVE", reason="Test requires dark archive")
def test_limits():
    dark_cal_ids = dark_cal.get_dark_cal_ids()
    first = cxotime.CxoTime(list(dark_cal_ids.keys())[0]) - 1 * cxotime.units.day
    with pytest.raises(MissingDataError, match="No dark cal found before"):
        dark_cal.get_dark_cal_id(first, "before")
    last = cxotime.CxoTime(list(dark_cal_ids.keys())[-1]) + 1 * cxotime.units.day
    with pytest.raises(MissingDataError, match="No dark cal found after"):
        dark_cal.get_dark_cal_id(last, "after")
