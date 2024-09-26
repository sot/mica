# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Files for ACA dark current access and processing
"""

MICA_FILES = {
    "dark_cals_dir": "aca_dark",
    "dark_cal_dir": "aca_dark/{{dark_cal.id}}",
    "dark_image": "aca_dark/{{dark_cal.id}}/image",
    "dark_props": "aca_dark/{{dark_cal.id}}/properties",
}
