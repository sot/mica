# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Files for ACA dark current access and processing
"""
SKA_FILES = {'dark_cals_dir': 'data/aca_dark_cal',
             'dark_cal_dir': 'data/aca_dark_cal/{{dark_cal.id}}',
             'dark_image': 'data/aca_dark_cal/{{dark_cal.id}}/imd',
             'dark_info': 'data/aca_dark_cal/{{dark_cal.id}}/info'}

MICA_FILES = {'dark_cals_dir': 'aca_dark',
              'dark_cal_dir': 'aca_dark/{{dark_cal.id}}',
              'dark_image': 'aca_dark/{{dark_cal.id}}/image',
              'dark_props': 'aca_dark/{{dark_cal.id}}/properties'}
