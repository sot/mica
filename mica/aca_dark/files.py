"""
Files for ACA dark current access and processing
"""
import os

import pyyaks.context    # Template rendering to provide context values

from mica.common import MICA_ARCHIVE

DARK_CAL = pyyaks.context.ContextDict('dark_cal')

SKA_FILES = pyyaks.context.ContextDict('ska_files', basedir='/proj/sot/ska/data/aca_dark_cal')
SKA_FILES.update({'dark_cals_dir': '',
                  'dark_cal_dir': '{{dark_cal.id}}',
                  'dark_image': '{{dark_cal.id}}/imd',
                  'dark_info': '{{dark_cal.id}}/info'})

MICA_FILES = pyyaks.context.ContextDict('mica_files',
                                        basedir=os.path.join(MICA_ARCHIVE, 'aca_dark'))
MICA_FILES.update({'dark_cals_dir': '',
                   'dark_cal_dir': '{{dark_cal.id}}',
                   'dark_image': '{{dark_cal.id}}/image',
                   'dark_props': '{{dark_cal.id}}/properties'})
