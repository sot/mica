"""
Files for ACA dark current access and processing
"""

import pyyaks.context    # Template rendering to provide context values

DARK_CAL = pyyaks.context.ContextDict('dark_cal')

SKA_FILES = pyyaks.context.ContextDict('ska_files', basedir='/proj/sot/ska/data/aca_dark_cal')
SKA_FILES.update({'dark_cal_dir': '{{dark_cal.id}}',
                  'image': '{{dark_cal.id}}/imd',
                  'info': '{{dark_cal.id}}/info'})

# Temporarily set default mica archive location to /tmp for safety.  In main() this gets
# set to os.path.join(opt.data_root, 'archive', 'aca_dark').
# One must set --data-root=/data/aca/archive/aca_dark for production.
MICA_FILES = pyyaks.context.ContextDict('mica_files', basedir='/tmp')
MICA_FILES.update({'dark_cal_dir': '{{dark_cal.id}}',
                   'image': '{{dark_cal.id}}/image',
                   'properties': '{{dark_cal.id}}/properties'})
