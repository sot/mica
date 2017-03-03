import os
import tempfile
from Chandra.Time import DateTime
from Ska.Shell import bash
import mica.common

# Override MICA_ARCHIVE with a temporary directory
TESTDIR = tempfile.mkdtemp()
mica.common.MICA_ARCHIVE = TESTDIR
# import mica.starcheck.starcheck after setting MICA_ARCHIVE
import mica.starcheck.starcheck

# Just ingest files from the last couple of weeks or so
# This still uses the silly find files newer than this other file method, so
# set the time stamp on that reference file
if not os.path.exists(os.path.join(TESTDIR, 'starcheck')):
    os.makedirs(os.path.join(TESTDIR, 'starcheck'))
bash("touch -d {} {}".format(DateTime(-15).iso, mica.starcheck.starcheck.FILES['touch_file']))
# And just check that the update script didn't raise any exceptions
mica.starcheck.starcheck.update()
# Cleanup manually
bash("rm -r {}".format(TESTDIR))
