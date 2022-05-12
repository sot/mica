# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import shutil
import tempfile
from Chandra.Time import DateTime
import mica.common

# Override MICA_ARCHIVE with a temporary directory
TESTDIR = tempfile.mkdtemp()
mica.common.MICA_ARCHIVE = TESTDIR
# import mica.starcheck.starcheck after setting MICA_ARCHIVE
import mica.starcheck.process

# Just ingest files from the last couple of weeks or so
# And just check (implicitly) that the update script didn't raise any exceptions
config = mica.starcheck.process.DEFAULT_CONFIG
config['start'] = DateTime() - 30
mica.starcheck.process.update(config)
# Cleanup manually
shutil.rmtree(TESTDIR)
