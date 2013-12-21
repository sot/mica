"""
Common definitions for Mica subpackages
"""
import os


class MissingDataError(Exception):
    pass

FLIGHT_MICA_ARCHIVE = '/data/aca/archive'

# The MICA_ARCHIVE env. var can be a colon-delimited path, which allows
# packages using pyyaks context files to see files within multiple trees.
# This can be convenient for testing and development.  Normally the flight
# path is put last so that the test path is preferred if available.
MICA_ARCHIVE_PATH = os.environ.get('MICA_ARCHIVE', FLIGHT_MICA_ARCHIVE)

# Most of the existing subpackages just expect a single path, which should
# be the test path (first one) if defined that way.
MICA_ARCHIVE = MICA_ARCHIVE_PATH.split(':')[-1]
