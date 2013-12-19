"""
Common definitions for Mica subpackages
"""
import os


class MissingDataError(Exception):
    pass

FLIGHT_MICA_ARCHIVE = '/data/aca/archive'
MICA_ARCHIVE = os.environ.get('MICA_ARCHIVE', FLIGHT_MICA_ARCHIVE)
