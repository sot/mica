#!/usr/bin/env python
"""
Script to update Ska file archive aspect L1 products.  Module
also provides methods to retrieve the directory (or directories)
for an obsid.

This uses the obsid_archive module with a configuration specific
to the aspect L1 products.

"""
import os
import logging
import obsid_archive


# these columns are available in the headers of the fetched telemetry
# for this product (ASP L1) and will be included in the file lookup table
ARCHFILES_HDR_COLS = ('tstart', 'tstop', 'caldbver', 'content',
                      'ascdsver', 'revision', 'date')

mica_archive = os.environ.get('MICA_ARCHIVE') or '/data/aca/archive'
#config = ConfigObj('asp1.conf')
CONFIG = dict(data_root=os.path.join(mica_archive, 'asp1'),
              temp_root=os.path.join(mica_archive, 'temp'),
              sql_def='archfiles_asp_l1_def.sql',
              apstat_table='aspect_1',
              apstat_id='aspect_1_id',
              label='asp_l1',
              small='asp1{fidprops}',
              small_glob='*fidpr*',
              small_ver_regex='pcadf\d+N(\d{3})_',
              full='asp1',
              filecheck=False,
              cols=ARCHFILES_HDR_COLS,
              content_types=['ASPQUAL', 'ASPSOL', 'ACADATA', 'GSPROPS',
                             'GYRODATA', 'KALMAN', 'ACACAL', 'ACACENT',
                             'FIDPROPS', 'GYROCAL', 'ACA_BADPIX'])


def get_options():
    import argparse
    desc = \
"""
Run the update process to get new ASP L1 telemetry, save it in the Ska
file archive, and include it in the file lookup database.  This is intended
to be run as a cron task, and in regular processing, the update will fetch
and ingest all telemetry since the task's last run.  Options also provided
to fetch and ingest specific obsids and versions.

See the ``CONFIG`` in the asp_l1.py file and the config description in
obsid_archive for more information on the asp l1 default config if parameters
without command-line options need to be changed.
"""
    parser = argparse.ArgumentParser(description=desc)
    defaults = dict(CONFIG)
    parser.set_defaults(**defaults)
    parser.add_argument("--obsid",
                        type=int,
                        help="specific obsid to process")
    parser.add_argument("--version",
                        default='last',
                        help="specific processing version to retrieve")
    parser.add_argument("--firstrun",
                        action='store_true',
                        help="for archive init., ignore rev in aspect_1 table")
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--temp-root",
                        help="parent temp directory")
    parser.add_argument("--filecheck",
                        action="store_true",
                        help="for provisional data, download files and check"
                        + " that all are present.  If unset, proceed if dir"
                        + " exists")
    parser.add_argument("--rebuild",
                        action="store_true",
                        help="Allow update to rebuild archive from obsid 1")
    opt = parser.parse_args()
    return opt

# set up an archive object with default config for use by the other
# get_* methods
archive = obsid_archive.ObsArchive(CONFIG)


def get_dir(obsid):
    """
    Get ASP L1 directory for default/released products for an obsid.

      >>> from mica.archive import asp_l1
      >>> asp_l1.get_dir(2121)
      '/data/aca/archive/asp1/02/02121'

    :param obsid: obsid
    :returns: directory
    :rtype: string
    """
    return archive.get_dir(obsid)


def get_obs_dirs(obsid):
    """
    Get all ASP L1 directories for an obsid in the Ska file archive.

      >>> from mica.archive import asp_l1
      >>> obsdirs = asp_l1.get_obs_dirs(6000)

    obsdirs will look something like::

      {'default': '/data/aca/archive/asp1/06/06000',
      2: '/data/aca/archive/asp1/06/06000_v02',
      3: '/data/aca/archive/asp1/06/06000_v03',
      'last': '/data/aca/archive/asp1/06/06000',
      'revisions': [2, 3]}

    :param obsid: obsid
    :returns: map of obsid version to directories
    :rtype: dictionary
    """
    return archive.get_obs_dirs(obsid)


def get_files(obsid=None, start=None, stop=None,
              revision=None, content=None):
    """
    List asp_l1 files for an obsid or a time range.

      >>> from mica.archive import asp_l1
      >>> obs_files = asp_l1.get_files(6000)
      >>> obs_gspr = asp_l1.get_files(6000, content=['GSPROPS'])
      >>> range_fidpr = asp_l1.get_files(start='2012:001',
      ...                                stop='2012:030',
      ...                                content=['FIDPROPS'])

    :param obsid: obsid
    :param start: time range start (Chandra.Time compatible)
    :param stop: time range stop (Chandra.Time compatible)
    :param revision: revision integer or 'last'
                     defaults to current released version
    :param content: archive CONTENT type
                    defaults to all available ASP1 types
    :returns: full path of files matching query
    """
    return archive.get_files(obsid=obsid, start=start, stop=stop,
                             revision=revision, content=content)


def main():
    """
    Run the update process to get new ASP L1 telemetry, save it in the Ska
    file archive, and include it in the file lookup database.
    """
    opt = get_options()
    config = vars(opt)
    archive = obsid_archive.ObsArchive(config)
    archive.logger.setLevel(logging.INFO)
    archive.logger.addHandler(logging.StreamHandler())
    archive.update()

if __name__ == '__main__':
    main()
