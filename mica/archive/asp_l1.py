#!/usr/bin/env python
"""
Script to update Ska file archive aspect L1 products.  Module
also provides methods to retrieve the directory (or directories)
for an obsid.

This uses the obsid_archive module with a configuration specific
to the aspect L1 products.

"""

import logging
import obsid_archive


# these columns are available in the headers of the fetched telemetry
# for this product (ASP L1) and will be included in the file lookup table
archfiles_hdr_cols = ('tstart', 'tstop', 'caldbver', 'content',
                      'ascdsver', 'revision', 'date')

#config = ConfigObj('asp1.conf')
config = dict(data_root='/data/aca/archive/asp1',
              temp_root='/data/aca/archive/temp',
              sql_def='archfiles_asp_l1_def.sql',
              apstat_table='aspect_1',
              apstat_id='aspect_1_id',
              label='asp_l1',
              small='asp1{fidprops}',
              small_glob='*fidpr*',
              small_ver_regex='pcadf\d+N(\d{3})_',
              full='asp1',
              filecheck=False,
              cols=archfiles_hdr_cols,
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

See the ``config`` in the asp_l1.py file and the config description in
obsid_archive for more information on the asp l1 default config if parameters
without command-line options need to be changed.
"""
    parser = argparse.ArgumentParser(description=desc)
    defaults = dict(config)
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
    opt = parser.parse_args()
    return opt

# set up an archive object with default config for use by the other
# get_* methods
archive = obsid_archive.ObsArchive(config)


def get_dir(obsid):
    """
    Get ASP L1 directory for default/released products for an obsid.

    :param obsid: obsid
    :returns: directory
    :rtype: string
    """
    return archive.get_dir(obsid)


def get_obs_dirs(obsid):
    """
    Get all ASP L1 directories for an obsid in the Ska file archive.

    :param obsid: obsid
    :returns: map of obsid version to directories
    :rtype: dictionary
    """
    return archive.get_obs_dirs(obsid)

def get_files(obsid, revision=None, content=None):
    return archive.get_files(obsid, revision=revision, content=content)

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
