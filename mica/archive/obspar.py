#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Script to update Ska file archive obspars.  Module
also provides methods to retrieve the directory (or directories)
for an obsid.

This uses the obsid_archive module with a configuration specific
to the obspar products.

"""
import os
import logging
from glob import glob
import re

from mica.archive import obsid_archive
from mica.archive.obsid_archive import parse_obspar
from mica.common import MICA_ARCHIVE

# these keys are available in the obspar and will be included in the
# file lookup table (archfiles)
ARCHFILES_HDR_COLS = [
'filename',
'obsid',
'title',
'observer',
'ao',
'object',
'ss_object',
'obs_id',
'obi_num',
'seq_num',
'instrume',
'grating',
'detector',
'detnam',
'si_mode',
'optical_monitor',
'raster_scan',
'dither',
'dither_y_amp',
'dither_y_freq',
'dither_y_phase',
'dither_z_amp',
'dither_z_freq',
'dither_z_phase',
'ra_pnt',
'dec_pnt',
'roll_pnt',
'ra_targ',
'dec_targ',
'y_det_offset',
'z_det_offset',
'radial_offset',
'defocus',
'sim_z_offset',
'pre_id',
'uninterrupted',
'seg_max_num',
'ra_nom',
'dec_nom',
'roll_nom',
'date_obs',
'date_end',
'tstart',
'tstop',
'sched_start',
'sched_stop',
'sched_exp_time',
'obs_mode',
'maneuver',
'maneuver_v1',
'maneuver_v2',
'maneuver_v3',
'maneuver_angle',
'maneuver_ref',
'mjdref',
'timezero',
'timeunit',
'timesys',
'timversn',
'datamode',
'readmode',
'ccdi0_on',
'ccdi1_on',
'ccdi2_on',
'ccdi3_on',
'ccds0_on',
'ccds1_on',
'ccds2_on',
'ccds3_on',
'ccds4_on',
'ccds5_on',
'dropped_chip_count',
'onchip_sum',
'sumrow',
'sumcol',
'subarray',
'startrow',
'rowcnt',
'subarray_frame_time',
'duty_cycle',
'dtycycle',
'exptimea',
'exptimeb',
'most_efficient',
'eventfilter',
'phamin',
'pharange',
'bias_request',
'mode',
'mission',
'telescop',
'sim_x',
'sim_y',
'sim_z',
'foc_len',
'py_shutter_position',
'range_switch_level',
'antico_enable',
'upper_level_disc',
'timing_mode',
'trigger_level',
'u_blank_hi',
'u_blank_low',
'v_blank_low',
'v_blank_hi',
'uld_enable',
'zero_block',
'blank_enable',
'width_enable',
'spect_mode',
'my_shutter_position',
'width_threshold',
'mt_a',
'mt_aop',
'mt_e',
'mt_epoch',
'mt_i',
'mt_ma',
'mt_raan',
'timeref',
'tassign',
'origin',
'ascdsver',
'obi_data',
'revision',
'obspar_ver',
'obspar_type',
'obspar_stat']

CONFIG = dict(data_root=os.path.join(MICA_ARCHIVE, 'obspar'),
              temp_root=os.path.join(MICA_ARCHIVE, 'tempobs'),
              sql_def='obspar_def.sql',
              apstat_table='obidet_0_5',
              apstat_id='obidet_0_5_id',
              label='obspar',
              small='obspar',
              small_glob='axaff*par*',
              small_ver_regex='axaff\d{5}_\d{3}N(\d{3}).*',
              full='obspar',
              filecheck=True,
              cols=ARCHFILES_HDR_COLS,
              content_types=['OBSPAR'])


def get_options():
    import argparse
    desc = \
"""
Run the update process to get new obspars, save them in the Ska
file archive, and include records of them in the file lookup database.
This is intended to be run as a cron task, and in regular processing,
the update will fetch and ingest all telemetry since the task's last run.
Options also provided to fetch and ingest specific obsids and versions.

See the ``CONFIG`` in the obspar.py file and the config description in
obsid_archive for more information on the obspar default config if parameters
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
                        help="for archive init., ignore revisions")
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--temp-root",
                        help="parent temp directory")
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
    Get obspar directory for default/released products for an obsid.

      >>> from mica.archive import obspar
      >>> obspar.get_dir(2121)
      '/proj/sot/ska/data/mica/archive/obspar/02/02121'

    :param obsid: obsid
    :returns: directory
    :rtype: string
    """
    return archive.get_dir(obsid)


def get_obs_dirs(obsid):
    """
    Get all obspar directories for an obsid in the Ska file archive.

      >>> from mica.archive import obspar
      >>> obsdirs = obspar.get_obs_dirs(6000)

    obsdirs will look something like::

      {'default': '/proj/sot/ska/data/mica/archive/obspar/06/06000',
      2: '/proj/sot/ska/data/mica/archive/obspar/06/06000_v02',
      3: '/proj/sot/ska/data/mica/archive/obspar/06/06000_v03',
      'last': '/proj/sot/ska/data/mica/archive/obspar/06/06000',
      'revisions': [2, 3]}

    :param obsid: obsid
    :returns: map of obsid version to directories
    :rtype: dictionary
    """
    return archive.get_obs_dirs(obsid)


def get_obspar_file(obsid, version='default'):
    """
    Get location of requested obsid's obspar.

      >>> from mica.archive import obspar
      >>> obspar.get_obspar_file(7000)
      '/proj/sot/ska/data/mica/archive/obspar/07/07000/axaff07000_000N002_obs0a.par.gz'
      >>> obspar.get_obspar_file(14262, version=1)
      '/proj/sot/ska/data/mica/archive/obspar/14/14262_v01/axaff14262_001N001_obs0a.par.gz'

    :param obsid: obsid
    :param version: processing version/revision

    :returns: path of obspar or None
    """
    dirs = archive.get_obs_dirs(obsid)
    if dirs is not None and version in dirs:
        obsparfiles = glob(os.path.join(dirs[version], 'axaff*par*'))
        if len(obsparfiles):
            return obsparfiles[0]
    return None


def get_obspar(obsid, version='default'):
    """
    Get the obspar for obsid.  Return as a dict.

      >>> from mica.archive import obspar
      >>> obspar.get_obspar(7001)['detector']
      'ACIS-I'

    :param obsid: obsid
    :param version: processing version/revision

    :returns: dictionary of obspar
    """

    obsparfile = get_obspar_file(obsid, version)
    if obsparfile is None:
        return None
    obspar = {'num_ccd_on': 0}
    for row in parse_obspar(obsparfile):
        obspar.update({row['name']: row['value']})
        if re.match(r'^ccd[is]\d_on$', row['name']) and row['value'] == 'Y':
            obspar['num_ccd_on'] += 1
    return obspar


def get_files(obsid=None, start=None, stop=None, revision=None):
    """
    List obspar files for an obsid or a time range.

      >>> from mica.archive import obspar
      >>> obs_files = obspar.get_files(6000)
      >>> range = obspar.get_files(start='2012:001',
      ...                          stop='2012:030')


    :param obsid: obsid
    :param start: time range start (Chandra.Time compatible)
    :param stop: time range stop (Chandra.Time compatible)
    :param revision: revision integer or 'last'
                     defaults to current released version
    :returns: full path of files matching query
    """
    return archive.get_files(obsid=obsid, start=start, stop=stop,
                             revision=revision)


def get_obsids(start=None, stop=None, revision=None):
    """
    List obsids with obspars in a time range.

      >>> from mica.archive import obspar
      >>> obsids = obspar.get_obsids('2015:001', '2015:003')

    To retrieve the complete set of obsids used over a time range,
    use kadi (kadi.events.obsid) instead of this obspar-based method.

    :param start: time range start (Chandra.Time compatible)
    :param stop: time range stop (Chandra.Time compatible)
    :param revision: None, revision integer, or 'last'
                     None is equivalent to limiting the search
                     to those with 'default'/released obspars.
    :returns: list of obsids
    """
    files = archive._get_file_records(start=start, stop=stop,
                                      revision=revision)
    if len(files) == 0:
        return []
    return files['obsid'].tolist()


def main():
    """
    Run the update process to get new obspars, save them in the Ska
    file archive, and include new entries in the file lookup database.
    """
    opt = get_options()
    config = vars(opt)
    archive = obsid_archive.ObsArchive(config)
    archive.logger.setLevel(logging.INFO)
    archive.logger.addHandler(logging.StreamHandler())
    archive.update()


if __name__ == '__main__':
    main()

