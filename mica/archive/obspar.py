#!/usr/bin/env python
"""
Script to update Ska file archive obspars.  Module
also provides methods to retrieve the directory (or directories)
for an obsid.

This uses the obsid_archive module with a configuration specific
to the obspar products.

"""
import os
import logging
import obsid_archive
from glob import glob

# these keys are available in the obspar and will be included in the
# file lookup table (archfiles)
archfiles_hdr_cols = [
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

config = dict(data_root='/data/aca/archive/obspar',
              temp_root='/data/aca/archive/tempobs',
              sql_def='obspar_def.sql',
              apstat_table='obidet_0_5',
              apstat_id='obidet_0_5_id',
              label='obspar',
              small='obspar',
              small_glob='axaff*par*',
              small_ver_regex='axaff\d{5}_\d{3}N(\d{3}).*',
              full='obspar',
              cols=archfiles_hdr_cols)


def get_options():
    import argparse
    desc = \
"""
Run the update process to get new obspars, save them in the Ska
file archive, and include records of them in the file lookup database.
This is intended to be run as a cron task, and in regular processing,
the update will fetch and ingest all telemetry since the task's last run.
Options also provided to fetch and ingest specific obsids and versions.

See the ``config`` in the obspar.py file and the config description in
obsid_archive for more information on the obspar default config if parameters
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
                        help="for archive init., ignore revisions")
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--temp-root",
                        help="parent temp directory")
    opt = parser.parse_args()
    return opt

# set up an archive object with default config for use by the other
# get_* methods
archive = obsid_archive.ObsArchive(config)


def get_dir(obsid):
    """
    Get obspar directory for default/released products for an obsid.

    :param obsid: obsid
    :returns: directory
    :rtype: string
    """
    return archive.get_dir(obsid)

def get_obs_dirs(obsid):
    """
    Get all obspar directories for an obsid in the Ska file archive.

    :param obsid: obsid
    :returns: map of obsid version to directories
    :rtype: dictionary
    """
    return archive.get_obs_dirs(obsid)

def get_obspar(obsid, version='default'):
    """
    Get location of requested obsid's obspar.

    :param obsid: obsid
    :param version: processing version/revision

    :returns: path of obspar or None
    """
    dirs = archive.get_obs_dirs(obsid)
    if version in dirs:
        obsparfiles = glob(os.path.join(dirs[version], 'axaff*par*'))
        if len(obsparfiles):
            return obsparfiles[0]
    return None


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
