#!/usr/bin/env python

import os
import logging
import obsid_archive

#from configobj import ConfigObj
#config = ConfigObj("obspar.conf")
config = dict(data_root='/data/aca/archive/obspar',
              temp_root='/data/aca/archive/tempobs',
              sql_def='obspar_def.sql',
              apstat_table='obidet_0_5',
              apstat_id='obidet_0_5_id',
              label='obspar',
              small='obspar',
              small_glob='axaff*par*',
              small_ver_regex='axaff\d{5}_\d{3}N(\d{3}).*',
              full='obspar')


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

def get_options():
#    from optparse import OptionParser
#    parser = OptionParser()
    import argparse
    parser = argparse.ArgumentParser(
        description="Fetch aspect level 1 products and make a file archive")
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
    parser.add_argument("--proctype",
                        default="asp1")
    opt = parser.parse_args()
    return opt 


def get_dir(obsid):
    archive = obsid_archive.ObsArchive(dict(config))
    return archive.get_dir(obsid)

def get_obs_dirs(obsid):
    archive = obsid_archive.ObsArchive(dict(config))
    return archive.get_obs_dirs(obsid)

def main():
    opt = get_options()
    config = dict(opt.__dict__, cols=archfiles_hdr_cols)
    archive = obsid_archive.ObsArchive(config)
    archive.logger.setLevel(logging.INFO)
    archive.logger.addHandler(logging.StreamHandler())
    archive.update()

if __name__ == '__main__':
    main()
