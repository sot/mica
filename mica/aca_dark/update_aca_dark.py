#!/usr/bin/env python

from __future__ import division
import os
import shutil
import argparse
import json
import re

import numpy as np
from astropy.io import ascii

import kadi.events
from Chandra.Time import DateTime
import pyyaks.task       # Pipeline definition and execution
import pyyaks.logger     # Output logging control
from Ska.engarchive import fetch_sci as fetch

from mica.archive import aca_hdr3
from mica.common import MissingDataError
from .files import SKA_FILES, MICA_FILES, DARK_CAL

logger = None
ZODI_PROPS = None


def get_opt(args=None):
    parser = argparse.ArgumentParser(description='Update the aca_dark database')
    parser.add_argument("--start",
                        help="Processing start date default=now - 30 days")
    parser.add_argument("--log-level",
                        type=int,
                        default=pyyaks.logger.INFO,
                        help=("Logging level"))
    parser.add_argument("--data-root",
                        default=".",
                        help="Root data directory (default='.')")

    args = parser.parse_args(args)
    return args


def get_dark_cal_id(date):
    """
    Get the dark cal ID in the Ska dark current files corresponding to ``date``.  It is
    assumed the dark cal dirs are labeled by YYYYDOY within ``root``.

    This routine allows for plus/minus one day of slop.

    :param date: Date of dark cal
    :returns: str in yeardoy format
    """
    date0 = DateTime(date)
    for delta_day in (-1, 0, 1):
        date_try = (date0 + delta_day).date
        yeardoy = date_try[:4] + date_try[5:8]
        dark_cal_dir = os.path.join(SKA_FILES['dark_cals_dir'].abs, yeardoy)
        if os.path.exists(dark_cal_dir):
            return yeardoy

    raise MissingDataError('No dark calibration directory for {}'.format(date))


@pyyaks.task.task()
def get_id():
    """
    Get the YYYYDOY identifier for the dark cal.
    """
    DARK_CAL['id'] = get_dark_cal_id(DARK_CAL['start'].val)
    logger.verbose('Dark cal starting at {} has id={}'
                   .format(DARK_CAL['start'], DARK_CAL['id']))


@pyyaks.task.task()
@pyyaks.task.depends(depends=(SKA_FILES['dark_image.fits'],),
                     targets=(MICA_FILES['dark_image.fits'],))
def copy_dark_image():
    """
    Copy dark cal image from Ska to Mica
    """
    outdir = MICA_FILES['dark_cal_dir'].abs
    if not os.path.exists(outdir):
        logger.info('Making output dark cal directory {}'.format(outdir))
        os.makedirs(outdir)

    infile = SKA_FILES['dark_image.fits'].abs
    outfile = MICA_FILES['dark_image.fits'].abs
    logger.info('Copying {} to {}'.format(infile, outfile))
    shutil.copy(infile, outfile)


def get_ccd_temp(tstart, tstop):
    """
    Get the best estimate of CCD temperature between tstart and tstop

    :param tstart: start time in CXC seconds
    :param tstop: stop time in CXC seconds
    """
    # Just guess for pre-2000 data
    if tstart < DateTime('2000:001').secs:
        return -10, 'GUESS'

    # Try the HDR3 archive
    ccd_temp = aca_hdr3.Msid('ccd_temp', tstart, tstop)
    if len(ccd_temp.vals) > 4:
        return np.mean(ccd_temp.vals), 'HDR3'

    # Insufficient HDR3 data available, fall back to AACCCDPT and interpolate
    time0 = (tstart + tstop) / 2.0
    ccd_temp = fetch.Msid('AACCCDPT', time0 - 20000, time0 + 20000)
    if len(ccd_temp) > 100:
        x = ccd_temp.times - time0
        r = np.polyfit(x, ccd_temp.vals, 2)
        y_fit = np.polyval(r, 0.0)
        return y_fit, 'AACCCDPT'

    # Nothing worked.
    raise Exception('Unable to determine CCD temperature')


def get_zodi_props(dark_id):
    """
    Get zodiacal light information for given ``dark_id``.

    Parameters are: ('date', 'ra', 'dec', 'el', 'eb', 'sun_el', 'l_l0', 'zodib').

    :param dark_id: dark cal ID as YYYYDOY
    :returns: table Row object
    """
    global ZODI_PROPS
    dark_cals_dir = SKA_FILES['dark_cals_dir'].abs
    if ZODI_PROPS is None:
        dark_dir_files = [f for f in os.listdir(dark_cals_dir)
                          if re.search(r'[12]\d{6}$', f)]
        last_dark_id = sorted(dark_dir_files)[-1]
        filename = os.path.join(dark_cals_dir, last_dark_id, 'Result', 'zodi.csv')
        ZODI_PROPS = ascii.read(filename, delimiter=',', guess=False)

    date = '{}:{}'.format(dark_id[:4], dark_id[4:])
    for zodi_prop in ZODI_PROPS:
        if zodi_prop['date'] == date:
            return zodi_prop
    else:
        raise MissingDataError('No Zodiacal properties found for {}'.format(date))


@pyyaks.task.task()
@pyyaks.task.depends(depends=(MICA_FILES['dark_image.fits'],),
                     targets=(MICA_FILES['dark_props.json'],))
def make_properties():
    """
    Compute basic observation properties and store in a local JSON file.
    """
    props = {key: DARK_CAL[key].val for key in DARK_CAL.keys()}

    # First get the individual replicas and set properties
    props['replicas'] = []

    start = DateTime(DARK_CAL['start'].val)
    stop = DateTime(DARK_CAL['stop'].val)
    replicas = kadi.events.dark_cal_replicas.filter(start=start - 0.5, stop=start + 1)
    for i_replica, replica in enumerate(replicas):
        # Get the CCD temperature within +/- 10 minutes of replica
        ccd_temp, temp_source = get_ccd_temp(replica.tstart - 600, replica.tstop + 600)

        # Determine if replica was successfully downlinked.  Prior to 2001:080 the file
        # structure was inconsistent so it is hard to know algorithmically, but for these
        # cases all replicas were good.  Otherwise look for VC2_Replica<N>_SFDU
        if start.date < '2001:080':
            downlink = True
        else:
            dark_files = os.listdir(SKA_FILES['dark_cal_dir'].abs)
            replica_match = 'VC2_Replica{}_SFDU'.format(i_replica + 1)
            downlink = any(replica_match in filename for filename in dark_files)

        replica_props = {'count': i_replica + 1,
                         'ccd_temp': ccd_temp,
                         'ccd_temp_source': temp_source,
                         'start': start.date,
                         'stop': stop.date,
                         'downlink': downlink}
        props['replicas'].append(replica_props)

    # Add ccd temperature, zodiacal light brightness, attitude, sun position
    props['ccd_temp'] = np.mean([repl['ccd_temp'] for repl in props['replicas']])
    zodi_props = get_zodi_props(DARK_CAL['id'].val)
    props.update({key: zodi_props[key].tolist() for key in zodi_props.colnames})

    outfile = MICA_FILES['dark_props.json'].abs
    logger.info('Writing dark cal replica properties to {}'.format(outfile))
    with open(outfile, 'w') as fh:
        json.dump(props, fh, sort_keys=True, indent=4, separators=(',', ': '))


def main():
    """
    Update all the Mica dark cal directories.
    """
    global logger

    opt = get_opt()
    logger = pyyaks.logger.get_logger(level=opt.log_level)

    MICA_FILES.basedir = os.path.join(opt.data_root, 'archive', 'aca_dark')

    # Get dark cals after the processing start time (or now - 30 days) from kadi
    start = DateTime() - 30 if opt.start is None else DateTime(opt.start)
    dark_events = kadi.events.dark_cals.filter(start=start)

    for dark_event in dark_events:
        # Clear previous DARK_CAL values and set DARK_CAL from the kadi dark cal event
        DARK_CAL.clear()
        for attr in ('start', 'stop', 'dur'):
            DARK_CAL[attr] = getattr(dark_event, attr)

        # Start the processing loop
        process_msg = 'Processing dark cal at date {}'.format(dark_event.start)
        pyyaks.task.start(message=process_msg)

        # Do the actual pipeline processing events
        get_id()
        copy_dark_image()
        make_properties()

        pyyaks.task.end(message=process_msg)

if __name__ == '__main__':
    main()
