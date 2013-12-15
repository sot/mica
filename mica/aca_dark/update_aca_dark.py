#!/usr/bin/env python

from __future__ import division
import os
import shutil
import argparse

import kadi.events
from Chandra.Time import DateTime
import pyyaks.task       # Pipeline definition and execution
import pyyaks.logger     # Output logging control
import pyyaks.context    # Template rendering to provide context values


DARK_CAL = pyyaks.context.ContextDict('dark_cal')

SKA_DARK_CAL = '/proj/sot/ska/data/aca_dark_cal'
SKA_FILES = pyyaks.context.ContextDict('ska_files', basedir=SKA_DARK_CAL)
SKA_FILES.update({'image': '{{dark_cal.id}}/imd',
                  'info': '{{dark_cal.id}}/info'})

# Temporarily set default mica archive location to /tmp for safety.  In main() this gets
# set to os.path.join(opt.data_root, 'archive', 'aca_dark').
# One must set --data-root=/data/aca/archive/aca_dark for production.
MICA_FILES = pyyaks.context.ContextDict('mica_files', basedir='/tmp')
MICA_FILES.update({'dark_cal_dir': '{{dark_cal.id}}',
                   'image': '{{dark_cal.id}}/image',
                   'properties': '{{dark_cal.id}}/properties'})

logger = None


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


def set_dark_cal_id(date, rootdir=SKA_DARK_CAL):
    """
    Get the dark cal ID in ``rootdir`` corresponding to ``date``.  It is assumed
    the dark cal dirs are labeled by YYYYDOY within ``root``.

    This routine allows for plus/minus one day of slop.

    :param date: Date of dark cal
    :returns: str in yeardoy format
    """
    date0 = DateTime(date)
    for delta_day in (-1, 0, 1):
        date_try = (date0 + delta_day).date
        yeardoy = date_try[:4] + date_try[5:8]
        dark_cal_dir = os.path.join(rootdir, yeardoy)
        if os.path.exists(dark_cal_dir):
            return yeardoy

    raise ValueError('No dark calibration directory for {}'.format(date))


@pyyaks.task.task()
def get_id():
    """
    Get the YYYYDOY identifier for the dark cal.
    """
    DARK_CAL['id'] = set_dark_cal_id(DARK_CAL['start'].val)
    logger.verbose('Dark cal starting at {} has id={}'
                   .format(DARK_CAL['start'], DARK_CAL['id']))


@pyyaks.task.task()
@pyyaks.task.depends(depends=(SKA_FILES['image.fits'],),
                     targets=(MICA_FILES['image.fits'],))
def copy_dark_image():
    """
    Copy dark cal image from Ska to Mica
    """
    outdir = MICA_FILES['dark_cal_dir'].abs
    if not os.path.exists(outdir):
        logger.info('Making output dark cal directory {}'.format(outdir))
        os.makedirs(outdir)

    infile = SKA_FILES['image.fits'].abs
    outfile = MICA_FILES['image.fits'].abs
    logger.info('Copying {} to {}'.format(infile, outfile))
    shutil.copy(infile, outfile)


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


if __name__ == '__main__':
    main()
