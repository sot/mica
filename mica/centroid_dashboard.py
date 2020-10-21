#!/usr/bin/env python

import os
import numpy as np
from glob import glob

# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pylab as plt

from kadi import events
from chandra_aca.centroid_resid import CentroidResiduals
from chandra_aca.transform import yagzag_to_pixels
from chandra_aca.plot import plot_stars
from mica.starcheck import get_mp_dir, get_starcat, get_att
from agasc import get_star

from Ska.Matplotlib import plot_cxctime
from Ska.engarchive import fetch
from Quaternion import Quat
from Chandra.Time import DateTime
from Ska.Numpy import interpolate

from astropy.table import Table, vstack
from astropy.io import ascii

import argparse
import pyyaks.logger

GUIDE_METRICS_OBSID = 'guide_metrics_obsid.dat'
GUIDE_METRICS_SLOT = 'guide_metrics_slot.dat'
MICA_REPORTS = "https://icxc.harvard.edu/aspect/mica_reports/"
MICA_PORTAL = "http://kadi.cfa.harvard.edu/mica/?obsid_or_date="
CD_ROOT = "/proj/sot/ska/www/ASPECT_ICXC/centroid_reports"

# Set up logging
loglevel = pyyaks.logger.INFO
logger = pyyaks.logger.get_logger(name='centroid_dashboard',
                                  level=loglevel,
                                  format="%(asctime)s %(message)s")

# Update guide metrics file with new obsids between NOW and (NOW - NDAYS) days
NDAYS = 7


def get_opt(args=None):
    parser = argparse.ArgumentParser(description='Centroid dashboard')
    parser.add_argument("--obsid",
                        help="Processing obsid (default=None)")
    parser.add_argument("--start",
                        help=f"Processing start date (default=NOW - {NDAYS} days)")
    parser.add_argument("--stop",
                        help="Processing stop date (default=NOW)")
    parser.add_argument("--data-root",
                        default=CD_ROOT,
                        help=f"Root directory for data files (default='{CD_ROOT}')")
    add_bool_arg(parser=parser,
                 name='make_plots',
                 help_="make plots if the option is included")
    add_bool_arg(parser=parser,
                 name='save',
                 help_="save plots if the option is included")
    add_bool_arg(parser=parser,
                 name='force',
                 help_="force processing if the option is included")
    return parser.parse_args(args)


def add_bool_arg(parser, name, default=True, help_=''):
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(f'--{name}',
                       dest=name,
                       action='store_true',
                       help=help_)
    group.add_argument(f'--no-{name}',
                       dest=name,
                       action='store_false',
                       help=f"Don't {help_}")
    parser.set_defaults(**{name: default})


def get_dr_dp_dy(refs, atts):
    att_errors = {}
    drs = []
    dps = []
    dys = []

    for ref_q, att_q in zip(refs, atts):
        dq = Quat(ref_q).dq(att_q)
        drs.append(dq.roll0 * 3600)
        dps.append(dq.pitch * 3600)
        dys.append(dq.yaw * 3600)

    att_errors['dr'] = np.array(drs)
    att_errors['dp'] = np.array(dps)
    att_errors['dy'] = np.array(dys)

    return att_errors


def get_observed_att_errors(obsid, crs=None, on_the_fly=False):
    """
    Get OBC pitch, yaw, roll errors with respect to the ground
    aspect solution

    :param obsid: obsid
    :param crs: dictionary with keys 'ground' and 'obc'. Values are also
                dictionaries keyd by slot number containing corresponding
                CentroidResiduals objects.
    :param on_the_fly:

    on_the_fly determins whether centroids residuals are provided
    or need to computed for the requested obsid
    """

    try:
        events.dwells.filter(obsid=obsid)[0]
    except Exception as err:
        logger.info(f'ERROR: {err}')
        return None

    att_errors = {'dr': [], 'dy': [], 'dp': []}

    # flag = 0 OK, ground aspect solution
    # flag = 1 'ground' is None for all slots, use obc solution
    # flag = 2 no common times between ground vs obc times
    # flag = 3 times could not be matched for yet another reason
    flag = 0

    if on_the_fly:
        crs = get_crs_per_obsid(obsid)
    else:
        if crs is None:
            raise Exception('Provide crs if on_the_fly is False')

    if all([item is None for item in crs['ground'].values()]):
        # No good ground solution in any of the slots
        # Use obc aspect solution
        flag = 1
        crs_ref = crs['obc']
        logger.info(
            f'No ground aspect solution for {obsid}. Using obc aspect solution for reference')
    else:
        crs_ref = crs['ground']

    try:
        # Adjust time axis if needed
        # TODO: what if slot 3 is not tracking a star?

        ref_att_times = crs_ref[3].att_times
        obc_att_times = crs['obc'][3].att_times

        ii = np.in1d(ref_att_times, obc_att_times)
        ref_att_times_adjusted = ref_att_times[ii]
        att_ref = crs_ref[3].atts[ii]

        if len(ref_att_times_adjusted) == 0:
            # no common times for obc and grnd solutions
            flag = 2
            raise ValueError('No common time vals for obc and ground att times')

        idx = list(obc_att_times).index(ref_att_times_adjusted[0])
        obc_att_times_adjusted = obc_att_times[idx: len(ref_att_times_adjusted) + idx]
        att_obc = crs['obc'][3].atts[idx: len(ref_att_times_adjusted) + idx]

        if not np.all(obc_att_times_adjusted == ref_att_times_adjusted):
            flag = 3
            raise ValueError('Could not align obc and ref aspect solution times')

    except Exception as err:
        logger.info(f'ERROR get_observed_att_errors: {err}')
        return {'obsid': obsid,
                'flag': flag}

    # Compute attitude errors
    att_errors = get_dr_dp_dy(att_ref, att_obc)

    out = {'obsid': obsid,
           'time': obc_att_times_adjusted,
           'dr': att_errors['dr'],
           'dy': att_errors['dy'],
           'dp': att_errors['dp'],
           'flag': flag,
           'crs': crs}

    return out


def get_crs_per_obsid(obsid):
    """
    Get OBC centroid residuals per obsid for all slots, with respect
    to both the OBC and ground aspect solution

    :param obsid: obsid
    """
    crs = {'ground': {}, 'obc': {}}
    att_sources = ['obc'] if int(obsid) > 40000 else ['ground', 'obc']

    cat = get_starcat(obsid)

    if cat is not None:

        ok = (cat['type'] == 'BOT') | (cat['type'] == 'GUI')
        cat = cat[ok]

        cols = ['obsid', 'idx', 'slot', 'id', 'type', 'sz', 'mag', 'yang', 'zang']
        crs['cat'] = cat[cols]

        slots = cat['slot']
        crs['slots'] = np.array(slots)

        for att_source in att_sources:
            for slot in slots:
                try:
                    cr = CentroidResiduals.for_slot(obsid=obsid,
                                                    slot=slot,
                                                    att_source=att_source,
                                                    centroid_source='obc')
                    crs[att_source][slot] = cr
                except Exception:
                    crs[att_source][slot] = None
                    logger.info(
                        f'Could not compute crs for {obsid} slot {slot} (att_source={att_source})')

    return crs


def get_observed_metrics(obsid):
    """
    Fetch manvr angle, one shot updates and aberration corrections,
    calculate centroid residuals and observed OBC roll error with
    respect to the ground (obc) solution for science (ER) observations
    as a function of time. Calculate also 50th and 95th percentile of
    the roll error, and log the preceding and next obsid.

    :param obsid: obsid
    """

    # One shot
    manvr = events.manvrs.filter(obsid=obsid)[0]
    one_shot = manvr.one_shot
    one_shot_pitch = manvr.one_shot_pitch
    one_shot_yaw = manvr.one_shot_yaw

    # Manvr angle
    manvr_angle = manvr.angle

    # Next obsid
    manvr_next = manvr.get_next()
    if manvr_next:
        obsid_next = manvr_next.get_obsid()
    else:
        obsid_next = -9999

    # Preceding obsid
    manvr_preceding = manvr.get_previous()
    if manvr_preceding:
        obsid_preceding = manvr_preceding.get_obsid()
    else:
        obsid_preceding = -9999

    # Attitude errors
    att_errors = get_observed_att_errors(obsid, on_the_fly=True)
    att_flag = att_errors['flag']

    if att_flag > 1:
        return ({'obsid': obsid,
                 'obsid_preceding': obsid_preceding,
                 'obsid_next': obsid_next,
                 'att_flag': att_flag,
                 'dwell': True},
                None)

    if att_errors is None:
        return ({'obsid': obsid,
                 'obsid_preceding': obsid_preceding,
                 'obsid_next': obsid_next,
                 'att_flag': att_flag,
                 'dwell': False},
                None)

    # dr statistics
    drs = att_errors['dr']
    dr50 = float(np.percentile(np.abs(drs), 50))
    dr95 = float(np.percentile(np.abs(drs), 95))
    mean_date = DateTime(0.5 * (att_errors['time'][0] + att_errors['time'][-1])).date

    # Centroid residuals
    crs = att_errors['crs']

    # Ending roll error prior to the manvr
    att_errors_preceding = get_observed_att_errors(obsid_preceding, on_the_fly=True)

    if att_errors_preceding is None:
        ending_roll_err = -9999
    else:
        ending_roll_err = att_errors_preceding['dr'][-1]

    # Aberration correction
    aber_flag = 0
    path_ = get_mp_dir(obsid)[0]
    if path_ is None:
        logger.info(f'No mp_dir for {obsid}. Skipping aber correction')
        aber_y = -9999
        aber_z = -9999
        aber_flag = 1
    else:
        mp_dir = f"/data/mpcrit1/mplogs/{path_}"
        manerr = glob(f'{mp_dir}/*ManErr.txt')

        if len(manerr) == 0:
            logger.info(f'No ManErr file for {obsid}. Skipping aber correction')
            aber_y = -9999
            aber_z = -9999
            aber_flag = 2
        else:
            dat = ascii.read(manerr[0], header_start=2, data_start=3)
            ok = dat['obsid'] == obsid

            if np.sum(ok) > 1:
                logger.info(f'More than one entry per {obsid}. Skipping aber correction')
                aber_y = -9999
                aber_z = -9999
                aber_flag = 3
            else:
                aber_y = dat['aber-Y'][ok][0]
                aber_z = dat['aber-Z'][ok][0]

    if aber_flag == 0:
        one_shot_aber_corrected = np.sqrt((one_shot_pitch - aber_y)**2
                                          + (one_shot_yaw - aber_z)**2)
    else:
        one_shot_aber_corrected = -9999

    out_obsid = {'obsid': obsid,
                 'mean_date': mean_date,
                 'dr50': dr50,
                 'dr95': dr95,
                 'one_shot': one_shot,
                 'one_shot_pitch': one_shot_pitch,
                 'one_shot_yaw': one_shot_yaw,
                 'manvr_angle': manvr_angle,
                 'obsid_preceding': obsid_preceding,
                 'ending_roll_err': ending_roll_err,
                 'aber_y': aber_y,
                 'aber_z': aber_z,
                 'aber_flag': aber_flag,
                 'one_shot_aber_corrected': one_shot_aber_corrected,
                 'obsid_next': obsid_next,
                 'att_errors': att_errors,
                 'att_flag': att_flag,
                 'dwell': True}

    out_slot = {'obsid': obsid, 'slots': {k: {} for k in range(8)}}

    cat = crs['cat']
    d = events.dwells.filter(obsid=obsid)[0]

    for slot in range(8):
        ok = cat['slot'] == slot
        out = {}
        if len(cat[ok]) > 0:
            out['id'] = cat['id'][ok][0]
            out['type'] = cat['type'][ok][0]
            out['mag'] = cat['mag'][ok][0]
            out['yang'] = cat['yang'][ok][0]
            out['zang'] = cat['zang'][ok][0]
            out['slot'] = slot

            if att_flag == 0:
                # Ground solution exists
                val = crs['ground'][slot]
            else:
                # Use obc solution
                val = crs['obc'][slot]

            if len(val.dyags) > 0 and len(val.dzags) > 0:
                out['std_dy'] = np.std(val.dyags)
                out['std_dz'] = np.std(val.dzags)
                out['rms_dy'] = np.sqrt(np.mean(val.dyags ** 2))
                out['rms_dz'] = np.sqrt(np.mean(val.dzags ** 2))
                out['median_dy'] = np.median(val.dyags)
                out['median_dz'] = np.median(val.dzags)
                drs = np.sqrt((val.dyags ** 2) + (val.dzags ** 2))
                for dist in ['1.5', '3.0', '5.0']:
                    out[f'f_within_{dist}'] = np.count_nonzero(drs < float(dist)) / len(drs)
            else:
                for metric in ['std_dy', 'std_dz', 'rms_dy', 'rms_dz', 'median_dy', 'median_dz']:
                    out[metric] = -9999
                for metric in ['f_within_5.0', 'f_within_3.0', 'f_within_1.5']:
                    out[metric] = 0

            mags = fetch.Msid(f'aoacmag{slot}', start=d.start, stop=d.stop)
            out['median_mag'] = np.median(mags.vals)

            out_slot['slots'][slot] = out

    return out_obsid, out_slot


def get_n_kalman(start, stop):
    """
    Get the AOKALSTR data with number of kalman stars reported by OBC
    """
    start = DateTime(start).date
    stop = DateTime(stop).date
    dat = fetch.Msid('aokalstr', start, stop)
    dat.interpolate(1.025)
    return dat


def get_cd_dir(obsid, data_root):
    """
    Check if the centroid dashbord directory exists for the requested obsid,
    and create it if needed
    """
    if obsid == -1:
        return ""

    cd_obsid_root = os.path.join(data_root, np.str(obsid)[:2], f"{obsid}")

    if not os.path.exists(cd_obsid_root):
        os.makedirs(cd_obsid_root)
        logger.info(f'Creating directory {cd_obsid_root}')

    return cd_obsid_root


def plot_att_errors_per_obsid(obsid, plot_dir, coord='dr', att_errors=None,
                              save=False, on_the_fly=False):
    """
    Make png plot of att errors vs time per obsid.

    :param obsid: obsid
    :param att_errors: dictionary with keys including at minimum a coordinate ('dr',
                       'dy' or 'dp' for roll, yaw and pitch) and 'time' (default 'dr')
    :param on_the_fly: default False, if True then ignore param att_errors and derive
                       attitude errors for the requested obsid.
    """

    if coord not in ('dr', 'dp', 'dy'):
        raise ValueError('Coordinate for att error should be dr, dp or dy')

    if on_the_fly:
        att_errors = get_observed_att_errors(obsid, on_the_fly=on_the_fly)
        if att_errors is None:
            return None
    else:
        if att_errors is None:
            raise ValueError('Need to provide att_errors if on_the_fly is False')

    errs = att_errors[coord]
    dates = DateTime(att_errors['time'])

    plt.figure(figsize=(8, 2.5))

    # Skip the first 5 min for observations with duration > 5 min
    dur = dates.secs[-1] - dates.secs[0]
    if dur > 5 * 60:
        ok = dates.secs > dates.secs[0] + 5 * 60
    else:
        ok = np.ones_like(dates.secs, dtype=bool)

    plt.plot(dates.secs[ok] - dates.secs[ok][0], errs[ok],
             '-', lw=2, color='k')

    ylims = plt.ylim()

    if max(ylims) > 100:
        plt.ylim(-max(ylims) - 10, max(ylims) + 10)
    else:
        plt.ylim(-100, 100)

    plt.ylabel(f'{coord} (arcsec)')
    plt.xlabel('Time (sec)')
    plt.grid(ls=':')

    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.25, top=0.95)

    if save:
        outroot = os.path.join(plot_dir, f'observed_{coord}s_{obsid}')
        logger.info(f'Writing plot file {outroot}.png')
        plt.savefig(outroot + '.png')
        plt.close()

    return att_errors


def plot_crs_per_obsid(obsid, plot_dir, crs=None, save=False, on_the_fly=False):
    """
    Make png plot of OBC centroid residuals in each slot. Residuals
    computed using ground attitude solution for science observations
    and OBC attitude solution for ER observations.

    :param crs: dictionary with keys 'ground' and 'obc', dictionary values
                are also dictionaries keyed by slot number containing
                corresponding CentroidResiduals objects
    :param on_the_fly: if True then ignore crs and derive centroid
                       residuals for the requested obsid (default False)
    """

    if on_the_fly:
        crs = get_crs_per_obsid(obsid)
    else:
        if crs is None:
            raise ValueError('Need to provide crs if on_the_fly is False')

    crs_obc = crs['obc']
    crs_grnd = crs['ground']

    cs = {'yag': 'k', 'zag': 'slategray'}

    fig = plt.figure(figsize=(8, 7))

    n = len(crs['slots'])
    legend = False

    for ii, slot in enumerate(crs['slots']):

        plt.subplot(n, 1, ii + 1)

        for coord in ['yag', 'zag']:
            resids_obc = getattr(crs_obc[slot], f'd{coord}s')
            times_obc = getattr(crs_obc[slot], f'{coord}_times')

            if slot not in crs_grnd or crs_grnd[slot] is None:
                resids_ref = resids_obc
                times_ref = times_obc
            else:
                resids_ref = getattr(crs_grnd[slot], f'd{coord}s')
                times_ref = getattr(crs_grnd[slot], f'{coord}_times')

            if len(times_obc) == 0 or len(times_ref) == 0:
                continue

            resids_obc_interp = interpolate(resids_obc, times_obc, times_ref)
            ok = np.abs(resids_obc_interp) <= 5

            plt.plot(times_ref - times_ref[0],
                     np.ma.array(resids_ref, mask=~ok),
                     color=cs[coord],
                     alpha=0.9,
                     label=f'd{coord}')

            if np.sum(~ok) > 0:
                plt.plot(times_ref - times_ref[0],
                         np.ma.array(resids_ref, mask=ok),
                         color='crimson',
                         alpha=0.9)

        plt.grid(ls=':')
        ylims = plt.ylim()
        if max(np.abs(ylims)) < 5:
            plt.ylim(-6, 6)
            plt.yticks([-5, 0, 5], ["-5", "0", "5"])
        else:
            plt.ylim(-12, 12)
            plt.yticks([-10, -5, 0, 5, 10], ["-10", "-5", "0", "5", "10"])
        plt.xlabel('Time (sec)')
        plt.ylabel(f'Slot {slot}\n(arcsec)')
        axs = fig.gca()
        axs.get_yaxis().set_label_coords(-0.065, 0.5)
        handles, labels = axs.get_legend_handles_labels()
        if len(labels) > 0 and not legend:
            plt.legend(loc=1)
            legend = True

    plt.subplots_adjust(left=0.1, right=0.95, hspace=0,
                        bottom=0.08, top=0.98)

    if save:
        outroot = os.path.join(plot_dir, f'crs_time_{obsid}')
        logger.info(f'Writing plot file {outroot}.png')
        plt.savefig(outroot + '.png')
        plt.close()

    return crs


def plot_n_kalman(obsid, plot_dir, save=False):
    """
    Fetch and plot number of Kalman stars as function of time for
    the requested obsid.
    """
    d = events.dwells.filter(obsid=obsid)[0]
    start = d.start
    stop = d.stop
    n_kalman = get_n_kalman(start, stop)

    plt.figure(figsize=(8, 2.5))

    t0 = n_kalman.times[0]

    # The Kalman vals are strings, so these can be out of order on y axis
    # if not handled as ints.
    plot_cxctime(n_kalman.times, n_kalman.vals.astype(int), color='k')
    plot_cxctime([t0, t0 + 1000], [0.5, 0.5], lw=3, color='orange')

    plt.text(DateTime(t0).plotdate, 0.7, "1 ksec")
    plt.ylabel(f'# Kalman stars')
    ylims = plt.ylim()
    plt.ylim(-0.2, ylims[1] + 0.2)
    plt.grid(ls=':')

    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.25, top=0.95)

    if save:
        outroot = os.path.join(plot_dir, f'n_kalman_{obsid}')
        logger.info(f'Writing plot file {outroot}.png')
        plt.savefig(outroot + '.png')
        plt.close()


def plot_crs_visualization(obsid, plot_dir, crs=None, factor=20, save=False, on_the_fly=False):
    """
    Plot visualization of OBC centroid residuals with respect to ground (obc)
    aspect solution for science (ER) observations in the yang/zang plain.

    :param obsid: obsid
    :param crs: dictionary with keys 'ground' and 'obc'. Dictionary values
                are dictionaries keyed by slot number containing corresponding
                CentroidResiduals objects. If ground or obc centroid residuals
                cannot be computed for a given slot, the value is None.
    :param on_the_fly: default False, if True then ignore param crs and calculate
                       centroid residuals for the requested obsid
    """

    # catalog
    cat = get_starcat(obsid)

    # keep only BOT and GUI entries
    ok = (cat['type'] == 'BOT') | (cat['type'] == 'GUI')
    cat = cat[ok]
    cat['idx'] = cat['slot']  # so that the plot is numbered by slot

    # attitude
    att = get_att(obsid)

    # stars
    cols = ['RA_PMCORR', 'DEC_PMCORR', 'MAG_ACA', 'MAG_ACA_ERR',
            'CLASS', 'ASPQ1', 'ASPQ2', 'ASPQ3', 'VAR', 'POS_ERR']
    stars = Table(names=cols)
    for star in cat:
        row = []
        s = get_star(star['id'])
        for col in cols:
            row.append(s[col])
        stars.add_row(row)

    fig = plot_stars(att, cat, stars)
    cs = ['orange', 'forestgreen', 'steelblue', 'maroon', 'gray']

    if on_the_fly:
        crs = get_crs_per_obsid(obsid)
    else:
        if crs is None:
            raise ValueError('Need to provide crs if on_the_fly is False')

    if obsid > 40000:  # ERs
        crs_ref = crs['obc']
    else:
        crs_ref = crs['ground']

    ax = fig.axes[0]

    for slot in cat['slot']:
        ok = cat['slot'] == slot
        yag = cat['yang'][ok]
        zag = cat['zang'][ok]
        yp, zp = yagzag_to_pixels(yag, zag)
        # 1 px -> factor px; 5 arcsec = 5 * factor arcsec
        try:
            """
            Minus sign for y-coord to reflect sign flip in the pixel
            to yag conversion and yag scale going from positive to negative
            """
            yy = yp - crs_ref[slot].dyags * factor
            zz = zp + crs_ref[slot].dzags * factor
            ax.plot(yy, zz, alpha=0.3, marker=',', color=cs[slot - 3])
            ax.plot([-1000, -1020], [2700, 2700], color='k', lw=3)
            circle = plt.Circle((yp, zp), 5 * factor,
                                color='darkorange', fill=False)
            ax.add_artist(circle)
        except Exception:
            pass

    plt.text(-511, 530, "ring radius = 5 arcsec (scaled)", color='darkorange')

    if save:
        outroot = os.path.join(plot_dir, f'crs_vis_{obsid}')
        logger.info(f'Writing plot file {outroot}.png')
        plt.savefig(outroot + '.png')
        plt.close()

    return crs


def plot_observed_metrics(obsids, plot_dir, coord='dr', att_errors=None, factor=20,
                          save=False, on_the_fly=False):
    """
    Generate plots of attitude errors, centroid residuals and number of Kalman stars
    for the centroid dashboard.

    :param obsids: obsid or a set of obsids
    :param att_errors: dictionary containing at minimum the following keys 'time',
                       coordinate ('dr', 'dp' or 'dy), 'crs'
    :param factor: scaling for the centroid residuals visualization
                   in the yang/zang plain
    :param on_the_fly: default False, if True then derive the attitude errors
                       (including centroid residuals) for the requested obsid
    """

    if type(obsids) not in (int, set):
        raise TypeError('obsids should be an int or a set')

    if type(obsids) is int:
        obsids = {obsids}

    for obsid in list(obsids):
        kwargs = {'save': save, 'on_the_fly': on_the_fly}
        errs = plot_att_errors_per_obsid(obsid, coord=coord, att_errors=att_errors,
                                         plot_dir=plot_dir,
                                         **kwargs)
        crs = errs['crs']

        # To prevent another computation of crs
        if on_the_fly:
            kwargs['on_the_fly'] = False

        plot_crs_per_obsid(obsid, plot_dir, crs=crs, **kwargs)
        plot_crs_visualization(obsid, plot_dir, factor=factor, crs=crs, **kwargs)
        plot_n_kalman(obsid, plot_dir, save=save)


def read_metrics_from_file(filename):
    """
    Read in processed guide metrics (dr95, dr50, manvr_angle, ending dr,
    one shot updates, aber corrections) from file
    """
    if os.path.exists(filename):
        logger.info(f'Reading {filename}')
        dat_old = Table.read(filename, format='ascii.ecsv', guess=False)
        processed_obsids = set(dat_old['obsid'])
    else:
        logger.info(f'File {filename} does not exist')
        dat_old = None
        processed_obsids = set()

    return dat_old, processed_obsids


class NoObsidError(ValueError):
    pass


def update_observed_metrics(obsid=None, start=None, stop=None, data_root=None, force=False,
                            factor=20, make_plots=False, save=False):
    """
    Update the ``GUIDE_METRICS_OBSID`` and ``GUIDE_METRICS_SLOT`` tables
    (ascii ECSV format) in place to reflect information about observed
    guide metrics: dr95, dr50, manvr angle, ending roll error, one shot
    updates, aberration corrections; and mean yag/zag, mean dyag/dzag
    centroid errors.

    :param factor: scaling for the centroid residual visualization in
                   the yag/zag plane (default 20)
    :param make_plots: if True then generate plots for centroid dashboard
                       (default False)
    :param save: if True then save the plots (default False)
    """

    if obsid is None:
        # Default is between NOW and (NOW - NDAYS) days
        start = DateTime(start) - (NDAYS if start is None else 0)
        stop = DateTime(stop)
        # Get obsids, both science and ERs
        obsids = [evt.obsid for evt in events.obsids.filter(start, stop)]
    else:
        obsids = [np.int(obsid)]

    if data_root is None:
        data_root = CD_ROOT

    obsid_metrics_file = os.path.join(data_root, GUIDE_METRICS_OBSID)
    slot_metrics_file = os.path.join(data_root, GUIDE_METRICS_SLOT)
    # Read in existing files if they exists and make a set of already-processed obsids
    dat_obsid_old, processed_obsids = read_metrics_from_file(obsid_metrics_file)
    dat_slot_old, tmp = read_metrics_from_file(slot_metrics_file)

    rows_obsid = []
    rows_slots = []
    for obsid in obsids:

        logger.info(f'Obsid={obsid}')
        obs_dir = get_cd_dir(obsid, data_root)

        if obsid in processed_obsids:
            if not force:
                logger.info(f'Skipping obsid {obsid}: already processed')
                continue
            else:
                if obsid in dat_obsid_old['obsid']:
                    idx = list(dat_obsid_old['obsid']).index(obsid)
                    dat_obsid_old.remove_row(idx)

                if obsid in dat_slot_old['obsid']:
                    ok = dat_slot_old['obsid'] == obsid
                    dat_slot_old.remove_rows(ok)

        try:
            metrics_obsid, metrics_slot = get_observed_metrics(obsid)

            if not metrics_obsid['dwell']:
                logger.info(f'Skipping obsid {obsid}: not a dwell?')
                info = 'Not a dwell'
                make_special_case_html(metrics_obsid, obs_dir, info=info)
                continue

            if metrics_obsid['att_flag'] > 1:
                logger.info(f'Skipping obsid {obsid}: problem matching obc/ground times')
                info = 'Problem matching obc/ground att times'
                make_special_case_html(metrics_obsid, obs_dir, info=info)
                continue

            if obsid < 40000 and metrics_obsid['att_flag'] == 1:
                logger.info(f'Skipping science obsid {obsid}: no ground aspect solution')
                info = 'No ground aspect solution for science obsid'
                make_special_case_html(metrics_obsid, obs_dir, info=info)
                continue

            if make_plots:
                kwargs = {'factor': factor, 'save': save}
                plot_observed_metrics(obsid,
                                      plot_dir=obs_dir,
                                      coord='dr',
                                      att_errors=metrics_obsid['att_errors'],
                                      **kwargs)

        except NoObsidError:  # not yet in archive
            logger.info(f'Skipping obsid {obsid}: not yet in archive')
            continue
        except Exception as err:
            logger.info(f'ERROR: {err}')
            continue

        # Process entries for 'per obsid' metrics

        keys_obsid = ('obsid', 'mean_date',
                      'att_flag', 'dr50', 'dr95',
                      'aber_y', 'aber_z', 'aber_flag',
                      'one_shot_pitch', 'one_shot_yaw',
                      'one_shot', 'one_shot_aber_corrected',
                      'manvr_angle', 'ending_roll_err',
                      'obsid_preceding', 'obsid_next')

        row_obsid = {k: metrics_obsid[k] for k in keys_obsid}
        rows_obsid.append(row_obsid)

        # Process entries for 'per slot' metrics
        row_slots = []
        for slot in range(8):
            out = {}
            slot_data = metrics_slot['slots'][slot]
            if bool(slot_data):
                out['obsid'] = obsid
                out['slot'] = slot
                out['mean_date'] = metrics_obsid['mean_date']
                keys_slot = ('id', 'type', 'mag', 'yang', 'zang',
                             'median_mag', 'median_dy', 'median_dz')
                out.update({k: slot_data[k] for k in keys_slot})
                # Needed to build html
                row_slots.append(out)
                # Needed to update 'per slot' data file
                rows_slots.append(out)

        # Build html page for this obsid
        make_html(row_obsid, row_slots, obs_dir)

        # Update the 'per_obsid' table
        if rows_obsid:
            sort_cols = ['mean_date']
            update_data_table(rows_obsid, dat_obsid_old, obsid_metrics_file, sort_cols)

        # Update the 'per_slot' table
        if rows_slots:
            sort_cols = ['mean_date', 'slot']
            update_data_table(rows_slots, dat_slot_old, slot_metrics_file, sort_cols)


def update_data_table(rows, dat_old, metrics_file, sort_cols):
    """
    Update data files in place

    :param rows: list of dicts with newly processed'per obsid' or 'per slot' entries
    :param dat_old: old entries
    :param metrics_file: file with metrics to be updated in place
    :param sort_cols: columns used for sorting
    """
    dat = Table(rows=rows, names=rows[0])
    if dat_old is not None:
        dat = vstack([dat_old, dat])

    logger.info(f'Writing {metrics_file}')

    dat.sort(sort_cols)
    dat.write(metrics_file, format='ascii.ecsv', overwrite=True)


def make_html(row_obsid, rows_slot, obs_dir):
    """
    Build html page

    :param metrics_obsid: dict with `per obsid` metrics
    :param metrics_slot: list of dicts with `per slot` metrics
    """

    # Build 'per obsid' part of the webpage
    obsid = row_obsid['obsid']
    obsid_preceding = row_obsid['obsid_preceding']
    obsid_next = row_obsid['obsid_next']
    ending_roll_err = row_obsid['ending_roll_err']
    one_shot_pitch = row_obsid['one_shot_pitch']
    one_shot_yaw = row_obsid['one_shot_yaw']
    one_shot = row_obsid['one_shot']
    one_shot_aber_corrected = row_obsid['one_shot_aber_corrected']
    manvr_angle = row_obsid['manvr_angle']
    dr50 = row_obsid['dr50']
    dr95 = row_obsid['dr95']
    aber_y = row_obsid['aber_y']
    aber_z = row_obsid['aber_z']
    mean_date = row_obsid['mean_date']

    cd_preceding_root = os.path.join('../../', np.str(obsid_preceding)[:2], f"{obsid_preceding}")
    cd_next_root = os.path.join('../../', np.str(obsid_next)[:2], f"{obsid_next}")

    star_path_root = os.path.join(MICA_REPORTS, np.str(obsid)[:2], f"{obsid}")

    if obsid_preceding == -9999:
        preceding_obsid_link = ""
    else:
        preceding_obsid_link = f"<a href='{cd_preceding_root}/index.html'>Previous</a>"

    if obsid_next == -9999:
        next_obsid_link = ""
    else:
        next_obsid_link = f"<a href='{cd_next_root}/index.html'>Next</a>"

    string = f"""
<html>
<head>
<title>Obsid {obsid}</title>

<script type="text/javascript" src="/aspect/overlib.js"></script>
<style type="text/css">
  body {{ min-width:900px;
         margin-left: 20px;
         background-color: lightblue;
}}
  div#fullwidth {{
    clear: both;
    position: static;
    width: 1200px;
    left: 0px;
  }}
  div#leftmain {{
    float: left;
    width: 520px;
    left: 0px;
    top: 0px;
  }}
  div#leftsmall {{
    float: left;
    width: 500px;
  }}
  div.border {{
    background-color: white;
    padding-left: 10px;
    padding-bottom: 5px;
    border-style: solid;
    border-width: 5px;
    border-radius: 20px;
  }}
  div#rightlarge {{
    float: right;
    width: 415px;
  }}
  div#rightsmall {{
    float: right;
    width: 635px;
    left: 650px;
    top: 0px;
    padding-right: 5px;
  }}

button {{
    background:none!important;
    border:none;
    padding:0!important;
    color:#069;
    text-decoration:underline;
    cursor:pointer;
}}
</style>
</head>

<body>

<div id="fullwidth">

    <div id="leftmain">

        <div class="border" style="padding: 10px; margin-bottom: .5cm">
            <span id="label" style="font-size:150%; font-weight:bold;">
                {preceding_obsid_link}
                -- Obsid {obsid} --
                {next_obsid_link}
            </span>
        </div>

        <div class="border" style="float: left; margin-bottom: .5cm;">

            <div id="leftsmall">
                <h2>Performance Details</h2>

                <pre>
OBSID <a href="{MICA_PORTAL}{obsid}">{obsid}</a>         Mean date: {DateTime(mean_date).caldate}

  One shot Pitch: {one_shot_pitch:7.2f} arcsec    Aber-y = {aber_y:6.2f} arcsec
  One shot Yaw:   {one_shot_yaw:7.2f} arcsec    Aber-z = {aber_z:6.2f} arcsec
  One shot = {one_shot:.2f} arcsec
  One shot with aberration correction = {one_shot_aber_corrected:.2f} arcsec
  Manvr angle = {manvr_angle:.2f} deg

  dr50 = {dr50:.2f} arcsec
  dr95 = {dr95:.2f} arcsec

Preceding Observation

  OBSID <a href="{MICA_PORTAL}{obsid_preceding}">{obsid_preceding}</a>
  Ending roll error = {ending_roll_err:.2f} arcsec

Next Observation

  OBSID <a href="{MICA_PORTAL}{obsid_next}">{obsid_next}</a>
                </pre>
            </div>
        </div>

        <div class="border" style="float: left; margin-bottom: .5cm;">
            <div id="leftsmall">
                <h2>Catalog and visualization of centroid residuals</h2>
<table border=1 style="font-size:11px" align="center">
<tr>
<th align='right'>SLOT</th>
<th align='right'>ID</th>
<th align='right'>TYPE</th>
<th align='right'>MAG</th>
<th align='right'>YANG</th>
<th align='right'>ZANG</th>
<th align='right'>MEDIAN MAG</th>
<th align='right'>MEDIAN DY</th>
<th align='right'>MEDIAN DZ</th>
</tr>
"""

    # Build 'per slot' part of the webpage

    t_slot = Table(rows_slot)
    t_slot.sort('slot')

    for row in t_slot:
        string += f"""<tr>
<td align='right'>{row['slot']}</td>
"""
        if row['id'] < 100:
            id_ = ""
        else:
            id_ = (f"<a href='{star_path_root}/star_{row['id']}.html' "
                   + f"'style='text-decoration: none;'>{row['id']}</a>")

        string += f"""<td align='right'>{id_}</td>
<td align='right'>{row['type']}</td>
<td align='right'>{row['mag']}</td>
<td align='right'>{row['yang']}</td>
<td align='right'>{row['zang']}</td>
<td align='right'>{row['median_mag']:.3f}</td>
<td align='right'>{row['median_dy']:.2f}</td>
<td align='right'>{row['median_dz']:.2f}</td>
</tr>
"""
    string += f"""</table>

<img src="crs_vis_{obsid}.png" width="490" height="490">
            </div>
        </div>

    </div>

    <div id="rightsmall" class="border" style="padding-top: 10px; margin-bottom: .5cm">
        <img src="observed_drs_{obsid}.png" width="635" height="198">
    </div>
    <div id="rightsmall" class="border" style="padding-top: 10px; margin-bottom: .5cm">
        <img src="crs_time_{obsid}.png" width="635" height="556">
    </div>
    <div id="rightsmall" class="border" style="padding-top: 10px; margin-bottom: .5cm">
        <img src="n_kalman_{obsid}.png" width="635" height="198">
    </div>

</div>

</body>
</html>
"""

    logger.info(f'Writing index.html for obsid {obsid}')
    with open(f"{obs_dir}/index.html", "w") as outfile:
        outfile.write(string)


def make_special_case_html(metrics_obsid, obs_dir, info=''):
    """
    Build html page

    :param metrics_obsid: dict with limited `per obsid` metrics
    """

    # Build 'per obsid' part of the webpage
    obsid = metrics_obsid['obsid']
    obsid_preceding = metrics_obsid['obsid_preceding']
    obsid_next = metrics_obsid['obsid_next']

    cd_preceding_root = os.path.join('../../', np.str(obsid_preceding)[:2], f"{obsid_preceding}")
    cd_next_root = os.path.join('../../', np.str(obsid_next)[:2], f"{obsid_next}")

    if obsid_preceding == -9999:
        preceding_obsid_link = ""
    else:
        preceding_obsid_link = f"<a href='{cd_preceding_root}/index.html'>Previous</a>"

    if obsid_next == -9999:
        next_obsid_link = ""
    else:
        next_obsid_link = f"<a href='{cd_next_root}/index.html'>Next</a>"

    string = f"""
<html>
<head>
<title>Obsid {obsid}</title>

<script type="text/javascript" src="/aspect/overlib.js"></script>
<style type="text/css">
  body {{ min-width:900px;
         margin-left: 20px;
         background-color: lightblue;
}}
  div#fullwidth {{
    clear: both;
    position: static;
    width: 1200px;
    left: 0px;
  }}
  div#leftmain {{
    float: left;
    width: 520px;
    left: 0px;
    top: 0px;
  }}
  div#leftsmall {{
    float: left;
    width: 500px;
  }}
  div.border {{
    background-color: white;
    padding-left: 10px;
    padding-bottom: 5px;
    border-style: solid;
    border-width: 5px;
    border-radius: 20px;
  }}
  div#rightlarge {{
    float: right;
    width: 415px;
  }}
  div#rightsmall {{
    float: right;
    width: 635px;
    left: 650px;
    top: 0px;
    padding-right: 5px;
  }}

button {{
    background:none!important;
    border:none;
    padding:0!important;
    color:#069;
    text-decoration:underline;
    cursor:pointer;
}}
</style>
</head>

<body>

<div id="fullwidth">

    <div id="leftmain">

        <div class="border" style="padding: 10px; margin-bottom: .5cm">
            <span id="label" style="font-size:150%; font-weight:bold;">
                {preceding_obsid_link}
                -- Obsid {obsid} --
                {next_obsid_link}
            </span>
        </div>

        <div class="border" style="float: left; margin-bottom: .5cm;">

            <div id="leftsmall">
                <h2>Performance Details</h2>

                <pre>
OBSID <a href="{MICA_PORTAL}{obsid}">{obsid}</a>

  {info}

Preceding Observation

  OBSID <a href="{MICA_PORTAL}{obsid_preceding}">{obsid_preceding}</a>

Next Observation

  OBSID <a href="{MICA_PORTAL}{obsid_next}">{obsid_next}</a>
                </pre>
            </div>
        </div>

    </div>
</div>

</body>
</html>
"""

    logger.info(f'Writing index.html for obsid {obsid}')
    with open(f"{obs_dir}/index.html", "w") as outfile:
        outfile.write(string)


def main():
    opt = get_opt()
    logger.info('Centroid dashboard, started')
    update_observed_metrics(obsid=opt.obsid, start=opt.start, stop=opt.stop,
                            force=opt.force, data_root=opt.data_root,
                            make_plots=opt.make_plots, save=opt.save)
    logger.info('Centroid dashboard, ended')


if __name__ == '__main__':
    main()
