from __future__ import division
import os
import re
import pyfits
from glob import glob
import numpy as np
from Ska.Table import read_table
from scipy.signal import medfilt as medfilt

import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

# borrowed from telem_archive
import csv
import gzip


def parse_obspar(file):
    convert = {'i': int,
               'r': float,
               's': str}
    try:
        lines = gzip.open(file).readlines()
    except IOError:
        lines = open(file).readlines()
    obs_read = csv.DictReader(lines,
                              fieldnames=('name', 'type', 'hidden', 'value',
                                          'def1', 'def2', 'descr'),
                              dialect='excel')

    for row in obs_read:
        row['value'] = convert[row['type']](row['value'])
        row['name'] = row['name'].replace('-', '_')
        yield row

    return


def get_obspar(obsparfile):
    """Get the obspar for obsid starting at tstart.  Return as a dict."""

    obspar = {'num_ccd_on': 0}
    for row in parse_obspar(obsparfile):
        obspar.update({row['name']: row['value']})
        if re.match(r'^ccd[is]\d_on$', row['name']) and row['value'] == 'Y':
            obspar['num_ccd_on'] += 1

    return obspar


r2a = 206264.81
d2r = 0.017453293
M0 = 10.32
C0 = 5263.0

enc_rad_1 = 0.68
enc_rad_2 = 0.99
dyz_big_lim = 0.7
frac_dyz_lim = 0.05

data_cols = ['qual', 'dy', 'dz', 'mag', 'time',
             'ang_y_sm', 'ang_y', 'ang_z_sm', 'ang_z']


def rms(data, median):
    return np.sqrt(np.mean((data - median) ** 2))


def frac_bad(data, median, limit):
    big_offset = (data - median) > limit
    frac_bad = len(np.flatnonzero(big_offset)) / len(data)
    return frac_bad


class Obi(object):
    def __init__(self, obsparfile, obsdir):
        self.obspar = get_obspar(obsparfile)
        self.obsdir = obsdir
        self.find_aspect_intervals()
        self.process_aspect_intervals()
        self.slot = dict()
        self.concat_slot_data()
        self.check_intervals()
        self.agg_slot_data()

        #self.plots = [self.plot_slot(slot) for slot in range(0, 8)]


    def find_aspect_intervals(self):
        obsdir = self.obsdir
        asol_files = sorted(glob(os.path.join(obsdir, 'pcad*asol*')))
        self.aiids = []
        for file in asol_files:
            hdulist = pyfits.open(file)
            header = hdulist[1].header
            # skip files that aren't in the obspar range
            if ((self.obspar['tstart'] >= header['TSTOP'])
                or (self.obspar['tstop'] <= header['TSTART'])):
                continue
            aiid_match = re.search('(pcadf\d+N\d{3})', file)
            if aiid_match:
                self.aiids.append(aiid_match.group(0))

    def process_aspect_intervals(self):
        self.aspect_intervals = []
        for aiid in self.aiids:
            self.aspect_intervals.append(AspectInterval(aiid, self.obsdir))

    def concat_slot_data(self):
        all_slot = self.slot
        for slot in range(0, 8):
            cslot = dict()
            for d in data_cols:
                slotval = [i.deltas[slot][d] for i in self.aspect_intervals]
                cslot[d] = np.concatenate(slotval)
            all_slot[slot] = cslot


    def agg_slot_data(self):
        all_slot = self.slot
        for slot_id in all_slot:
            slot = all_slot[slot_id]
            # these should only be calculated over good data, right?
            for axdir in ['dy', 'dz', 'mag']:
                data = slot[axdir][slot['qual'] == 0]
                smean = np.mean(data)
                slot['%s_mean' % axdir] = smean
                med = np.median(data)
                slot['%s_med' % axdir] = med
                srms = rms(data, med)
                slot['%s_rms' % axdir] = srms
                bad_frac = frac_bad(data, med, dyz_big_lim)
                slot['frac_%s_big' % axdir] = bad_frac
            slot['rad_off'] = np.sqrt(slot['dy_med']**2
                                      + slot['dz_med']**2)
            slot['n_pts'] = len(slot['dy'])
            if slot['type'] == 'fid':
                slot['mean_y'] = np.mean(slot['ang_y_sm'])
                slot['mean_z'] = np.mean(slot['ang_z_sm'])
            else:
                slot['mean_y'] = np.mean(slot['ang_y'])
                slot['mean_z'] = np.mean(slot['ang_z'])



    def check_intervals(self):
        all_slot = self.slot
        # check for weirdness across intervals
        slot_list = {}
        slot_status = {}
        for t in ('gsprop', 'fidprop'):
            slot_list[t] = getattr(self.aspect_intervals[0], t).slot
            slot_status[t] = getattr(self.aspect_intervals[0], t).id_status
            for ai in self.aspect_intervals[1:]:
                if len(slot_list[t]) != len(getattr(ai, t).slot):
                    raise ValueError(
                        "differing %s slots across aspect intervals" % t)
                if (len(slot_list[t]) == len(getattr(ai, t).slot)
                    & (not np.all(slot_list[t] == getattr(ai, t).slot))):
                    raise ValueError(
                        "differing %s slots across aspect intervals" % t)
                if (len(slot_list[t]) == len(getattr(ai, t).slot)
                    & (not np.all(slot_status[t] == getattr(ai, t).id_status))):
                    raise ValueError(
                        "differing %s status across aspect intervals" % t)
        for gs in slot_list['gsprop']:
            all_slot[gs]['type'] = 'guide'
        for fs in slot_list['fidprop']:
            all_slot[fs]['type'] = 'fid'

        self.guide_slots = slot_list['gsprop']
        self.fid_slots = slot_list['fidprop']
        

    def plot_slot(self, slot):
        y = None
        z = None
        xy_range = None
        (qual, dy, dz, mag, time) = [self.slot[slot][x] for x in data_cols]
        ai_starts = [interv.deltas[slot]['time'][0]
                     for interv in self.aspect_intervals]
        time0 = time[0]
        dy0 = np.median(dy)
        dz0 = np.median(dz)

        fig = plt.figure(num=slot, figsize=(12, 10))
        #fid_plot = np.abs(np.max(y) - np.min(y)) > 1e-6
        fid_plot = slot < 3
        ok = qual == 0
        bad = qual != 0

        if fid_plot and not xy_range:
            xy_range = 0.1
            circ_rad = 0.05
        else:
            xy_range = 1.0
            circ_rad = 0.7
        #    xmarg = [10, 10]
        #    ystyl = 8
        #else:
        #    xmarg = [10, 3]
        #    ystyl = 1

        if time0 is not None:
            plottime = (time - time0) / 1000.
        else:
            plottime = time / 1000.

        labelfontsize = 10

        ayz = fig.add_axes([.05, .25, .20, .20])
        ayz.plot(dy[ok], dz[ok], 'g.')
        ayz.plot(dy[bad], dz[bad], 'r.')
        ayz.set_aspect('equal', 'datalim')
        plt.setp(ayz.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(ayz.get_xticklabels(), fontsize=labelfontsize)
        circle = plt.Circle((dy0, dz0), radius=circ_rad)
        ayz.add_patch(circle)
        ayz.set_xlabel('Y offset (arcsec)')
        ayz.set_ylabel('Z offset (arcsec)')

        ayzf = fig.add_axes([.05, .7, .20, .20])
        ayzf.plot(dy[ok], dz[ok], 'g.')
        ayzf.plot(dy[bad], dz[bad], 'r.')
        ayzf.set_aspect('equal')
        plt.setp(ayzf.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(ayzf.get_xticklabels(), fontsize=labelfontsize)
        ayzf.set_xlabel('Y offset (arcsec)')
        ayzf.set_ylabel('Z offset (arcsec)')
        if xy_range is not None:
            ayzf.set_xlim([dy0 - xy_range, dy0 + xy_range])
            ayzf.set_ylim([dz0 - xy_range, dz0 + xy_range])

        ay = fig.add_axes([.32, .7, .62, .25])
        ay.plot(plottime[ok], dy[ok], color='green')
        ay.plot(plottime[bad], dy[bad], 'r.')
        plt.setp(ay.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(ay.get_xticklabels(), visible=False)
        for t in ai_starts:
            s_t = (t - time0) / 1000.
            ay.plot([s_t, s_t], ay.get_ylim(), color='blue',
                    linestyle='dashed')
        ay.grid()

        if y is not None:
            ay2 = ay.twinx()
            ay2.plot(plottime, y)
            ay2y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            ay2.yaxis.set_major_formatter(ay2y_formatter)
            plt.setp(ay2.get_yticklabels(), fontsize=labelfontsize,
                     color='blue')

        az = fig.add_axes([.32, .4, .62, .25], sharex=ay)
        az.plot(plottime[ok], dz[ok], color='green')
        az.plot(plottime[bad], dz[bad], 'r.')
        plt.setp(az.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(az.get_xticklabels(), visible=False)
        for t in ai_starts:
            s_t = (t - time0) / 1000.
            az.plot([s_t, s_t], az.get_ylim(), color='blue',
                    linestyle='dashed')
        az.grid()

        if z is not None:
            az2 = az.twinx()
            az2.plot(plottime, z)
            az2y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            az2.yaxis.set_major_formatter(az2y_formatter)
            plt.setp(az2.get_yticklabels(), fontsize=labelfontsize,
                     color='blue')

        am = fig.add_axes([.32, .1, .62, .25], sharex=ay)
        plt.setp(am.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(am.get_xticklabels(), fontsize=labelfontsize)
        am.plot(plottime[ok], mag[ok], color='green')
        am.plot(plottime[bad], mag[bad], 'r.')
        am.set_ylim(am.get_ylim()[::-1])
        for t in ai_starts:
            s_t = (t - time0) / 1000.
            am.plot([s_t, s_t], am.get_ylim(), color='blue',
                    linestyle='dashed')
        am.grid()
        return fig


class AspectInterval(object):
    def __init__(self, aiid, obsdir, opt=None):
        self.aiid = aiid
        self.obsdir = obsdir
        self.opt = opt if opt else {'obc': None, 'alg': 8, 'noacal': None}
        self.read_in_data()
        self.deltas = {}
        self.calc_fid_deltas()
        self.calc_guide_deltas()

    def get_prop(self, propname, propstring):
        datadir = self.obsdir
        print 'Reading %s stars' % propname
        gsfile = glob(os.path.join(
                datadir, "%s_%s1.fits*" % (self.aiid, propstring)))[0]
        prop_all = read_table(gsfile)
        good = prop_all['id_status'] == 'GOOD      '
        if len(np.flatnonzero(good)) == 0:
            raise ValueError("No good %s stars" % propname)
        prop = prop_all[good]

        info = []
        hdulist = pyfits.open(os.path.join(datadir, gsfile))
        header = hdulist[1].header
        proptable = read_table(os.path.join(datadir, gsfile))
        for gs in proptable:
            info.append(dict(slot=gs['slot'],
                             id_status=gs['id_status'],
                             tstart=header['TSTART'],
                             tstop=header['TSTOP']))
        return (prop, info, header)

    def read_in_data(self):
        aiid = self.aiid
        datadir = self.obsdir
        opt = self.opt

        (self.gsprop, self.gspr_info, self.h_gspr) \
            = self.get_prop('guide', 'gspr')
        (self.fidprop, self.fidpr_info, self.h_fidpr) \
            = self.get_prop('fid', 'fidpr')

        print 'Reading aspect solution'
        if opt['obc']:
            asol = read_table(glob(
                    os.path.join(datadir, "%s_osol1.fits*" % aiid))[0])
        else:
            asol = read_table(glob(
                    os.path.join(datadir, "%s_asol1.fits*" % aiid))[0])
        self.asol = asol

        print 'Reading aspect quality'
        self.aqual = read_table(glob(
                os.path.join(datadir, "%s_aqual1.fits*" % aiid))[0])

        #if opt['noacal']:
        #    aca_misalign = np.array([[1.0,0,0], [0,1,0],[0,0,1]])
        #else:
        print 'Reading ACA and FTS align file'
        acal = read_table(glob(
                os.path.join(datadir, "%s_acal1.fits*" % aiid))[0])
        self.aca_misalign = acal['aca_misalign'].reshape(3, 3)
        self.fts_misalign = acal['fts_misalign'].reshape(3, 3)
        self.acal = acal

        print 'Reading Centroids'
        cen = read_table(glob(
                os.path.join(datadir, "%s_acen1.fits*" % aiid))[0])
        # do we want the other algs?
        self.cen = cen[(cen['alg'] == opt['alg'])
                  & (cen['time'] >= asol[0]['time'])
                  & (cen['time'] <= asol[-1]['time'])]
        #          & (cen['status'] == 0)
        self.cenhdulist = pyfits.open(glob(
                os.path.join(datadir, "%s_acen1.fits*" % aiid))[0])
        self.integ_time = self.cenhdulist[1].header['INTGTIME']

    def calc_fid_deltas(self):
        asol = self.asol
        h_fidpr = self.h_fidpr
        aca_misalign = self.aca_misalign
        fts_misalign = self.fts_misalign
        fidprop = self.fidprop
        cen = self.cen
        fidpr_info = self.fidpr_info
        integ_time = self.integ_time

        #mm2a = 20.0
        #dy0 = np.median(asol['dy']) * mm2a
        #dz0 = np.median(asol['dz']) * mm2a
        #dt0 = np.median(asol['dtheta']) * 3600
        #dr = 2.0

        print 'Calculating fid solution quality'

        lsi0_stt = [h_fidpr['LSI0STT%d' % x] for x in [1, 2, 3]]
        stt0_stf = [h_fidpr['STT0STF%d' % x] for x in [1, 2, 3]]
        rrc0_fc_x = h_fidpr['RRC0FCX']

        M = np.dot(aca_misalign, fts_misalign)

        #fid_dy_rms = np.zeros(len(fidprop))
        #fid_dz_rms = np.zeros(len(fidprop))
        #fid_dy_med = np.zeros(len(fidprop))
        #fid_dz_med = np.zeros(len(fidprop))

        rot_x = np.zeros([3, 3])
        rot_x[0, 0] = 1
        for fid in fidprop:
            print "Processing fid %s in slot %d " % (
                fid['id_string'], fid['slot'])
            #slot = fid['slot']
            p_lsi = fid['p_lsi']
            ok = cen['slot'] == fid['slot']
            ceni = cen[ok]
            #n_ceni = len(ceni)
            #dyag = np.zeros(n_ceni)
            #dzag = np.zeros(n_ceni)
            #y = ceni['ang_y_sm'] * 3600
            #z = ceni['ang_z_sm'] * 3600
            dy = np.interp(ceni['time'], asol['time'], asol['dy'])
            dz = np.interp(ceni['time'], asol['time'], asol['dz'])
            dtheta = (np.interp(ceni['time'], asol['time'], asol['dtheta'])
                      * d2r)
            #dtheta[0] = 1.8455539e-05
            #dy[0] = 0.46696303
            #dz[0] = 0.68532118
            p_stf = p_lsi + lsi0_stt + stt0_stf
            for j in range(0, len(ceni)):
                s_th = np.sin(dtheta[j])
                c_th = np.cos(dtheta[j])
                rot_x[1, 1] = c_th
                rot_x[2, 1] = s_th
                rot_x[1, 2] = -s_th
                rot_x[2, 2] = c_th

                p_fc = np.dot(rot_x.transpose(), p_stf)
                p_fc = p_fc + [0., dy[j], dz[j]]
                d_fc = p_fc
                d_fc[0] = d_fc[0] - rrc0_fc_x
                d_fc = -d_fc
                d_aca = np.dot(M, d_fc)
                yag = np.arctan2(d_aca[1], d_aca[0]) * r2a
                zag = np.arctan2(d_aca[2], d_aca[0]) * r2a
                dy[j] = ceni[j]['ang_y_sm'] * 3600 - yag
                dz[j] = ceni[j]['ang_z_sm'] * 3600 - zag

            qual = ceni['status'].copy()
            slot_fidpr = [pr for pr in fidpr_info if pr['slot'] == fid['slot']]
            if not slot_fidpr:
                raise ValueError("No FIDPR info found for slot %d"
                                 % fid['slot'])
            for pr in slot_fidpr:
                if pr['id_status'] != 'GOOD      ':
                    bad = ((ceni['time'] >= pr['tstart'])
                           & (ceni['time'] <= pr['tstop']))
                    n_bad = len(np.flatnonzero(bad))
                    if n_bad:
                        qual[bad] = 1
                    else:
                        err = "Bad fidpr interval not contained in ceni range"
                        raise ValueError(err)

            mag = medfilt(M0 - 2.5
                          * np.log10(ceni['counts'][ceni['counts'] > 10.0]
                                     / integ_time
                                     / C0), 3)

            self.deltas[fid['slot']] = dict(time=ceni['time'],
                                            dy=dy,
                                            dz=dz,
                                            mag=mag,
                                            qual=qual,
                                            ang_y_sm=ceni['ang_y_sm'],
                                            ang_z_sm=ceni['ang_z_sm'],
                                            ang_y=ceni['ang_y'],
                                            ang_z=ceni['ang_z'],
                                            )

    def calc_guide_deltas(self):
        from Quaternion import Quat

        asol = self.asol
        cen = self.cen
        gsprop = self.gsprop
        aca_misalign = self.aca_misalign
        integ_time = self.integ_time

        print 'Interpolating quaternions'
        q_att = np.array([np.interp(cen['time'],
                                    asol['time'],
                                    asol['q_att'][:, ax])
                          for ax in range(0, 4)]).transpose()
        for star in gsprop:
            print 'Processing star in slot %(slot)d' % star
            ok = cen['slot'] == star['slot']
            ceni = cen[ok]
            q_atti = q_att[ok]
            print 'Found %d centroids ' % len(ceni)
            # use ang_y or ang_y_sm?
            dy = np.zeros_like(ceni['ang_y'])
            dz = np.zeros_like(ceni['ang_z'])
            for i in range(0, len(ceni)):
                att = Quat(q_atti[i])
                T = att.transform
                # IDL was:
                #d_aca = aca_misalign ## T ## gsprop[i].pos_eci
                #  ;         MNC2ACA     *  ECI2MNC  *    ECI
                d_aca = np.dot(np.dot(aca_misalign,
                                      T.transpose()), star['pos_eci'])
                yag = np.arctan2(d_aca[1], d_aca[0]) * r2a
                zag = np.arctan2(d_aca[2], d_aca[0]) * r2a
                dy[i] = ceni[i]['ang_y'] * 3600 - yag
                dz[i] = ceni[i]['ang_z'] * 3600 - zag
            mag = medfilt(M0 - 2.5
                          * np.log10(ceni['counts'][ceni['counts'] > 10.0]
                                     / integ_time
                                     / C0), 3)

            qual = ceni['status'].copy()
            slot_fidpr = [pr for pr in self.gspr_info
                          if pr['slot'] == star['slot']]
            if not slot_fidpr:
                err = "No GSPR info found for slot %d" % star['slot']
                raise ValueError(err)
            for pr in slot_fidpr:
                if pr['id_status'] != 'GOOD      ':
                    bad = ((ceni['time'] >= pr['tstart'])
                           & (ceni['time'] <= pr['tstop']))
                    n_bad = len(np.flatnonzero(bad))
                    if n_bad:
                        qual[bad] = 1
                    else:
                        err = "Bad gspr interval not contained in ceni range"
                        raise ValueError(err)

            self.deltas[star['slot']]= dict(dy=dy,
                                            dz=dz,
                                            time=ceni['time'],
                                            mag=mag,
                                            qual=qual,
                                            ang_y_sm=ceni['ang_y_sm'],
                                            ang_z_sm=ceni['ang_z_sm'],
                                            ang_y=ceni['ang_y'],
                                            ang_z=ceni['ang_z'],
                                            )
    def calc_sim_offset(self):
        mm2a = 20.0
        max_d_dyz = 0.2
        max_drift = 3.0
        max_abs_sim = 10
        abs_sim_dy0 = 8.0
        abs_sim_dz0 = 8.0
        n_med = 21
        medf_dy = medfilt(self.asol.dy * mm2a, kernel_size=n_med)
        medf_dz = medfilt(self.asol.dz * mm2a, kernel_size=n_med)
        d_dy = abs(medf_dy[1:] - medf_dy[:-1])
        d_dz = abs(medf_dz[1:] - medf_dz[:-1])
        self.sim['time'] = self.asol['time'][:-1]
        self.sim['d_dy'] = d_dy
        self.sim['d_dz'] = d_dz


        
#import unittest
#from nose.tools import eq_ as 
import operator as op
from functools import partial

high_limits = dict(fid=dict(dy_rms=0.05,
                            dz_rms=0.05,
                            dy_med=0.4,
                            dz_med=0.4,
                            rad_off=0.4),
                   guide=dict(dy_rms=0.4,
                              dz_rms=0.4,
                              rad_off=1.5,
                              frac_dy_big=0.05,
                              frac_dz_big=0.05))
low_limits = dict(fid=dict(),
                  guide=dict(dy_rms=0.25,
                             dz_rms=0.25))


med_cols = ['dy_med', 'dz_med', 'rad_off']
rms_cols = ['dy_rms', 'dz_rms', 'dy_rms_low', 'dz_rms_low']
frac_cols = ['frac_dy_big', 'frac_dz_big']



class ObiTest(object):
    def __init__(self, obi=None):
        self.obi = obi
        self.tests = []


    def slot_checks(self):
        obi_slot_check = {}
        for slot_id in self.obi.slot:
            slot_check = {}
            slot = self.obi.slot[slot_id]
            hlim = high_limits[slot['type']]
            llim = low_limits[slot['type']]
            for check in hlim:
                stest = dict(name='%s' % (check),
                             success=slot[check] < hlim[check],
                             have=slot[check],
                             limit=hlim[check],
                             param=check)
                slot_check[stest['name']] = stest
            for check in llim:
                stest = dict(name='%s_low' % (check),
                             success=slot[check] < llim[check],
                             have=slot[check],
                             limit=llim[check],
                             param=check)
                slot_check[stest['name']] = stest
            obi_slot_check[slot_id] = slot_check
        self.checks = obi_slot_check

    def slot_table(self):
        checks = self.checks
        

#    def set_disposition(self, disp):
#        weight = {'OK' : 1,
#                  'ReproNoDD' : 2
#                  'Hold' : 3}
#        if disp in weight:
#            if (not overall_status in self 
#                or weight[self.overall_status] < weight[disp]):
#                self.overall_status = disp
#        else:
#            raise ValueError('status must be OK, ReproNoDD, or Hold')
#
#
#    def disposition_checks(self):
#        obi = self.obi
#        set_disposition('OK')
#        status = dict(slot=dict(),
#                      warn=[])
#        for slot_id in range(0, 8):
#            checks = self.checks[slot_id]
#            status['slot'][slot_id] = ''
#            for c in checks:
#                lowmatch = re.match('.*_low$', c)
#                if not lowmatch and not checks[c]['success']:
#                    status['slot'][slot_id] = 'fail check %s' % c
#                    status['warn'].append('slot %d fail check %s' % (slot_id, c))
#        all_dy_rms = np.array([obi.slot[s]['dy_rms'] for s in obi.guide_slots])
#        all_dz_rms = np.array([obi.slot[s]['dz_rms'] for s in obi.guide_slots])
#        n04 = (np.count_nonzero(all_dy_rms > 0.4) + 
#               np.count_nonzero(all_dz_rms > 0.4))
#        n025 = (np.count_nonzero(all_dy_rms > 0.25) + 
#                np.count_nonzero(all_dz_rms > 0.25))
#        if (n04 > 2):
#            status['warn'].append('2 or more stars with RMS > 0.4')
#        if (n025 > 2)
#            status['warn'].append('2 or more stars with RMS > 0.4')
#        if (n04 == 1 & n025 == 1):
#            set_disposition('ReproNoDD')
#            status['warn'].append('repro without bad slot')
#        
#        if (np.count_nonzero(all_dy_rms > 0.25) + 
#            np.count_nonzero(all_dz_rms > 0.25)) > 2:
#            status['warn'].append('2 or more stars with RMS > 0.25')
#        if (np.count_nonzero(all_dy_rms > 0.4)  + 
#            np.count_nonzero(all_dz_rms > 0.25)) > 2:
#            status['warn'].append('2 or more stars with RMS > 0.4')
#
#            
#        self.status = status

        #for fslot in obi.fid_slots:
        #    checks = self.checks[fslot]
        #    for c in checks:
        #        lowmatch = re.match('.*_low$', c)
        #        if not lowmatch and checks[c]['val'] == False:
        #            status['overall'] = False
        #            status['slot'][fslot] = False

                    
        
    

                
        
        
        
#    def fid_checks(self):
#        obi = self.obi
#        slot_list = np.array([s for s in obi.slot if obi.slot[s]['type'] == 'fid'])
#
#        stest = dict(val=len(slot_list) > 0,
#                     name="has fid lights",
#                     id="test_has_fid_lights",
#                     slots=slot_list.tolist())
#        self.tests.append(stest)
#        
#        dy_med = np.array([np.median(obi.slot[s]['dy']) for s in slot_list])
#        dz_med = np.array([np.median(obi.slot[s]['dz']) for s in slot_list])
#        med_limit = 0.4
#        med_idx = np.unique(np.hstack([np.flatnonzero(dy_med > med_limit),
#                                      np.flatnonzero(dz_med > med_limit)]))
#        med_slots = slot_list[med_idx]
#        stest = dict(val=len(med_slots),
#                     name="one or more fid slots mean residual > %f" % med_limit,
#                     id="test_n_fid_med",
#                     slots=med_slots.tolist())
#        self.tests.append(stest)
#
#        dy_rms = np.array([np.std(obi.slot[s]['dy']) for s in slot_list])
#        dz_rms = np.array([np.std(obi.slot[s]['dz']) for s in slot_list])
#        rms_limit = 0.05
#        rms_idx = np.unique(np.hstack([np.flatnonzero(dy_rms > rms_limit),
#                                      np.flatnonzero(dz_rms > rms_limit)]))
#        rms_slots = slot_list[rms_idx]
#        stest = dict(val=len(rms_slots),
#                     name="one or more fid slots RMS dev > %f" % rms_limit,
#                     id="test_n_fid_rms",
#                     slots=rms_slots.tolist())
#        self.tests.append(stest)
#
#
#
#    def star_checks(self):
#        obi = self.obi
#        slot_list = np.array([s for s in obi.slot if obi.slot[s]['type'] == 'guide'])
#        dy_rms = np.array([np.std(obi.slot[s]['dy']) for s in slot_list])
#        dz_rms = np.array([np.std(obi.slot[s]['dz']) for s in slot_list])
#        high_limit = 0.4
#        high_idx = np.unique(np.hstack([np.flatnonzero(dy_rms > high_limit),
#                                      np.flatnonzero(dz_rms > high_limit)]))
#        high_slots = slot_list[high_idx]
#        low_limit = 0.25
#        low_idx = np.unique(np.hstack([np.flatnonzero(dy_rms > low_limit),
#                                      np.flatnonzero(dz_rms > low_limit)]))
#        low_slots = slot_list[low_idx]
#        stest = dict(val=len(high_slots) > 2,
#                     name="two or more star slots rms > %f" % high_limit,
#                     id="test_n_star_high_rms",
#                     slots=high_slots.tolist())
#        self.tests.append(stest)
#        stest = dict(val=len(low_slots) > 2,
#                     name="two or more star slots rms > %f" % low_limit,
#                     id="test_n_star_low_rms",
#                     slots=low_slots.tolist())
#        self.tests.append(stest)
#        stest = dict(val=((len(high_slots) == 1) & (len(low_slots) == 1)), 
#                     name="one star with rms > %f" % high_limit,
#                     id="test_one_star_high_rms",
#                     slots=high_slots.tolist())
#        self.tests.append(stest)
#
#        max_rad = 1.5
#        dy_med = np.array([np.median(obi.slot[s]['dy']) for s in slot_list])
#        dz_med = np.array([np.median(obi.slot[s]['dz']) for s in slot_list])
#        rad_offset_idx = np.flatnonzero(np.sqrt(dy_med**2 + dz_med**2) > max_rad)
#        rad_slots = slot_list[rad_offset_idx]
#        stest = dict(val=len(rad_slots) == 0,
#                     name="all stars rad offset < %f" % max_rad,
#                     id="test_no_star_max_rad",
#                     slots=rad_slots.tolist())
#        self.tests.append(stest)
#        stest = dict(val=len(rad_slots) == 1,
#                     name="one star radial offset  > %f" % max_rad,
#                     id="test_one_star_max_rad",
#                     slots=rad_slots.tolist())
#        self.tests.append(stest)
#        stest = dict(val=len(rad_slots) > 1,
#                     name="radial offset of star(s) > %f" % max_rad,
#                     id="test_multi_star_max_rad",
#                     slots=rad_slots.tolist())
#        self.tests.append(stest)
#
#        dyz_big_lim = 0.7
#        frac_dyz_lim = 0.05
#        dy_big = np.array([obi.slot[s]['dy']
#                                - np.median(obi.slot[s]['dy']) 
#                                > dyz_big_lim for s in slot_list])
#        dz_big = np.array([obi.slot[s]['dz'] 
#                                - np.median(obi.slot[s]['dz']) 
#                                > dyz_big_lim for s in slot_list])
#        raise ValueError
