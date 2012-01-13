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


class Obi:
    def __init__(self, obsparfile, obsdir):
        self.obspar = get_obspar(obsparfile)
        self.obsdir = obsdir
        self.find_aspect_intervals()
        self.process_aspect_intervals()
        for slot in range(0, 8):
            self.plot_slot(slot)

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

    def plot_slot(self, slot):
        y = None
        z = None
        xy_range = None
        qual = np.concatenate(
            [interv.deltas[slot]['qual'] for interv in self.aspect_intervals])
        dy = np.concatenate(
            [interv.deltas[slot]['dy'] for interv in self.aspect_intervals])
        dz = np.concatenate(
            [interv.deltas[slot]['dy'] for interv in self.aspect_intervals])
        mag = np.concatenate(
            [interv.deltas[slot]['mag'] for interv in self.aspect_intervals])
        time = np.concatenate(
            [interv.deltas[slot]['time'] for interv in self.aspect_intervals])
        ai_starts = [interv.deltas[slot]['time'][0]
                     for interv in self.aspect_intervals]

        time0 = time[0]
        dy0 = np.median(dy)
        dz0 = np.median(dz)

        fig = plt.figure(num=slot, figsize=(12, 10))
        #fid_plot = np.abs(np.max(y) - np.min(y)) > 1e-6
        #fid_plot = slot < 3
        ok = qual == 0
        bad = qual != 0

        #if fid_plot:
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

        ayz = fig.add_axes([.05, .7, .20, .20])
        ayz.plot(dy[ok], dz[ok], 'g.')
        ayz.plot(dy[bad], dz[bad], 'r.')
        ayz.set_aspect('equal', 'datalim')
        plt.setp(ayz.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(ayz.get_xticklabels(), fontsize=labelfontsize)
        ayz.set_xlabel('Y offset (arcsec)')
        ayz.set_ylabel('Z offset (arcsec)')

        ayzf = fig.add_axes([.05, .25, .20, .20])
        ayzf.plot(dy[ok], dz[ok], 'g.')
        ayzf.plot(dy[bad], dz[bad], 'r.')
        ayzf.set_aspect('equal', 'datalim')
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


class AspectInterval:
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
        if np.count_nonzero(good) == 0:
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

        r2a = 206264.81
        d2r = 0.017453293
        M0 = 10.32
        C0 = 5263.0

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
                                            qual=qual
                                            )

    def calc_guide_deltas(self):
        from Quaternion import Quat

        M0 = 10.32
        C0 = 5263.0
        r2a = 180. * 3600. / 3.14159265
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

            self.deltas[star['slot']]= {'dy': dy,
                                        'dz': dz,
                                        'time': ceni['time'],
                                        'mag': mag,
                                        'qual': qual}


        
