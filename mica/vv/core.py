from __future__ import division
import os
import re
import pyfits
import pickle
import json
import shelve
import csv
import gzip
import logging
from glob import glob
import numpy as np
import numpy.ma as ma
import tempfile
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import medfilt as medfilt
from scipy.stats import scoreatpercentile

import Ska.Numpy
from Chandra.Time import DateTime
from Ska.Table import read_table
from Ska.astro import sph_dist
from Ska.engarchive import fetch

class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, 'tolist'):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


class InconsistentAspectIntervals(ValueError):
    pass

# integer code version for lightweight database tracking
VV_VERSION = 3

logger = logging.getLogger('vv')


def parse_obspar(file):
    """
    Parse obspar.
    """
# borrowed from telem_archive
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


R2A = 206264.81
D2R = 0.017453293
M0 = 10.32
C0 = 5263.0

fidDrEncFrac1 = 0.68
fidDrEncFrac2 = 0.95 
starDrEncFrac1 = 0.68
starDrEncFrac2 = 0.95

enc_rad_1 = 0.68
enc_rad_2 = 0.99
dyz_big_lim = 0.7
frac_dyz_lim = 0.05

DATA_COLS = ['qual', 'dy', 'dz', 'dr', 'mag', 'time', 'yag', 'zag',
             'ang_y_sm', 'ang_y', 'ang_z_sm', 'ang_z']

save_asol_header = ['ASCDSVER', 'CALDBVER', 'DATE', 'OBI_NUM', 'OBS_ID',
                    'REVISION', 'TSTART', 'TSTOP']

#def rms(data, median):
#    return np.sqrt(np.mean((data - median) ** 2))


def frac_bad(data, median, limit):
    """
    Count data points more than limit from median
    """
    if len(data) == 0:
        return None
    big_offset = abs(data - median) > limit
    frac_bad = len(np.flatnonzero(big_offset)) / len(data)
    return frac_bad


class Obi(object):
    def __init__(self, obsparfile, obsdir, temproot=None):
        self._info = None
        self._errors = []
        self._isdefault = None
        self.obspar = get_obspar(obsparfile)
        self.obsdir = obsdir
        self._process()

    def _process(self):
        self._find_aspect_intervals()
        self._process_aspect_intervals()
        self._concat_slot_data()
        try:
            self._check_over_intervals()
        except InconsistentAspectIntervals as e:
            self._errors.append(e)
        self._sim_data()
        self.slot_report = dict()
        self._label_slots()
        self._agg_slot_data()
        self._get_info()

    def save_plots_and_resid(self):
        self._save_info_json()
        self._save_info_pkl()
        #for slot in self.all_slot_data:
        #    self.plot_slot(slot, save=True, close=True)
        for slot in self.all_slot_data:
            self.plot_slot(slot, save=True, close=True, singles=True)

    def set_dbh(self, dbhandle, slot_table='vv_slots'):
        self.db = dbhandle
        self.slot_table = slot_table

    def set_tbl(self, table):
        self.table = table

    def _save_slot_pkl(self, file=None):
        """
        save the slot residuals as a pkl
        """
        if file is None:
            file = os.path.join(self.tempdir, 'vv_slot.pkl')
        pickle.dump(self.all_slot_data, open(file, 'w'))

    def _set_default(self, isdefault):
        self._isdefault = isdefault

    def _get_default(self):
        return self._isdefault

    isdefault = property(_get_default, _set_default)

    def _just_slot_data(self):
        """
        get just the aggregate values for a slot
        """
        slist = []
        for slot_id in range(0, 8):
            slot = self.slot[slot_id]
            if not 'n_pts' in slot:
                continue
            # ignore the arrays
            save = dict((k, v) for k, v in slot.iteritems()
                        if type(v) not in [np.ndarray, np.ma.core.MaskedArray])
            slist.append(save)
        return slist

    def info(self):
        if self._info is None:
            self._get_info()
        return self._info

    def _aiid_info(self, save_cols=save_asol_header):
        """
        grab some values from the asol header for the top level
        """
        slist = []
        for ai in self.aspect_intervals:
            asol_header = ai.asol_header
            save = dict((k, asol_header[k]) for k in save_cols)
            slist.append(save)
        return slist

# this should probably be handled in mica.archive.asp_l1
    @staticmethod
    def _asp1_lookup(obsid, obi, revision):
        import Ska.DBI
        apstat = Ska.DBI.DBI(dbi='sybase',
                             server='sqlsao',
                             database='axafapstat')
        # take these from the first aspect solution file header
        aspect_1 = apstat.fetchall("""SELECT * FROM aspect_1
                                      WHERE obsid = {obsid}
                                      AND obi = {obi}
                                      AND revision = {revision}
                                   """.format(obsid=obsid,
                                              obi=obi,
                                              revision=revision))
        apstat.conn.close()
        if len(aspect_1) > 1:
            raise ValueError(
                "More than one entry found for obsid/obi/rev in aspect_1")
        if len(aspect_1) == 0:
            logger.warn("obsid / revision not in axafapstat.aspect_1")
            return (None, None)
        return aspect_1[0]['aspect_1_id'], aspect_1[0]['ap_date']

    def _get_info(self):
        """
        get labels for top level
        """
        ai_list = self._aiid_info()
        obsid = ai_list[0]['OBS_ID']
        revision = ai_list[0]['REVISION']
        obi = ai_list[0]['OBI_NUM']
        for ai in ai_list:
            if ai['OBS_ID'] != obsid:
                raise ValueError
            if ai['REVISION'] != revision:
                raise ValueError
        fidprops = [[dict(zip(r.dtype.names,r)) for r in getattr(ai, 'fidprop')]
                    for ai in self.aspect_intervals]
        gsprops = [[dict(zip(r.dtype.names,r)) for r in getattr(ai, 'gsprop')]
                   for ai in self.aspect_intervals]
        self._info = {'obsid': int(obsid),
                     'revision': revision,
                     'tstart': self.obspar['tstart'],
                     'tstop': self.obspar['tstop'],
                      'fidprops': fidprops,
                      'gsprops': gsprops,
                     'sim_z': self.obspar['sim_z'],
                     'sim_z_offset': self.obspar['sim_z_offset'],
                     'ra_pnt': self.obspar['ra_pnt'],
                     'dec_pnt': self.obspar['dec_pnt'],
                     'roll_pnt': self.obspar['roll_pnt'],
                     'instrument': self.obspar['detnam'],
                     'intervals': ai_list,
                     'slots': self.slot_report,
                     'sim': self.sim_report,
                      'vv_version': VV_VERSION,
                      'errors': [str(e) for e in self._errors]}

        try:
            aspect_1_id, ap_date = self._asp1_lookup(obsid, obi, revision)
            self._info['aspect_1_id'] =  aspect_1_id
            self._info['ap_date'] =  str(ap_date)
        except:
            logger.warn("Could not determine aspect_1_id/date from database")

        # we don't care about the DateTimeType for ap_date,
        # so just cast to a string


    def _save_info_json(self, file=None):
        if file is None:
            file = os.path.join(self.tempdir, 'vv_report.json')
        save = self.info()
        jfile = open(file, 'w')
        jfile.write(json.dumps(save, sort_keys=True, indent=4,
                               cls=NumpyAwareJSONEncoder))
        jfile.close()
        logger.info("Saved JSON to {}".format(file))

    def _save_info_pkl(self, file=None):
        if file is None:
            file = os.path.join(self.tempdir, 'vv_report.pkl')
        pfile = open(file, 'w')
        pickle.dump(self.info(), pfile)
        pfile.close()

    def shelve_info(self, file):
        if self.info()['aspect_1_id'] is None:
            logger.warn("Shelving not implemented for obsids without aspect_1_ids")
            return
        s = shelve.open(file)
        s["%s_%s" % (self.info()['obsid'], self.info()['revision'])] \
            = self.info()
        s.close()
        logger.info("Saved to shelve file {}".format(file))

    def slots_to_db(self):
        if self.info()['aspect_1_id'] is None:
            logger.warn("Database save not implemented for obsids without aspect_1_ids")
            return
        save = self.info()
        in_db = self.db.fetchall("""select * from %s where
                                    obsid = %d and revision = %d"""
                                 % (self.slot_table,
                                    save['obsid'], save['revision']))
        if len(in_db):
            self.db.execute("""delete from %s where
                               obsid = %d and revision = %d"""
                            % (self.slot_table,
                               save['obsid'], save['revision']))
        for slot in save['slots']:
            slot.update(obsid=save['obsid'])
            slot.update(revision=save['revision'])
            self.db.insert(slot, self.slot_table)

    @staticmethod
    def _get_ccd_temp(tstart, tstop):
        temps = fetch.MSID('AACCCDPT', tstart, tstop)
        if not len(temps.vals):
            return None
        else:
            return np.mean(temps.vals) - 273.15

    def slots_to_table(self):
        save = self.info()
        if save['aspect_1_id'] is None:
            logger.warn("Table save not implemented for obsids without aspect_1_ids")
            return
        mean_aacccdpt = self._get_ccd_temp(save['tstart'], save['tstop'])
        if mean_aacccdpt is None:
            logger.warn(
                "No AACCCDPT data for {}.  Skipping HDF5 ingest".format(
                    save['obsid']))
            return
        # add obsid and revision to the slot dict
        for slot_str in save['slots']:
            slot = save['slots'][slot_str]
            slot.update(dict((k, save[k])
                             for k in self.table.dtype.names
                             if k in save))
            slot.update(dict(mean_aacccdpt=mean_aacccdpt))
            slot.update(isdefault=self.isdefault)
            if (save['slots'][slot_str]['type'] == 'FID'):
                save['slots'][slot_str]['used'] = 1
            else:
                save['slots'][slot_str]['used'] = save['slots'][slot_str]['cel_loc_flag']
        # make a recarray
        save_rec = np.rec.fromrecords(
            [[save['slots'][slot_str].get(k) for k in self.table.dtype.names]
             for slot_str in save['slots'] if 'n_pts' in save['slots'][slot_str]],
            dtype=self.table.dtype)

        have_obsid_coord = self.table.getWhereList('(obsid == %d)'
                                                   % (save['obsid']),
                                                   sort=True)
        # if there are previous records for this obsid
        if len(have_obsid_coord):
            logger.info("obsid %d is in table" % save['obsid'])
            obsid_rec = self.table.readCoordinates(have_obsid_coord)
            # if the current entry is default, mark other entries as
            # not-default
            if self.isdefault:
                obsid_rec['isdefault'] = 0
                self.table.modifyCoordinates(have_obsid_coord, obsid_rec)
            # if we already have the revision, update in place
            if np.any(obsid_rec['revision'] == save['revision']):
                rev_coord = self.table.getWhereList(
                    '(obsid == %d) & (revision == %d)'
                    % (save['obsid'], save['revision']),
                    sort=True)
                if len(rev_coord) != len(save_rec):
                    raise ValueError(
                        "Could not update; different number of slots")
                logger.info("updating obsid %d rev %d in place"
                       % (save['obsid'], save['revision']))
                self.table.modifyCoordinates(rev_coord, save_rec)
            else:
                self.table.append(save_rec)
        else:
            self.table.append(save_rec)
        self.table.flush()

    def _aiid_from_asol(self, asol_file, obsdir):
        hdulist = pyfits.open(asol_file)
        header = hdulist[1].header
            # skip files that aren't in the obspar range
        if ((self.obspar['tstart'] >= header['TSTOP'])
            or (self.obspar['tstop'] <= header['TSTART'])):
            return None
        aiid_match = re.search('(pcadf\d+[^_]*)_', asol_file)
        if aiid_match:
            return dict(id=aiid_match.group(1),
                        dir=obsdir)

    def _find_aspect_intervals(self):
        obsdir = self.obsdir
        asol_files = sorted(glob(os.path.join(obsdir, 'pcad*asol*')))
        self.aiids = []
        if len(asol_files):
            for file in asol_files:
                self.aiids.append(self._aiid_from_asol(file, obsdir))
        ASP_dirs = sorted(glob(os.path.join(obsdir, 'ASP_L1_*')))
        if len(ASP_dirs):
            for dir in ASP_dirs:
                max_out = max(sorted(glob(os.path.join(dir, 'out*'))))
                asol_files = sorted(glob(os.path.join(max_out, 'pcad*asol*')))
                if len(asol_files):
                    for file in asol_files:
                        self.aiids.append(self._aiid_from_asol(file, max_out))


    def _sim_data(self):
        ai_0 = self.aspect_intervals[0]
        sim_keys = ai_0.sim.keys()
        all_sim = dict()
        for d in sim_keys:
            if len(ai_0.sim[d]):
                if isinstance(ai_0.sim[d], np.ma.MaskedArray):
                    all_sim[d] = ma.concatenate([
                        ai.sim[d] for ai in self.aspect_intervals])
                else:
                    all_sim[d] = np.concatenate([
                        ai.sim[d] for ai in self.aspect_intervals])
        self.sim_data = all_sim
        self.sim_report = dict(max_medf_dy=np.max(all_sim['medf_dy']),
                               min_medf_dy=np.min(all_sim['medf_dy']),
                               max_medf_dz=np.max(all_sim['medf_dz']),
                               min_medf_dz=np.min(all_sim['medf_dz']),
                               max_d_dy=np.max(all_sim['d_dy']),
                               max_d_dz=np.max(all_sim['d_dz']))




    def _process_aspect_intervals(self):
        self.aspect_intervals = []
        for aiid in self.aiids:
            self.aspect_intervals.append(
                AspectInterval(aiid['id'], aiid['dir']))

    def _concat_slot_data(self):
        slot_data = dict()
        for slot in range(0, 8):
            cslot = dict()
            for d in DATA_COLS:
                slotval = [i.deltas[slot][d] for i in self.aspect_intervals
                           if slot in i.deltas]
                if len(slotval):
                    if isinstance(slotval[0], np.ma.MaskedArray):
                        cslot[d] = ma.concatenate(slotval)
                    else:
                        cslot[d] = np.concatenate(slotval)
            # only make a top-level slot key if there is some kind of data
            if len(cslot.keys()) and len(cslot['time']):
                slot_data[slot] = cslot
        self.all_slot_data = slot_data


    def _agg_slot_data(self):
        all_slot = self.all_slot_data
        for slot_id in range(0, 8):
            slot_report = self.slot_report[str(slot_id)]
            # get status info from guide or fid props
            # just use the first aspect interval and assume the rest are the
            # same (checked with _check_over_intervals)
            if slot_report['type'] == 'FID':
                slot_fidprop = self.aspect_intervals[0].fidprop[
                    self.aspect_intervals[0].fidprop['slot'] == slot_id][0]
                slot_report['id_status'] = slot_fidprop['id_status'].strip()
            else:
                slot_gsprop = self.aspect_intervals[0].gsprop[
                    self.aspect_intervals[0].gsprop['slot'] == slot_id][0]
                slot_report['id_status'] = slot_gsprop['id_status'].strip()
                slot_report['cel_loc_flag'] = slot_gsprop['cel_loc_flag']

            if slot_id not in all_slot:
                continue
            slot_data = all_slot[slot_id]
            if not ('dy' in slot_data and 'dz' in slot_data and 'mag' in slot_data):
                continue
            # these should only be calculated over good data, right?
            qual = dict(dy=slot_data['qual'],
                        dz=slot_data['qual'],
                        dr=slot_data['qual'])
            for axdir in ['dr', 'dy', 'dz', 'mag']:
                if axdir == 'mag':
                    data = slot_data['mag']
                else:
                    data = ma.array(slot_data[axdir])
                    data[qual[axdir] != 0] = ma.masked
                slot_report['%s_n_samples' % axdir] = len(data)
                slot_report['%s_bad_samples' % axdir] = len(np.flatnonzero(data.mask))
                smean = ma.mean(data)
                slot_report['%s_mean' % axdir] = smean
                med = ma.median(data)
                slot_report['%s_med' % axdir] = med.item()
                srms = ma.sqrt(ma.mean((data - med) ** 2))
                slot_report['%s_rms' % axdir] = srms
                bad_frac = frac_bad(data, med, dyz_big_lim)
                slot_report['frac_%s_big' % axdir] = bad_frac
            slot_report['frac_%s_big' % axdir] = bad_frac
            slot_report['rad_off'] = np.sqrt(slot_report['dy_med'] ** 2
                                      + slot_report['dz_med'] ** 2)
            slot_report['n_pts'] = len(slot_data['dy'])
            if slot_report['type'] == 'FID':
                slot_report['mean_y'] = np.mean(slot_data['ang_y_sm'])
                slot_report['mean_z'] = np.mean(slot_data['ang_z_sm'])
                slot_report['enc_rad1'] = scoreatpercentile(slot_data['dr'],
                                                            fidDrEncFrac1 * 100)
                slot_report['enc_rad2'] = scoreatpercentile(slot_data['dr'],
                                                            fidDrEncFrac2 * 100)
            else:
                slot_report['mean_y'] = np.mean(slot_data['ang_y'])
                slot_report['mean_z'] = np.mean(slot_data['ang_z'])
                slot_report['enc_rad1'] = scoreatpercentile(slot_data['dr'],
                                                            starDrEncFrac1 * 100)
                slot_report['enc_rad2'] = scoreatpercentile(slot_data['dr'],
                                                            starDrEncFrac2 * 100)




    def _label_slots(self):
        ai = self.aspect_intervals[0]
        self.guide_list = list(getattr(ai, 'gsprop').slot)
        for gs in getattr(ai, 'gsprop'):
            self.slot_report[str(gs.slot)] = dict(
                slot=gs.slot,
                type='GUIDE')
        if getattr(ai, 'fidprop') is not None:
            self.fid_list = list(getattr(ai, 'fidprop').slot)
            for fl in getattr(ai, 'fidprop'):
                self.slot_report[str(fl.slot)] = dict(
                    slot=fl.slot,
                    type='FID')
        else:
            self.fid_list = []

    def _check_over_intervals(self):
        """
        Check all aspect intervals and confirm that slots don't
        drop out, change status, or change type.
        """
        for t in ('gsprop', 'fidprop'):
            if getattr(self.aspect_intervals[0], t) is None:
                continue
            slot_id = getattr(self.aspect_intervals[0], t).slot
            slot_status = getattr(self.aspect_intervals[0], t).id_status
            slot_type = 'FID'
            if t == 'gsprop':
                slot_type = getattr(self.aspect_intervals[0], t).type
                slot_cel_loc_flag = getattr(self.aspect_intervals[0], t)['cel_loc_flag']
            for ai in self.aspect_intervals[1:]:
                if len(slot_id) != len(getattr(ai, t).slot):
                    raise InconsistentAspectIntervals(
                        "differing %s slots across aspect intervals" % t)
                if ((len(slot_id) == len(getattr(ai, t).slot)) &
                    (not np.all([slot_id == getattr(ai, t).slot]))):
                    raise InconsistentAspectIntervals(
                        "differing %s slots across aspect intervals" % t)
                if ((len(slot_id) == len(getattr(ai, t).slot)) &
                    (not np.all([slot_status == getattr(ai, t).id_status]))):
                    raise InconsistentAspectIntervals(
                        "differing %s status across aspect intervals" % t)
                if t == 'gsprop':
                    if ((len(slot_id) == len(getattr(ai, t).slot)) &
                        (not np.all([slot_type == getattr(ai, t).type]))):
                        raise InconsistentAspectIntervals(
                            "differing %s type across aspect intervals" % t)
                    if ((len(slot_id) == len(getattr(ai, t).slot)) &
                        (not np.all([slot_cel_loc_flag == getattr(ai, t)['cel_loc_flag']]))):
                        raise InconsistentAspectIntervals(
                            "differing %s cel_loc_flag across aspect intervals" % t)


    def plot_slot(self, slot_num, plotdir=None, save=False, close=False, singles=False):
        if plotdir is None and save:
            plotdir = self.tempdir
        y = None
        z = None
        xy_range = None
        fid_plot = (slot_num in self.fid_list)
        if slot_num not in self.all_slot_data:
            logger.info("Nothing to plot for slot %d" % slot_num)
            return None, None
        (qual, dy, dz, mag, time) = [
            self.all_slot_data[slot_num][x] for x in
            ['qual', 'dy', 'dz', 'mag', 'time']]
        if not fid_plot:
            (yag, zag, ang_y_sm, ang_z_sm) = [
                self.all_slot_data[slot_num][x] for x in
                ['yag', 'zag', 'ang_y_sm', 'ang_z_sm']]
        ai_starts = [interv.deltas[slot_num]['time'][0]
                     for interv in self.aspect_intervals
                     if len(interv.deltas[slot_num]['time'])]
        time0 = time[0]
        dy0 = np.median(dy)
        dz0 = np.median(dz)

        #fid_plot = np.abs(np.max(y) - np.min(y)) > 1e-6
        ok = qual == 0
        bad = qual != 0

        if fid_plot:
            y = self.all_slot_data[slot_num]['ang_y_sm']
            z = self.all_slot_data[slot_num]['ang_z_sm']

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
        timepad = 0.05 * (plottime[-1] - plottime[0])
        labelfontsize = 10
        axes = dict()

        max_y = np.max(dy[ok])
        min_y = np.min(dy[ok])
        max_z = np.max(dz[ok])
        min_z = np.min(dz[ok])
        ymid = (max_y + min_y) / 2
        zmid = (max_z + min_z) / 2
        max_extent = np.max([max_y - min_y, max_z - min_z]) / 2
        extent = max_extent * 1.10
        if extent < (xy_range / 2):
            extent = xy_range / 2
        plot_dict = dict()

        if singles:
            fig1 = plt.figure(figsize=(3,3))
            ayz = plt.subplot(1, 1, 1)
            plot_dict['ayz'] = fig1
        else:
            fig = plt.figure(figsize=(14, 10))
            ayz = fig.add_axes([.05, .7, .20, .20], aspect='equal')
            plot_dict['all'] = fig
        axes['yz'] = ayz
        ayz.plot(dy[ok], dz[ok], 'g.')
        ayz.plot(dy[bad], dz[bad], 'r.')
        ayz.grid()
        plt.setp(ayz.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(ayz.get_xticklabels(), fontsize=labelfontsize)
        plt.setp(ayz.get_xticklabels(), rotation=30)
        plt.setp(ayz.get_xticklabels(), horizontalalignment='right')
        ayz.set_xlabel('Y offset (arcsec)')
        ayz.set_ylabel('Z offset (arcsec)')
        # set limits to include all of the "ok" data
        # use the middle of the range of the data, not the median
        ayz.set_xlim(ymid - extent, ymid + extent)
        ayz.set_ylim(zmid - extent, zmid + extent)
        plt.draw()
        if singles:
            plt.tight_layout()
            if save:
                plotfile = os.path.join(self.tempdir,
                                        "slot_{}_yz.png".format(slot_num))
                plt.savefig(plotfile)
                logger.info("Saved plot {}".format(plotfile))
            if close:
                plt.close(fig1)
                

        if singles:
            fig2 = plt.figure(figsize=(3,3))
            ayzf = plt.subplot(1, 1, 1)
            plot_dict['ayzf'] = fig2
        else:
            ayzf = fig.add_axes([.05, .25, .20, .20], aspect='equal')
        axes['yz_fixed'] = ayzf
        ayzf.plot(dy[ok], dz[ok], 'g.')
        ayzf.plot(dy[bad], dz[bad], 'r.')
        ayzf.grid()
        plt.setp(ayzf.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(ayzf.get_xticklabels(), fontsize=labelfontsize)
        plt.setp(ayzf.get_xticklabels(), rotation=30)
        plt.setp(ayzf.get_xticklabels(), horizontalalignment='right')
        circle = plt.Circle((dy0, dz0), radius=circ_rad, fill=False)
        ayzf.add_patch(circle)
        ayzf.set_xlabel('Y offset (arcsec)')
        ayzf.set_ylabel('Z offset (arcsec)')
        # set limits to fixed range
        if xy_range is not None:
            ayzf.set_xlim([dy0 - xy_range, dy0 + xy_range])
            ayzf.set_ylim([dz0 - xy_range, dz0 + xy_range])
            plt.draw()
        if singles:
            plt.tight_layout()
            if save:
                plotfile = os.path.join(self.tempdir,
                                        "slot_{}_yzf.png".format(slot_num))
                plt.savefig(plotfile)
                logger.info("Saved plot {}".format(plotfile))
            if close:
                plt.close(fig2)

        if singles:
            if fid_plot:
                fig3 = plt.figure(figsize=(7.9,2.5))
            else:
                fig3 = plt.figure(figsize=(7,2.5))
            ay = plt.subplot(1, 1, 1)
            plot_dict['ay'] = fig3
        else:
            ay = fig.add_axes([.30, .7, .62, .25])
        axes['dy'] = ay
        ay.plot(plottime[ok], dy[ok], 'g.')
        ay.plot(plottime[bad], dy[bad], 'r.')
        plt.setp(ay.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(ay.get_xticklabels(), visible=False)
        ay.set_ylabel('Y offsets(dy) (arcsec)')
        ay.set_ylim(ymid - extent, ymid + extent)
        ay.set_xlim(plottime[0] - timepad, plottime[-1] + timepad)
        plt.draw()
        if len(ai_starts) > 1:
            for t in ai_starts:
                s_t = (t - time0) / 1000.
                ay.plot([s_t, s_t], ay.get_ylim(), color='blue',
                        linestyle='dashed')
        ay.grid()

        if y is not None:
            ay2 = ay.twinx()
            ay2.plot(plottime, y)
            ay2.autoscale(enable=False, axis='x')
            ay2.set_xlim(plottime[0] - timepad, plottime[-1] + timepad)
            ay2y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            ay2.yaxis.set_major_formatter(ay2y_formatter)
            plt.setp(ay2.get_yticklabels(), fontsize=labelfontsize,
                     color='blue')
            ay2.set_ylabel('centroid y angle', color='blue')
        if singles:
            plt.tight_layout()
            if save:
                plotfile = os.path.join(self.tempdir,
                                        "slot_{}_y.png".format(slot_num))
                plt.savefig(plotfile)
                logger.info("Saved plot {}".format(plotfile))
            if close:
                plt.close(fig3)

        if singles:
            if fid_plot:
                fig4 = plt.figure(figsize=(7.9,2.5))
            else:
                fig4 = plt.figure(figsize=(7,2.5))
            az = plt.subplot(1, 1, 1)
            plot_dict['az'] = fig4
        else:
            az = fig.add_axes([.30, .4, .62, .25], sharex=ay)
        axes['dz'] = az
        az.plot(plottime[ok], dz[ok], 'g.')
        az.plot(plottime[bad], dz[bad], 'r.')
        plt.setp(az.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(az.get_xticklabels(), visible=False)
        az.set_ylabel('Z offsets(dz) (arcsec)')
        az.set_ylim(zmid - extent, zmid + extent)
        az.set_xlim(plottime[0] - timepad, plottime[-1] + timepad)
        plt.draw()
        if len(ai_starts) > 1:
            for t in ai_starts:
                s_t = (t - time0) / 1000.
                az.plot([s_t, s_t], az.get_ylim(), color='blue',
                        linestyle='dashed')
        az.grid()

        if z is not None:
            az2 = az.twinx()
            az2.plot(plottime, z)
            az2.autoscale(enable=False, axis='x')
            az2.set_xlim(plottime[0] - timepad, plottime[-1] + timepad)
            az2y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            az2.yaxis.set_major_formatter(az2y_formatter)
            plt.setp(az2.get_yticklabels(), fontsize=labelfontsize,
                     color='blue')
            az2.set_ylabel('centroid z angle', color='blue')
        if singles:
            plt.tight_layout()
            if save:
                plotfile = os.path.join(self.tempdir,
                                        "slot_{}_z.png".format(slot_num))
                plt.savefig(plotfile)
                logger.info("Saved plot {}".format(plotfile))
            if close:
                plt.close(fig4)

        if singles:
            fig5 = plt.figure(figsize=(7,3))
            am = plt.subplot(1, 1, 1)
            plot_dict['am'] = fig5
        else:
            am = fig.add_axes([.30, .1, .62, .25], sharex=ay)
        axes['mag'] = am
        plt.setp(am.get_yticklabels(), fontsize=labelfontsize)
        plt.setp(am.get_xticklabels(), fontsize=labelfontsize)
        am.plot(plottime[ok], mag[ok], color='green')
        am.plot(plottime[bad], mag[bad], 'r.')
        am.set_ylim(am.get_ylim()[::-1])
        if len(ai_starts) > 1:
            for t in ai_starts:
                s_t = (t - time0) / 1000.
                am.plot([s_t, s_t], am.get_ylim(), color='blue',
                        linestyle='dashed')
        am.grid()
        am.set_ylabel('Magnitude(mag)')
        am.set_xlabel('Time(ksec)')
        am.autoscale(enable=False, axis='x')
        am.set_xlim(plottime[0] - timepad, plottime[-1] + timepad)
        plt.draw()
        if singles:
            plt.tight_layout()
            if save:
                plotfile = os.path.join(self.tempdir,
                                        "slot_{}_m.png".format(slot_num))
                plt.savefig(plotfile)
                logger.info("Saved plot {}".format(plotfile))
            if close:
                plt.close(fig5)

        if not singles:
            plt.suptitle('Obsid %d Slot %d Residuals' % (
                    self.info()['obsid'], slot_num))
            if save:
                plotfile = os.path.join(self.tempdir,
                                        "slot_{}.png".format(slot_num))
                plt.savefig(plotfile)
                logger.info("Saved plot {}".format(plotfile))
            if close:
                plt.close(fig)

        if not fid_plot and not save:
            cenifig = plt.figure(figsize=(12, 10))
            ceniy = cenifig.add_subplot(2, 1, 1, sharex=ay)
            ceniy.plot(plottime[ok], yag[ok], 'b.')
            ceniy.plot(plottime[ok], ang_y_sm[ok] * 3600, 'g.')
            ceniy.set_ylabel('Centroid Y (arcsec)')
            ceniy.set_xlabel('Time(ksec)')
            ceniz = cenifig.add_subplot(2, 1, 2, sharex=ay)
            ceniz.plot(plottime[ok], zag[ok], 'b.')
            ceniz.plot(plottime[ok], ang_z_sm[ok] * 3600, 'g.')
            ceniz.set_ylabel('Centroid Z (arcsec)')
            ceniz.set_xlabel('Time(ksec)')
            ceniy.autoscale(enable=False, axis='x')
            ceniz.autoscale(enable=False, axis='x')
            ceniy.set_xlim(plottime[0] - timepad, plottime[-1] + timepad)
            ceniz.set_xlim(plottime[0] - timepad, plottime[-1] + timepad)

            plt.suptitle('Slot %d Centroids (aspect sol in blue)' % slot_num)

        return plot_dict

AI_DEFAULT_CONF = {'obc': None,
                   'alg': 8,
                   'noacal': None}

class AspectInterval(object):
    def __init__(self, aiid, aspdir, opt=None):
        self.aiid = aiid
        self.aspdir = aspdir
        self.opt = opt if opt else AI_DEFAULT_CONF
        self._read_in_data()
        self._read_in_log()
        self._calc_guide_deltas()
        if self.fidprop is not None:
            self._calc_fid_deltas()
            self._calc_sim_offset()

    def _get_prop(self, propname, propstring):
        "Read gsprops or fidprops file"
        datadir = self.aspdir
        logger.debug('Reading %s stars' % propname)
        gsfile = glob(os.path.join(
                datadir, "%s_%s1.fits*" % (self.aiid, propstring)))[0]
        # don't filter for only good stars at this point
        prop = read_table(gsfile)
        info = []
        hdulist = pyfits.open(os.path.join(datadir, gsfile))
        header = hdulist[1].header
        proptable = read_table(os.path.join(datadir, gsfile))
        for gs in proptable:
            saveprop = dict(slot=gs['slot'],
                            id_status=gs['id_status'],
                            tstart=header['TSTART'],
                            tstop=header['TSTOP'])
            if 'type' in gs.dtype.names:
                saveprop['type'] = gs['type']
            info.append(saveprop)
        return (prop, info, header)

    def _read_ocat_stars(self):
        import Ska.DBI
        obsid = int(self.asol_header['OBS_ID'])
        obi = int(self.asol_header['OBI_NUM'])
        ocat_db = Ska.DBI.DBI(dbi='sybase', server='sqlsao', database='axafocat')
        stars = ocat_db.fetchall("select * from stars where "
                                 "obsid = {} and obi = {} "
                                 "and type != 0".format(obsid, obi))
        ocat_db.conn.close()
        if len(np.unique(stars['obi'])) > 1:
            raise ValueError("Multi-obi observation.  OCAT stars unhelpful to identify missing slot")
        return stars

    def _identify_missing_slot(self, slot):
        datadir = self.aspdir
        adat_files = glob(os.path.join(datadir, "pcadf*N???_adat{}1.fits*".format(slot)))
        if not len(adat_files):
            return None
        hdulist = pyfits.open(adat_files[0])
        header = hdulist[1].header
        if header['IMGTYPE'] == 0:
            return 'GUIDE'
        if header['IMGTYPE'] == 1:
            return 'FID'
        if header['IMGTYPE'] == 2:
            return 'MONITOR'
        raise ValueError("Slot {} could not be identified from image data".format(slot))


    def _read_in_data(self):
        aiid = self.aiid
        datadir = self.aspdir
        opt = self.opt

        logger.info('Processing aspect interval {}'.format(
                self.aiid))
        logger.debug('Reading aspect solution and header')
        #if opt['obc']:
        #    asol = read_table(glob(
        #            os.path.join(datadir, "%s_osol1.fits*" % aiid))[0])
        #else:
        asol_file = glob(
            os.path.join(datadir, "%s_asol1.fits*" % aiid))[0]
        asol = read_table(asol_file)
        # Add code to handle first processing of 16091 with
        # non-confirming asol file
        if ('dtheta' not in asol.dtype.names
            and 'DTHETA' in asol.dtype.names):
            asol = Ska.Numpy.add_column(asol, 'dtheta', asol['DTHETA'])
        hdulist = pyfits.open(asol_file)
        header = hdulist[1].header
        self.asol_header = header
        self.asol = asol

        logger.debug('Reading aspect quality')
        self.aqual = read_table(glob(
                os.path.join(datadir, "%s_aqual1.fits*" % aiid))[0])

        #if opt['noacal']:
        #    aca_misalign = np.array([[1.0,0,0], [0,1,0],[0,0,1]])
        #else:
        logger.debug('Reading ACA and FTS align file')
        acal = read_table(glob(
                os.path.join(datadir, "%s_acal1.fits*" % aiid))[0])
        self.aca_misalign = acal['aca_misalign'].reshape(3, 3)
        self.fts_misalign = acal['fts_misalign'].reshape(3, 3)
        self.acal = acal

        logger.debug('Reading Centroids')
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

        logger.debug('Reading gyro data')
        self.gdat = read_table(glob(
                os.path.join(datadir, "%s_gdat1.fits*" % aiid))[0])


        (self.gsprop, self.gspr_info, self.h_gspr) \
            = self._get_prop('guide', 'gspr')
        (self.fidprop, self.fidpr_info, self.h_fidpr) \
            = self._get_prop('fid', 'fidpr')

        missing_slots = [slot for slot in np.unique(self.cen['slot'])
                         if (slot not in self.gsprop['slot']
                             and slot not in self.fidprop['slot'])]

        # Nothing else to do if there are no missing slots
        if not len(missing_slots):
            return

        try:
            ocat_stars = self._read_ocat_stars()
        except:
            logger.warn('Could not get OCAT stars from database')
            logger.warn('Skipping checks for missing slots')
            return
        tstart = self.asol_header['TSTART']
        agasc_equinox = DateTime('2000:001:00:00:00.000')
        dyear = (DateTime(tstart) - agasc_equinox) / 365.25
        pm_to_degrees = dyear / (3600. * 1000.)
        missing_info = []
        for slot in missing_slots:
            #stype = self._identify_missing_slot(slot)
            #if stype is None:
            #    logger.warn("No image data to identify missing slot {}".format(slot))
            #    logger.warn("Skipping slot")
            #missing_types.append(stype)
            ocat_info = ocat_stars[ocat_stars['slot'] == slot][0]
            if ocat_info['type'] == 3:
                logger.info("Missing slot is MONITOR.  Skipping...")
                continue
            if ocat_info['type'] == 2:
                logger.warn("Missing slot is FID.")
                mock_prop = dict(
                                 id_status='OMITTED',
                                 slot=slot,
                                 id_string=ocat_info['id'],
                                 id_num=ocat_info['id'],
                                 ang_y_nom=ocat_info['y_ang'],
                                 ang_z_nom=ocat_info['z_ang'],
                                 mag_i_cmd=0,
                                 mag_i_avg=0,
                                 mag_i_min=0,
                                 mag_i_max=0,
                                 p_lsi=np.array([0,0,0]),
                                 )
                self.fidprop.resize(len(self.fidprop) + 1)
                self.fidprop[-1] = np.rec.fromrecords([[mock_prop[col]
                                                       for col in self.fidprop.dtype.names]],
                                                     dtype=self.fidprop.dtype)
                self.fidpr_info.append(dict(slot=slot,
                                            tstart=self.asol_header['TSTART'],
                                            tstop=self.asol_header['TSTOP'],
                                            id_status='OMITTED'))

            if ocat_info['type'] == 1:
                import agasc
                star_info = agasc.get_star(ocat_info['id'])
                mock_prop = dict(cel_loc_flag=0,
                                 id_status='OMITTED',
                                 agasc_id=ocat_info['id'],
                                 slot=slot,
                                 type='GUIDE     ',
                                 spectral_type='NONE      ',
                                 ra_offset_obs=0,
                                 dec_offset_obs=0,
                                 mag_aca_avg=0,
                                 mag_aca_min=0,
                                 mag_aca_max=0,
                                 aspq1_obs=0,
                                 spoil_radius=0,
                                 spoil_angle=0,
                                 spoil_mag_aca=0)
                for col in ['ra', 'dec', 'pos_err', 'pm_ra', 'pm_dec',
                            'plx', 'plx_err', 'mag_aca', 'mag_aca_err',
                            'class', 'mag', 'mag_err', 'mag_band', 'color1',
                            'color1_err', 'var', 'aspq1', 'aspq2', 'aspq3']:
                    mock_prop.update({col: star_info[col.upper()]})
                ra = star_info['RA']
                dec = star_info['DEC']
                if star_info['PM_RA'] != -9999:
                    ra = ra + star_info['PM_RA'] * pm_to_degrees
                if star_info['PM_DEC'] != -9999:
                    dec = dec + star_info['PM_DEC'] * pm_to_degrees
                mock_prop.update({'pos_eci': [np.cos(np.radians(ra)) * np.cos(np.radians(dec)),
                                              np.sin(np.radians(ra)) * np.cos(np.radians(dec)),
                                              np.sin(np.radians(dec))],
                                  'ra_corr': ra,
                                  'dec_corr': dec})

                self.gsprop.resize(len(self.gsprop) + 1)
                self.gsprop[-1] = np.rec.fromrecords([[mock_prop[col]
                                                       for col in self.gsprop.dtype.names]],
                                                     dtype=self.gsprop.dtype)
                self.gspr_info.append(dict(slot=slot,
                                           tstart=self.asol_header['TSTART'],
                                           tstop=self.asol_header['TSTOP'],
                                           type='GUIDE     ',
                                           id_status='OMITTED'))



    def _read_in_log(self):
        aiid = self.aiid
        datadir = self.aspdir
        try:
            id_end = re.match('\w+(\D\d+N\d{3})', aiid)
            logfiles = glob(os.path.join(datadir, "*%s.log*" % id_end.group(1)))
            logfile = logfiles[0]
            try:
                lines = gzip.open(logfile).readlines()
            except IOError:
                lines = open(logfile).readlines()
            self.log = lines
        except:
            logger.info("Did not find/read log file")


    def _calc_fid_deltas(self):
        asol = self.asol
        h_fidpr = self.h_fidpr
        aca_misalign = self.aca_misalign
        fts_misalign = self.fts_misalign
        fidprop = self.fidprop
        cen = self.cen
        fidpr_info = self.fidpr_info
        integ_time = self.integ_time

        logger.debug('Calculating fid solution quality')

        lsi0_stt = [h_fidpr['LSI0STT%d' % x] for x in [1, 2, 3]]
        stt0_stf = [h_fidpr['STT0STF%d' % x] for x in [1, 2, 3]]
        rrc0_fc_x = h_fidpr['RRC0FCX']

        M = np.dot(aca_misalign, fts_misalign)

        rot_x = np.zeros([3, 3])
        rot_x[0, 0] = 1
        for fid in fidprop:
            logger.info("Processing fid %s in slot %d " % (
                fid['id_string'], fid['slot']))
            p_lsi = fid['p_lsi']
            p_stf = p_lsi + lsi0_stt + stt0_stf

            ok = cen['slot'] == fid['slot']
            ceni = cen[ok]
            asol_cen_dy = np.interp(ceni['time'], asol['time'], asol['dy'])
            asol_cen_dz = np.interp(ceni['time'], asol['time'], asol['dz'])
            asol_cen_dtheta = (np.interp(ceni['time'],
                                         asol['time'], asol['dtheta'])
                               * D2R)
 
            rot_x = np.zeros([len(ceni['time']), 3, 3])
            s_th = np.sin(asol_cen_dtheta)
            c_th = np.cos(asol_cen_dtheta)
            rot_x[:, 0, 0] = 1.0
            rot_x[:, 1, 1] = c_th
            rot_x[:, 2, 1] = s_th
            rot_x[:, 1, 2] = -s_th
            rot_x[:, 2, 2] = c_th
            p_fc = np.dot(rot_x.transpose(0, 2, 1), p_stf)
            p_fc[:, 1] = p_fc[:, 1] + asol_cen_dy
            p_fc[:, 2] = p_fc[:, 2] + asol_cen_dz
            d_fc = p_fc
            d_fc[:, 0] = d_fc[:, 0] - rrc0_fc_x
            d_fc = -d_fc
            d_aca = np.dot(d_fc, M.transpose())
            yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
            zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
            dy = ceni['ang_y_sm'] * 3600 - yag
            dz = ceni['ang_z_sm'] * 3600 - zag
            dr = sph_dist(yag / 3600,
                          zag / 3600,
                          ceni['ang_y_sm'],
                          ceni['ang_z_sm']) * 3600
            slot_fidpr = [pr for pr in fidpr_info if pr['slot'] == fid['slot']]
            if not slot_fidpr:
                raise ValueError("No FIDPR info found for slot %d"
                                 % fid['slot'])
            mag = ma.zeros(len(ceni['counts']))
            mag[:] = ma.masked
            good_mag = medfilt(M0 - 2.5
                               * np.log10(ceni['counts'][ceni['counts'] > 10.0]
                                          / integ_time
                                          / C0), 3)
            mag[ceni['counts'] > 10] = good_mag
            self.deltas[fid['slot']] = dict(time=ceni['time'],
                                            dy=dy,
                                            dz=dz,
                                            dr=dr,
                                            yag=yag,
                                            zag=zag,
                                            mag=mag,
                                            qual=ceni['status'],
                                            ang_y_sm=ceni['ang_y_sm'],
                                            ang_z_sm=ceni['ang_z_sm'],
                                            ang_y=ceni['ang_y'],
                                            ang_z=ceni['ang_z'],
                                            )

    def _calc_guide_deltas(self):
        from mica.quaternion import Quat
        self.deltas = {}
        asol = self.asol
        cen = self.cen
        gsprop = self.gsprop
        aca_misalign = self.aca_misalign
        integ_time = self.integ_time

        logger.debug('Interpolating quaternions')
        q_att = np.array([np.interp(cen['time'],
                                    asol['time'],
                                    asol['q_att'][:, ax])
                          for ax in range(0, 4)]).transpose()
        for star in gsprop:
            logger.info('Processing {} star in slot {}'.format(
                    star['id_status'], star['slot']))
            ok = cen['slot'] == star['slot']
            ceni = cen[ok]
            logger.debug('Found %d centroids ' % len(ceni))
            if len(ceni) < 2:
                continue
            q_atts = Quat(q_att[ok])
            Ts = q_atts.transform
            # use ang_y or ang_y_sm?
            #inside = np.dot(aca_misalign, Ts.transpose(0,2,1)).transpose(1,0,2)
            d_aca = np.dot(np.dot(aca_misalign, Ts.transpose(0, 2, 1)),
                           star['pos_eci']).transpose()
            yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
            zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
            dy = ceni['ang_y'] * 3600 - yag
            dz = ceni['ang_z'] * 3600 - zag
            dr = sph_dist(yag / 3600,
                          zag / 3600,
                          ceni['ang_y'],
                          ceni['ang_z']) * 3600
            mag = ma.zeros(len(ceni['counts']))
            mag[:] = ma.masked
            good_mag = medfilt(M0 - 2.5
                               * np.log10(ceni['counts'][ceni['counts'] > 10.0]
                                          / integ_time
                                          / C0), 3)
            mag[ceni['counts'] > 10] = good_mag
            slot_gspr = [pr for pr in self.gspr_info
                          if pr['slot'] == star['slot']]
            if not slot_gspr:
                err = "No GSPR info found for slot %d" % star['slot']
                raise ValueError(err)
            self.deltas[star['slot']]= dict(dy=dy,
                                            dz=dz,
                                            dr=dr,
                                            yag=yag,
                                            zag=zag,
                                            time=ceni['time'],
                                            mag=mag,
                                            qual=ceni['status'],
                                            ang_y_sm=ceni['ang_y_sm'],
                                            ang_z_sm=ceni['ang_z_sm'],
                                            ang_y=ceni['ang_y'],
                                            ang_z=ceni['ang_z'],
                                            )

    def _calc_sim_offset(self):
        mm2a = 20.0
        abs_sim_dy0 = 10.0
        abs_sim_dz0 = 10.0
        n_med = 21
        medf_dy = medfilt(self.asol.dy * mm2a, kernel_size=n_med)
        medf_dz = medfilt(self.asol.dz * mm2a, kernel_size=n_med)
        d_dy = abs(medf_dy[1:] - medf_dy[:-1])
        d_dz = abs(medf_dz[1:] - medf_dz[:-1])
        self.sim = dict(time=self.asol['time'][:-1],
                        medf_dy=medf_dy,
                        medf_dz=medf_dz,
                        d_dy=d_dy,
                        d_dz=d_dz)



        
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
