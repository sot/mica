from __future__ import division
import os
import sys
import warnings
# Ignore compiler warnings that seem to be coming from a django.db
# interaction (via kadi.events)
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    import compiler
import re
import logging
import gzip
import jinja2
import datetime
import json
from glob import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from Ska.Shell import bash
import agasc
import Ska.DBI
from Chandra.Time import DateTime
from Ska.engarchive import fetch_sci
from kadi import events

from mica.archive import obspar
from mica.catalog import catalog
from mica.starcheck import starcheck
from mica.archive import asp_l1
import mica.vv
from mica.vv import get_vv, get_vv_files
from mica.version import version
from mica.common import MICA_ARCHIVE

WANT_VV_VERSION = 2
REPORT_VERSION = 2

plt.rcParams['lines.markeredgewidth'] = 0

logger = logging.getLogger('mica.report')
logger.setLevel(logging.INFO)
if not len(logger.handlers):
    logger.addHandler(logging.StreamHandler())

ACA_DB = None
DEFAULT_REPORT_ROOT = "/proj/web-icxc/htdocs/aspect/mica_reports"
DAILY_PLOT_ROOT = "http://occweb.cfa.harvard.edu/occweb/FOT/engineering/reports/dailies"

FILES = {'sql_def': 'report_processing.sql'}
DEFAULT_CONFIG = {
    'dbi': 'sqlite',
    'report_root': '/proj/web-icxc/htdocs/aspect/mica_reports',
    'server': os.path.join(MICA_ARCHIVE, 'report', 'report_processing.db3'),
    'update_mode': True,
    'retry_failure': False}


def get_options():
    parser = argparse.ArgumentParser(
        description="Create a mica report for an obsid")
    parser.add_argument("obsid",
                        help="obsid",
                        type=int)
    parser.add_argument("--report-root",
                        default=DEFAULT_REPORT_ROOT)
    opt = parser.parse_args()
    return opt


def get_starcheck(obsid):
    mp_dir, status, mp_date = starcheck.get_mp_dir(obsid)
    sc = starcheck.obsid(obsid, mp_dir)
    return (sc, mp_dir, status)


def starcheck_orig_link(obsid):
    """
    Return a link to the original starcheck products for an obsid
    """
    mp_dir, status, mp_date = starcheck.get_mp_dir(obsid)
    if mp_dir is None:
        return None
    starcheck_html = "{top}{mp_dir}starcheck.html#obsid{obsid}".format(
        top="https://icxc.harvard.edu/mp/mplogs",
        mp_dir=mp_dir,
        obsid=str(obsid))
    return starcheck_html


def starcheck_link(obsid):
    """
    Construct/return a link to the starcheck printer CGI.
    Making these kind of links scales well, as there are no database lookup required for each
    obsid until the link is clicked.
    """
    icxc_cgi_bin = "https://icxc.harvard.edu/cgi-bin/aspect/"
    starcheck_html = "{top}starcheck_print/starcheck_print.cgi?sselect=obsid;obsid1={obsid}".format(
        obsid=obsid,
        top=icxc_cgi_bin)
    return starcheck_html


def rec_to_dict_list(recarray):
    return [dict(zip(recarray.dtype.names, r)) for r in recarray]


def get_star_acq_stats(id):
    acqs = ACA_DB.fetchall("""select * from acq_stats_data
                              where agasc_id = %d
                              order by tstart""" % int(id))
    if not len(acqs):
        return None, {'n_acqs': 0, 'n_acq_noid': 0, 'avg_mag': None}
    n_acqs = len(acqs)
    n_acq_noid = len(np.flatnonzero(acqs['obc_id'] == 'NOID'))
    avg_mag = np.mean(acqs[acqs['mag_obs'] < 13.94]['mag_obs'])
    # make a list of dictionaries to make it easy to add values to the rows
    acqs = rec_to_dict_list(acqs)
    for acq in acqs:
        acq['sc_link'] = '<A HREF="{link}">{obsid}</A>'.format(
            link=starcheck_link(acq['obsid']),
            obsid=acq['obsid'])
    return acqs, {'n_acqs': n_acqs, 'n_acq_noid': n_acq_noid, 'avg_mag': avg_mag}


def get_star_trak_stats(id):
    traks = ACA_DB.fetchall("""select * from trak_stats_data
                               where id = %d
                               order by kalman_tstart""" % int(id))
    if not len(traks):
        return None, {'n_guis': 0, 'n_bad': 0, 'n_fail': 0,
                      'n_obc_bad': 0, 'avg_mag': None}

    n_guis = len(traks)
    n_bad = len(np.flatnonzero((traks['not_tracking_samples'] / traks['n_samples']) > 0.05))
    n_fail = len(np.flatnonzero((traks['not_tracking_samples'] / traks['n_samples']) == 1))
    n_obc_bad = len(np.flatnonzero((traks['obc_bad_status_samples'] / traks['n_samples']) > 0.05))
    avg_mag = np.mean(traks[traks['aoacmag_mean'] < 13.94]['aoacmag_mean'])

    # make a list of dictionaries to make it easy to add values to the rows
    traks = rec_to_dict_list(traks)
    mytraks = []
    for trak in traks:
        star = {
            'obsid': trak['obsid'],
            'tstart': "{:11.1f}".format(trak['kalman_tstart']),
            'mag_obs': "{:6.3f}".format(trak['aoacmag_median']),
            'mag_obs_rms': "{:.3f}".format(trak['aoacmag_rms']),
            'trak_frac': "{:.1f}".format(
                (100 * ((trak['n_samples']
                         - trak['not_tracking_samples'])
                        ) / trak['n_samples']))}
        for stat in ['obc_bad_status',
                     'def_pix', 'ion_rad', 'sat_pix',
                     'mult_star', 'quad_bound', 'common_col']:
            star[stat] = "{:.3f}".format(
                ((100 * trak["{}_samples".format(stat)])
                 / trak['n_samples']))
        star['sc_link'] = '<A HREF="{link}">{obsid}</A>'.format(
            link=starcheck_link(trak['obsid']),
            obsid=trak['obsid'])
        mytraks.append(star)
    return mytraks, {'n_guis': n_guis, 'n_bad': n_bad,
                     'n_fail': n_fail, 'n_obc_bad': n_obc_bad,
                     'avg_mag': avg_mag}


def get_obs_acq_stats(obsid):
    acqs = ACA_DB.fetchall("select * from acq_stats_data where obsid = %d" % obsid)
    if len(acqs):
        return Table(acqs)


def get_obs_trak_stats(obsid):
    guis = ACA_DB.fetchall("select * from trak_stats_data where obsid = %d" % obsid)

    if len(guis):
        return Table(guis)

def get_obs_temps(obsid, outdir):
    try:
        manvrs = events.manvrs.filter(obsid=obsid)
        dwells = events.dwells.filter(obsid=obsid)
    except ValueError:
        return None
    if len(manvrs) and len(dwells):
        ccd_temp = fetch_sci.MSID('AACCCDPT', manvrs[0].stop, dwells[0].stop)
        if len(ccd_temp.vals) == 0:
            return None
        return {'max': ccd_temp.vals.max(),
                'mean': ccd_temp.vals.mean()}



def target_summary(obsid):
    ocat_db = Ska.DBI.DBI(dbi='sybase', server='sqlsao', user='aca_ops', database='axafocat')
    ocat_info = ocat_db.fetchone("""select * from target inner join prop_info on
                                    target.proposal_id = prop_info.proposal_id
                                    and target.obsid = {}""".format(obsid))
    # If this target didn't have a proposal, just get whatever is there
    if ocat_info is None:
        ocat_info = ocat_db.fetchone("""select * from target where
                                    target.obsid = {}""".format(obsid))

    ocat_db.conn.close()
    del ocat_db
    return ocat_info


def guess_shortterm(mp_dir):
#file:///proj/web-icxc/htdocs/mp/schedules/cycle14/JUL2213B.html
    mp_short_top = '/proj/web-icxc/htdocs/mp/schedules'
    dir_match = re.match('/\d{4}/(\w{3}\d{4})/ofls(\w)', mp_dir)
    if not dir_match:
        return None
    dir_string = "{}{}.html".format(
        dir_match.group(1), dir_match.group(2).upper())
    shortterm = glob(os.path.join(mp_short_top,
                                  "cycle*",
                                  dir_string))
    if not len(shortterm):
        return None
    return re.sub('/proj/web-icxc/htdocs', 'https://icxc.harvard.edu', shortterm[0])


def guess_fot_summary(mp_dir):
    dir_match = re.match('/(\d{4})/((\w{3})\d{4})/ofls(\w)', mp_dir)
    if not dir_match:
        return None
    appr_loads = 'http://occweb.cfa.harvard.edu/occweb/FOT/mission_planning/PRODUCTS/APPR_LOADS'
    year = dir_match.group(1)
    week = dir_match.group(2) + dir_match.group(4).upper()
    mon = dir_match.group(3)
    dir_string = "{top}/{year}/{mon}/{week}/{week}.html".format(
        top=appr_loads, year=year, week=week, mon=mon)
    return dir_string
#http://occweb.cfa.harvard.edu/occweb/FOT/mission_planning/PRODUCTS/APPR_LOADS/2013/JUL/JUL1813A/JUL1813A.html


def official_vv(obsid):
    vv_db = Ska.DBI.DBI(dbi='sybase', server='sqlsao', user='jeanconn', database='axafvv')
    vv = vv_db.fetchone("""select vvid from vvreport where obsid = {obsid}
                           and creation_date = (
                              select max(creation_date) from vvreport where obsid = {obsid})
                        """.format(obsid=obsid))
    if vv is not None:
        vv_url = "https://icxc.cfa.harvard.edu/cgi-bin/vv/vv_report.cgi?vvid={}".format(vv['vvid'])
        del vv_db
        return vv_url, vv['vvid']
    else:
        del vv_db
        return None, None



def official_vv_notes(obsid):
    vv_db = Ska.DBI.DBI(dbi='sybase', server='sqlsao', user='jeanconn', database='axafvv',
                        numpy=False)
    all_vv = vv_db.fetchall("""select * from vvreport where obsid = {obsid}
                        """.format(obsid=obsid))
    if not len(all_vv):
        return None
    for report in all_vv:
        aspect_rev = vv_db.fetchone("""select * from vvreview where vvid = {vvid}
                        """.format(vvid=report['vvid']))
        report['aspect_review'] = aspect_rev

    del vv_db
    return all_vv


def obs_links(obsid, sequence=None, plan=None):
    links = {
        'obscat': {
            'label': "Target Param : {}".format(obsid),
            'link': "https://icxc.cfa.harvard.edu/cgi-bin/mp/target_param.cgi?{}".format(obsid)},
        'mp_dir': None,
        'shortterm': None,
        'fot_dir': None,
        'fot_daily': None,
        'starcheck_html': None,
        'vv': None}
    if sequence is not None:
        links['seq_sum'] = {
            'label': "Seq. Summary : {}".format(sequence),
            'link': "https://icxc.harvard.edu/cgi-bin/mp/target.cgi?{}".format(sequence)}
    mp_dir = None
    # if this is a science observation, only try to get a star catalog if it has a home
    # in the schedule either in the past or the near future
    if plan is not None:
        plan_date = DateTime("{:4d}:{:03d}:{:02d}:{:02d}:{:02d}.000".format(
                plan.year, plan.dayofyear, plan.hour, plan.minute, plan.second))
        if plan_date.secs < (DateTime() + 21).secs:
            mp_dir, status, mp_date = starcheck.get_mp_dir(obsid)
    else:
        mp_dir, status, mp_date = starcheck.get_mp_dir(obsid)

    if mp_dir is not None:
        dir_match = re.match('/\d{4}/(\w{3}\d{4})/ofls(\w)', mp_dir)
        mp_label = "{}{}".format(dir_match.group(1),
                                 dir_match.group(2).upper())
        mp_date_m = re.match('(\d{4}):(\d{3}):\d{2}:\d{2}:\d{2}\.\d{3}', mp_date)
        if mp_date_m:
            year = int(mp_date_m.group(1))
            doy = int(mp_date_m.group(2))

            dtm = datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)
            month = dtm.strftime("%b")
            dom = dtm.strftime("%d")
            links['fot_daily'] = {
                'label': "Daily Plots {}:{:03d}".format(year, int(doy)),
                'link': "{root}/{year}/{upper_month}/{lower_month}{day_of_month:02d}_{doy:03d}/".format(
                    root=DAILY_PLOT_ROOT, year=year, upper_month=month.upper(),
                    lower_month=month.lower(), day_of_month=int(dom), doy=int(doy))}
        vv, vvid = official_vv(obsid)
        links['shortterm'] = {'link': guess_shortterm(mp_dir),
                              'label': "Short Term Sched. {}".format(mp_label)}
        links['fot_dir'] = {'link': guess_fot_summary(mp_dir),
                            'label': "FOT Approved Sched. {}".format(mp_label)}
        links['starcheck_html'] = {'link': "{top}{mp_dir}starcheck.html#obsid{obsid}".format(top="https://icxc.harvard.edu/mp/mplogs", obsid=obsid, mp_dir=mp_dir),
                                   'label': "Starcheck obsid {}".format(obsid)}
        links['mp_dir'] = {'link': mp_dir,
                           'label': "Starcheck obsid {}".format(obsid)}
        links['vv'] = {'link': vv,
                       'label': "CXCDS V&V (id={})".format(vvid)}
    return links


def catalog_info(starcheck_cat, acqs=None, trak=None, vv=None):
    table = []
    sc = ['idx', 'slot', 'id', 'type', 'sz', 'mag', 'yang', 'zang']
    for row in starcheck_cat:
        sc_row = dict([(i, row[i]) for i in sc])
        sc_row['pass_notes'] = "{}{}".format(row['pass'] or '',
                                             row['notes'] or '')
        table.append(sc_row)

    if acqs is not None and len(acqs):
        for row in acqs:
            for sc_row in table:
                if ((sc_row['slot'] == row['slot'])
                    and ((sc_row['type'] == 'ACQ')
                         or (sc_row['type'] == 'BOT'))):
                    sc_row['mag_obs'] = row['mag_obs']
                    sc_row['obc_id'] = row['obc_id']
                    sc_row['mag_source'] = 'acq'
    if trak is not None and len(trak):
        for row in trak:
            for sc_row in table:
                if ((sc_row['slot'] == row['slot'])
                    and (sc_row['type'] != 'ACQ')):
                    sc_row['mag_obs'] = row['aoacmag_median']
                    sc_row['trak_frac'] = (100 * ((row['n_samples']
                                                  - row['not_tracking_samples'])
                                                  ) / row['n_samples'])
                    sc_row['obc_bad_frac'] = (100 * ((row['obc_bad_status_samples']
                                                      / row['n_samples'])))
                    sc_row['mag_source'] = 'trak'
    if vv is not None:
        for sc_row in table:
            if sc_row['type'] != 'ACQ':
                slot = sc_row['slot']
                if (str(slot) in vv['slots']) and ('n_pts' in vv['slots'][str(slot)]):
                    vvslot = vv['slots'][str(slot)]
                    sc_row['dr_rms'] = vvslot['dr_rms']
                    sc_row['dz_rms'] = vvslot['dz_rms']
                    sc_row['dy_rms'] = vvslot['dy_rms']
                    sc_row['dy_mean'] = vvslot['dy_mean']
                    sc_row['dz_mean'] = vvslot['dz_mean']
                    sc_row['id_status'] = vvslot.get('id_status')
                    sc_row['cel_loc_flag'] = vvslot.get('cel_loc_flag')


    # let's explicitly reformat everything
    format_spec = {
        'idx': "{:d}",
        'slot': "{:d}",
        'id': "{:d}",
        'type': "{}",
        'sz': "{}",
        'mag': "{:6.3f}",
        'yang': "{:d}",
        'zang': "{:d}",
        'pass_notes': "&nbsp;{}",
        'obc_id': "{}",
        'mag_obs': "{:6.3f}",
        'trak_frac': "{:.0f}",
        'obc_bad_frac': "{:.0f}",
        'mag_source': "{}",
        'id_status': "{}",
        'cel_loc_flag': "{}",
        'dr_rms': "{:6.3f}",
        'dz_rms': "{:6.3f}",
        'dy_rms': "{:6.3f}",
        'dy_mean': "{:6.3f}",
        'dz_mean': "{:6.3f}",
        }

    opt_elem = ['mag_obs', 'obc_id', 'trak_frac', 'obc_bad_frac', 'mag_source',
                'dr_rms', 'dz_rms', 'dy_rms', 'dy_mean', 'dz_mean',
                'pass_notes', 'id_status', 'cel_loc_flag']
    sc.extend(opt_elem)
    for sc_row in table:
        for k in sc:
            if k in sc_row and sc_row[k] is not None:
                sc_row[k] = format_spec[k].format(sc_row[k])
            else:
                sc_row[k] = '&nbsp;';
    return table


def get_star(id):
    try:
        agasc_info = agasc.get_star(id)
        agasc_list = [(key, agasc_info[key]) for key in agasc_info.dtype.names]
        return agasc_list
    except:
        return []


def star_info(id):
    agasc_info = get_star(id)
    acqs, agg_acq = get_star_acq_stats(id)
    traks, agg_trak = get_star_trak_stats(id)
    return {'agasc_info': agasc_info,
            'acqs': acqs,
            'agg_acq': agg_acq,
            'traks': traks,
            'agg_trak': agg_trak}

def get_aiprops(obsid):
    aiprops = ACA_DB.fetchall(
        "select * from aiprops where obsid = {} order by tstart".format(
            obsid))
    return aiprops

def main(obsid, config=None, report_root=None):

    if config is None:
         config = DEFAULT_CONFIG
    if report_root is None:
        report_root=config['report_root']

    global ACA_DB
    if ACA_DB is None or not ACA_DB.conn._is_connected:
        ACA_DB = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')


    strobs = "%05d" % obsid
    chunk_dir = strobs[0:2]
    topdir = os.path.join(report_root, chunk_dir)
    outdir = os.path.join(topdir, strobs)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    jinja_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(
            os.path.join(os.environ['SKA'], 'data', 'mica', 'templates')))
    jinja_env.line_comment_prefix = '##'
    jinja_env.line_statement_prefix = '#'

    logger.info("Making report for {}".format(obsid))
    logger.info("Getting target info from axafapstat")
    summary = target_summary(obsid)
    # Links
    logger.info("Looking up obsid links")

    all_progress = {'science': ['ocat',
                                 'long_term', 'short_term',
                                 'starcheck', 'observed',
                                 'engineering', 'aspect_1',
                                 'cxcds_vv', 'released'],
                    'er': ['starcheck', 'observed',
                           'engineering', 'cxcds_vv']}
    report_status = {}

    er = summary is None and obsid > 40000
    progress = all_progress['er' if er else 'science']
    if er:
        links = obs_links(obsid)
    else:
        if summary is None:
            raise ValueError("Obsid not found in target table")
        report_status['ocat'] = summary['status']
        links = obs_links(obsid, summary['seq_nbr'], summary['lts_lt_plan'])

    if not er and (summary['status'] in
                   ['canceled', 'unobserved', 'untriggered']):
        logger.info(
            "Obsid {obsid} has status {status}".format(
                obsid=obsid, status=summary['status']))

    if summary is not None:
        if summary['lts_lt_plan'] is not None:
            report_status['long_term'] = summary['lts_lt_plan']
        if summary['soe_st_sched_date']:
            report_status['short_term'] =  summary['soe_st_sched_date']

    last_sched = ''
    if not er:
        if summary['lts_lt_plan']:
            last_sched = "in LTS for {}".format(
                str(summary['lts_lt_plan']))
        if summary['soe_st_sched_date']:
            last_sched = "in ST sched for {}".format(
                str(summary['soe_st_sched_date']))

    ## Starcheck
    logger.info("Fetching starcheck catalog")
    try:
        if summary is not None and summary['lts_lt_plan'] is not None:
            plan = summary['lts_lt_plan']
            plan_date = DateTime("{:4d}:{:03d}:{:02d}:{:02d}:{:02d}.000".format(
                    plan.year, plan.dayofyear, plan.hour, plan.minute, plan.second))
            if plan_date.secs > (DateTime() + 21).secs:
                raise LookupError("No starcheck expected for {} lts date".format(str(plan)))
        mp_dir, status, mp_date = starcheck.get_mp_dir(obsid)
        obs_sc, mp_dir, status = get_starcheck(obsid)
        logger.info("Plotting starcheck catalog to {}".format(os.path.join(outdir, 'starcheck.png')))
        if obs_sc['obs'][0]['point_ra'] is None:
            raise LookupError("Observation has no pointing.")
        if len(obs_sc['catalog']) == 0:
            raise LookupError("Observation has no catalog")
        fig, cat, obs = catalog.plot(obsid, mp_dir)
        sc = starcheck.get_starcheck_catalog(obsid, mp_dir)
        fig.savefig(os.path.join(outdir, 'starcheck.png'))
        plt.close('all')
    except LookupError as detail:
        logger.info("No starcheck catalog.  Writing out OCAT info only")
        logger.info(detail)
        template = jinja_env.get_template('report.html')
        page = template.render(obsid=obsid,
                               target=summary,
                               links=links,
                               temps=None,
                               pred_temp=None,
                               cat_table=None,
                               er=er if er else None,
                               last_sched=last_sched,
                               obs=None,
                               version=version)
        full_report_file = os.path.join(outdir, 'index.html')
        logger.info("Writing out full report to {}".format(full_report_file))
        f = open(full_report_file, 'w')
        f.write(page)
        f.close()
        notes = {'report_version': REPORT_VERSION,
                 'vv_version': None,
                 'vv_revision': None,
                 'aspect_1_id': None,
                 'last_sched': last_sched,
                 'ocat_status': report_status.get('ocat'),
                 'long_term': str(report_status.get('long_term')),
                 'short_term': str(report_status.get('short_term')),
                 'starcheck': report_status.get('starcheck'),
                 'obsid': obsid,
                 'checked_date': DateTime().date}
        f = open(os.path.join(outdir, 'notes.json'), 'w')
        f.write(json.dumps(notes,
                           sort_keys=True,
                           indent=4))
        f.close()
        save_state_in_db(obsid, notes, config)
        return

    if not er and 'shortterm' in links:
        dir_match = re.match('/\d{4}/(\w{3}\d{4})/ofls(\w)', mp_dir)
        mp_label = "{}{}".format(dir_match.group(1),
                                 dir_match.group(2).upper())
        last_sched = 'in <A HREF="{}">{}</A> at {}'.format(
            links['shortterm']['link'], mp_label, str(obs_sc['obs'][0]['mp_starcat_time']))


    report_status['starcheck'] = mp_dir

    # engineering data available
    logger.info("Getting acq and trak stats")
    acqs = get_obs_acq_stats(obsid)
    trak = get_obs_trak_stats(obsid)
    temps = get_obs_temps(obsid, outdir)
    pred_temp = sc['pred_temp']
    if acqs or trak:
        last_sched = "eng. data available"

    er_status = None
    if er:
        stat_map = {'ran': 'ran on',
                    'approved': 'approved for',
                    'ran_pretimelines': 'ran on',
                    'planned': 'planned for'}
        er_status = "{} {}".format(stat_map[status], obs_sc['obs'][0]['mp_starcat_time'])
        run_obspar = None
        vv = None
        logger.info("Processing ER; no V&V available")
    else:
        # obspar ingested
        try:
            run_obspar = obspar.get_obspar(obsid, version='last')
        except:
            run_obspar = None

        # v&v available
        try:
            vv = get_vv(obsid, version='last')
        except LookupError:
            vv = None

        try:
            if vv is None or 'vv_version' not in vv or vv['vv_version'] < WANT_VV_VERSION:
                mica.vv.process.process(obsid, version="last")
                vv = get_vv(obsid, version='last')
            vv_files = get_vv_files(obsid, version='last')
            last_sched = "complete through mica v&v"
        except LookupError:
            logger.info("No V&V available")
            vv = None

    if vv is not None:
        report_status['mica_vv'] = True
        for file in vv_files:
            newfile = os.path.join(outdir, os.path.basename(file))
            if not os.path.exists(newfile):
                logger.info("linking {} into {}".format(file, outdir))
                bash("ln -s {} {}".format(file, outdir))
        asp_dir = asp_l1.get_obs_dirs(obsid)['last']
        asp_logs = sorted(glob(os.path.join(asp_dir, "asp*log*gz")))
        for log, interval in zip(asp_logs, vv['intervals']):
            logmatch = re.search('(.*log)\.gz', os.path.basename(log))
            if logmatch:
                newlogname = "{}.txt".format(logmatch.group(1))
                newlog = os.path.join(outdir, newlogname)
                if not os.path.exists(newlog):
                    logger.info("copying/gunzipping asp log {}".format(newlog))
                    logtext = gzip.open(log).readlines()
                    f = open(newlog, 'w')
                    f.writelines(logtext)
                    f.close()
                interval['loglink'] = newlogname

        aiprops = get_aiprops(obsid)
        aiprops_template = jinja_env.get_template('aiprops.html')
        aiprops_page = aiprops_template.render(obsid=obsid, aiprops=aiprops)
        aiprops_page_file = os.path.join(outdir, 'aiprops.html')
        logger.info("AIPROPS report to {}".format(aiprops_page_file))
        f = open(aiprops_page_file, 'w')
        f.write(aiprops_page)
        f.close()

        props_template = jinja_env.get_template('props.html')
        props_page = props_template.render(obsid=obsid, vv=vv)
        props_page_file = os.path.join(outdir, 'props.html')
        logger.info("GS/FIDPROPS report to {}".format(props_page_file))
        f = open(props_page_file, 'w')
        f.write(props_page)
        f.close()

        for slot in vv['slots']:
            if 'n_pts' not in vv['slots'][slot]:
                continue
            slot_template = jinja_env.get_template('vv_slots_single.html')
            slot_page = slot_template.render(obsid=obsid,
                                             vv=vv,
                                             slot=slot)
            slot_page_file = os.path.join(outdir, "slot_{}.html".format(slot))
            logger.info("VV SLOT report to {}".format(slot_page_file))
            f = open(slot_page_file, 'w')
            f.write(slot_page)
            f.close()

        official_notes=official_vv_notes(obsid)
        vv_template = jinja_env.get_template('vv.html')
        vv['has_errors'] = (('errors' in vv) and (len(vv['errors']))) or None
        vv_page = vv_template.render(obsid=obsid,
                                     vv=vv,
                                     obspar=run_obspar,
                                     official_vv_notes=official_notes,
                                     )
        vv_page_file = os.path.join(outdir, 'vv.html')
        logger.info("VV report to {}".format(vv_page_file))
        f = open(vv_page_file, 'w')
        f.write(vv_page)
        f.close()

    cat_table = catalog_info(obs_sc['catalog'], acqs, trak, vv)

    for row, cat_row in zip(obs_sc['catalog'], cat_table):
        if row['type'] != 'FID':
            if row['id'] is not None:
                s = star_info(row['id'])
                star_template = jinja_env.get_template('star.html')
                star_page = star_template.render(
                    star=s['agasc_info'],
                    acqs=s['acqs'],
                    traks=s['traks'],
                    agg_acq=s['agg_acq'],
                    agg_trak=s['agg_trak'])
                star_page_file = os.path.join(outdir, 'star_%d.html' % int(row['id']))
                logger.info("Writing out star info to {}".format(star_page_file))
                f = open(star_page_file, 'w')
                f.write(star_page)
                f.close()
                cat_row['idlink'] = (
                    '<A HREF="star_{id}.html" STYLE="text-decoration: none;"'
                    'ONMOUSEOVER="return overlib '
                    '(\'ACQ total:{n_acq} noid:{n_noid} <BR /> '
                    'GUI total:{n_gui} bad:{n_bad} fail:{n_fail} obc_bad:{n_obc_bad} '
                    '<BR /> Avg Mag {avg_mag:4.2f}\', WIDTH, 220);", ONMOUSEOUT="return nd();"> '
                    '{id}</A>'.format(id=int(row['id']),
                                      n_acq=s['agg_acq']['n_acqs'],
                                      n_noid=s['agg_acq']['n_acq_noid'],
                                      n_gui=s['agg_trak']['n_guis'],
                                      n_bad=s['agg_trak']['n_bad'],
                                      n_fail=s['agg_trak']['n_fail'],
                                      n_obc_bad=s['agg_trak']['n_obc_bad'],
                                      avg_mag=(s['agg_trak']['avg_mag']
                                               or s['agg_acq']['avg_mag']
                                               or 13.94)))
            else:
                cat_row['idlink'] = "&nbsp;"
        else:
            if 'id' in row:
                cat_row['idlink'] = int(row['id'])
            else:
                cat_row['idlink'] = ''
    template = jinja_env.get_template('report.html')

    page = template.render(cat_table=cat_table,
                           obs=obs,
                           sc=obs_sc,
                           vv=vv,
                           links=links,
                           target=summary,
                           temps=temps,
                           pred_temp=pred_temp,
                           er=er if er else None,
                           er_status=er_status,
                           last_sched=last_sched,
                           obsid=obsid,
                           version=version)
    full_report_file = os.path.join(outdir, 'index.html')
    logger.info("Writing out full report to {}".format(full_report_file))
    f = open(full_report_file, 'w')
    f.write(page)
    f.close()

    cat_file = os.path.join(outdir, 'star_table.json')
    f = open(cat_file, 'w')
    f.write(json.dumps(cat_table, sort_keys=True, indent=4))
    f.close()

    notes = {'report_version': REPORT_VERSION,
             'vv_version': None,
             'vv_revision': None,
             'aspect_1_id': None,
             'last_sched': last_sched,
             'ocat_status': report_status.get('ocat'),
             'long_term': str(report_status.get('long_term')),
             'short_term': str(report_status.get('short_term')),
             'starcheck': report_status.get('starcheck'),
             'obsid': obsid,
             'checked_date': DateTime().date}
    if vv:
        notes['vv_version'] = vv.get('vv_version')
        notes['vv_revision'] = vv.get('revision')
        notes['aspect_1_id'] = vv.get('aspect_1_id')
    f = open(os.path.join(outdir, 'notes.json'), 'w')
    f.write(json.dumps(notes,
                       sort_keys=True,
                       indent=4))
    f.close()
    save_state_in_db(obsid, notes, config)


def save_state_in_db(obsid, notes, config=None):

    if config is None:
         config = DEFAULT_CONFIG

    if (not os.path.exists(config['server'])
        or os.stat(config['server']).st_size == 0):
        if not os.path.exists(os.path.dirname(config['server'])):
            os.makedirs(os.path.dirname(config['server']))
        db_sql = os.path.join(os.environ['SKA_DATA'], 'mica', FILES['sql_def'])
        db_init_cmds = open(db_sql).read()
        db = Ska.DBI.DBI(config['dbi'], config['server'])
        db.execute(db_init_cmds)
    else:
        db = Ska.DBI.DBI(dbi=config['dbi'], server=config['server'])

    notes['report_status'] = notes['last_sched']
    del notes['last_sched']

    idcheck = db.fetchone(
        "select * from report_proc "
        "where obsid = '{}' "
        "and report_version = '{}'".format(obsid,
                                           REPORT_VERSION))

    if idcheck is None:
        db.insert(notes, 'report_proc')
    else:
        db.execute("UPDATE report_proc SET "
                   "checked_date = '{checked_date}', "
                   "ocat_status = '{ocat_status}', "
                   "report_status = '{report_status}', "
                   "vv_version = '{vv_version}', "
                   "vv_revision = '{vv_revision}', "
                   "aspect_1_id = '{aspect_1_id}', "
                   "report_version = '{report_version}', "
                   "long_term = '{long_term}', "
                   "short_term = '{short_term}', "
                   "starcheck = '{starcheck}' "
                   "where obsid = '{obsid}' "
                   "and report_version = '{report_version}'".format(**notes))
    db.conn.close()


def update(report_root=DEFAULT_REPORT_ROOT):

    global ACA_DB
    if ACA_DB is None or not ACA_DB.conn._is_connected:
        ACA_DB = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')

    recent_obs = ACA_DB.fetchall("select distinct obsid from cmd_states "
                             "where datestart > '{}'".format(DateTime(-7).date))
    for obs in recent_obs:
        process_obsids([obs['obsid']], report_root=report_root)

    ACA_DB.conn.close()


def process_obsids(obsids, config=None, report_root=None):
    if config is None:
         config = DEFAULT_CONFIG
    if report_root is None:
        report_root=config['report_root']

    for obsid in obsids:
        strobs = "%05d" % obsid
        chunk_dir = strobs[0:2]
        topdir = os.path.join(report_root, chunk_dir)
        outdir = os.path.join(topdir, strobs)
        if os.path.exists(outdir) and not config['update_mode']:
            logger.info("Skipping {}, output dir exists.".format(obsid))
            continue
        if not config['retry_failure'] and os.path.exists(os.path.join(outdir, "proc_err")):
            logger.info("Skipping {}, previous processing error.".format(obsid))
            continue
        if not os.path.exists(outdir):
            os.makedirs("{}".format(outdir))
        # Delete files from old failure if reprocessing
        for failfile in ['proc_err', 'trace.txt']:
            if os.path.exists(os.path.join(outdir, failfile)):
                os.unlink(os.path.join(outdir, failfile))
        try:
            main(obsid, config=config, report_root=report_root)
        except:
            import traceback
            etype, emess, trace = sys.exc_info()
            logger.info("Failed report on {}".format(obsid))
            # Make an empty file to record the error status
            f = open(os.path.join(outdir, 'proc_err'), 'w')
            f.close()
            # Write out the traceback too
            trace_file = open(os.path.join(outdir, 'trace.txt'), 'w')
            traceback.print_tb(trace, file=trace_file)
            trace_file.close()
            # Write out a notes jason file
            notes = {'report_version': REPORT_VERSION,
                     'obsid': obsid,
                     'checked_date': DateTime().date,
                     'last_sched': "{}".format(str(emess)),
                     'vv_version': None,
                     'vv_revision': None,
                     'aspect_1_id': None,
                     'ocat_status': None,
                     'long_term': None,
                     'short_term': None,
                     'starcheck': None}
            f = open(os.path.join(outdir, 'notes.json'), 'w')
            f.write(json.dumps(notes,
                               sort_keys=True,
                               indent=4))
            f.close()
            # Make a stub html page
            proc_date = DateTime().date
            jinja_env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(
                    os.path.join(os.environ['SKA'], 'data', 'mica', 'templates')))
            jinja_env.line_comment_prefix = '##'
            jinja_env.line_statement_prefix = '#'
            template = jinja_env.get_template('proc_error.html')
            page = template.render(obsid=obsid,
                                   proc_date=proc_date,
                                   version=version)
            full_report_file = os.path.join(outdir, 'index.html')
            logger.info("Writing out error stub report to {}".format(full_report_file))
            f = open(full_report_file, 'w')
            f.write(page)
            f.close()
            # Save the bad state in the database
            save_state_in_db(obsid, notes, config=None)


def fill_first_time(report_root=DEFAULT_REPORT_ROOT):

    global ACA_DB
    if ACA_DB is None:
        ACA_DB = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')

    obs = ACA_DB.fetchall("select * from observations_all order by kalman_idstart desc")
    for ob in obs:
        obsid = ob['obsid']
        print ob['date']
        strobs = "%05d" % obsid
        chunk_dir = strobs[0:2]
        topdir = os.path.join(report_root, chunk_dir)
        outdir = os.path.join(topdir, strobs)
        if os.path.exists(outdir):
            logger.info("Skipping {}, output dir exist.".format(obsid))
        if os.path.exists("{}.ERR".format(outdir)):
            logger.info("Skipping {}, output.ERR dir exist.".format(obsid))
        if not os.path.exists(outdir):
            try:
                os.makedirs("{}".format(outdir))
                main(obsid, config=None, report_root=report_root)
            except:
                os.makedirs("{}.ERR".format(outdir))
                etype, emess, traceback = sys.exc_info()
                notes = {'report_version': REPORT_VERSION,
                         'obsid': obsid,
                         'checked_date': DateTime().date,
                         'last_sched': "{} {}".format(etype, emess),
                         'vv_version': None,
                         'vv_revision': None,
                         'aspect_1_id': None,
                         'ocat_status': None,
                         'long_term': None,
                         'short_term': None,
                         'starcheck': None}
                save_state_in_db(obsid, notes, config=None)

    ACA_DB.conn.close()

if __name__ == '__main__':
    opt = get_options()
    main(opt.obsid, config=None, report_root=opt.report_root)
