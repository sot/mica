from __future__ import division
import os
import sys
import re
import logging
import gzip
import jinja2
import datetime
from glob import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.io import fits

from Ska.Shell import getenv, tcsh, bash
import agasc
import Ska.DBI
from Chandra.Time import DateTime

from mica.archive import obspar
from mica.catalog import catalog
from mica.starcheck import starcheck
from mica.archive import asp_l1
import mica.vv
from mica.vv import get_vv, get_vv_files, get_arch_vv
from mica.version import version

WANT_VV_VERSION = 2

plt.rcParams['lines.markeredgewidth'] = 0

logger = logging.getLogger('mica.report')
logger.setLevel(logging.INFO)
if not len(logger.handlers):
    logger.addHandler(logging.StreamHandler())

aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')
DEFAULT_REPORT_ROOT = "/proj/web-icxc/htdocs/aspect/mica_reports"
DAILY_PLOT_ROOT="http://occweb.cfa.harvard.edu/occweb/FOT/engineering/reports/dailies"

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
    starcheck_html="{top}{mp_dir}starcheck.html#obsid{obsid}".format(
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
    acqs = aca_db.fetchall("""select * from acq_stats_data
                              where agasc_id = %d
                              order by tstart""" % int(id))
    if not len(acqs):
        return None, dict(n_acqs=0, n_acq_noid=0, avg_mag=None)
    n_acqs = len(acqs)
    n_acq_noid = len(np.flatnonzero(acqs['obc_id'] == 'NOID'))
    avg_mag = np.mean(acqs[acqs['mag_obs'] < 13.94]['mag_obs'])
    # make a list of dictionaries to make it easy to add values to the rows
    acqs = rec_to_dict_list(acqs)
    for acq in acqs:
        acq['sc_link'] = '<A HREF="{link}">{obsid}</A>'.format(
            link=starcheck_link(acq['obsid']),
            obsid=acq['obsid'])
    return acqs, dict(n_acqs=n_acqs, n_acq_noid=n_acq_noid, avg_mag=avg_mag)


def get_star_trak_stats(id):
    traks = aca_db.fetchall("""select * from trak_stats_data
                               where id = %d
                               order by kalman_tstart""" % int(id))
    if not len(traks):
        return None, dict(n_guis=0, n_bad=0, n_fail=0, n_obc_bad=0, avg_mag=None)

    n_guis = len(traks)
    n_bad = len(np.flatnonzero((traks['not_tracking_samples'] / traks['n_samples']) > 0.05))
    n_fail = len(np.flatnonzero((traks['not_tracking_samples'] / traks['n_samples']) == 1))
    n_obc_bad = len(np.flatnonzero((traks['obc_bad_status_samples'] / traks['n_samples']) > 0.05))
    avg_mag = np.mean(traks[traks['aoacmag_mean'] < 13.94]['aoacmag_mean'])

    # make a list of dictionaries to make it easy to add values to the rows
    traks = rec_to_dict_list(traks)
    mytraks = []
    for trak in traks:
        star = dict(
            obsid=trak['obsid'],
            tstart="{:11.1f}".format(trak['kalman_tstart']),
            mag_obs="{:6.3f}".format(trak['aoacmag_median']),
            mag_obs_rms = "{:.3f}".format(trak['aoacmag_rms']),
            trak_frac="{:.1f}".format((100 * ((trak['n_samples']
                               - trak['not_tracking_samples'])
                              ) / trak['n_samples'])))
        for stat in ['obc_bad_status',
                     'def_pix', 'ion_rad', 'sat_pix',
                     'mult_star', 'quad_bound', 'common_col']:
            star.update({stat: "{:.3f}".format(
                        ((100 * trak["{}_samples".format(stat)])
                         / trak['n_samples']))})
        star['sc_link'] = '<A HREF="{link}">{obsid}</A>'.format(
            link=starcheck_link(trak['obsid']),
            obsid=trak['obsid'])
        mytraks.append(star)
    #if id == 896534664:
    #    raise ValueError
    return mytraks, dict(n_guis=n_guis, n_bad=n_bad, n_fail=n_fail, n_obc_bad=n_obc_bad, avg_mag=avg_mag)


def get_obs_acq_stats(obsid):
    acqs = aca_db.fetchall("select * from acq_stats_data where obsid = %d" % obsid)
    if len(acqs):
        return Table(acqs)


def get_obs_trak_stats(obsid):
    guis = aca_db.fetchall("select * from trak_stats_data where obsid = %d" % obsid)

    if len(guis):
        return Table(guis)


def target_summary(obsid):
    ocat_db = Ska.DBI.DBI(dbi='sybase', server='sqlsao', user='aca_ops', database='axafocat')
    ocat_info = ocat_db.fetchone("""select * from target inner join prop_info on
                                    target.proposal_id = prop_info.proposal_id
                                    and target.obsid = {}""".format(obsid))
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


def obs_links(obsid, sequence=None):
    mp_dir, status, mp_date = starcheck.get_mp_dir(obsid)
    links = dict(
        obscat=dict(
            label="Target Param : {}".format(obsid),
            link="https://icxc.cfa.harvard.edu/cgi-bin/mp/target_param.cgi?{}".format(obsid)),
        mp_dir=None,
        shortterm=None,
        fot_dir=None,
        fot_daily=None,
        starcheck_html=None,
        vv=None)
    if sequence is not None:
        links.update(dict(
                seq_sum=dict(
                    label="Seq. Summary : {}".format(sequence),
                    link= "https://icxc.harvard.edu/cgi-bin/mp/target.cgi?{}".format(sequence))))
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
            links.update(
                dict(fot_daily=dict(
                    label="Daily Plots {}:{}".format(year, doy),
                    link="{root}/{year}/{upper_month}/{lower_month}{day_of_month}_{doy}/".format(
                        root=DAILY_PLOT_ROOT, year=year, upper_month=month.upper(), 
                        lower_month=month.lower(), day_of_month=dom, doy=doy))))
        vv, vvid = official_vv(obsid)
        links.update(dict(
                shortterm=dict(link=guess_shortterm(mp_dir),
                               label="Short Term Sched. {}".format(mp_label)),
                fot_dir=dict(link=guess_fot_summary(mp_dir),
                             label="FOT Approved Sched. {}".format(mp_label)),
                starcheck_html=dict(
                    link="{top}{mp_dir}starcheck.html#obsid{obsid}".format(
                        top="https://icxc.harvard.edu/mp/mplogs",
                        obsid=obsid,
                        mp_dir=mp_dir),
                    label="Starcheck obsid {}".format(obsid)),
                vv=dict(link=vv,
                        label="CXCDS V&V (id={})".format(vvid))))
    return links


def catalog_info(starcheck_cat, acqs=None, trak=None, vv=None):
    table = []
    sc = ['idx', 'slot', 'id', 'type', 'sz', 'mag', 'yang', 'zang']
    for row in starcheck_cat:
        sc_row = dict([(i, row[i]) for i in sc])
        sc_row.update({"pass_notes": "{}{}".format(row['pass'] or '',
                                                   row['notes'] or '')})
        table.append(sc_row)

    if acqs is not None and len(acqs):
        for row in acqs:
            for sc_row in table:
                if ((sc_row['slot'] == row['slot'])
                    and ((sc_row['type'] == 'ACQ')
                         or (sc_row['type'] == 'BOT'))):
                    sc_row.update(dict(mag_obs=row['mag_obs'],
                                       obc_id=row['obc_id'],
                                       mag_source='acq'))
    if trak is not None and len(trak):
        for row in trak:
            for sc_row in table:
                if ((sc_row['slot'] == row['slot'])
                    and (sc_row['type'] != 'ACQ')):
                    sc_row.update(dict(mag_obs=row['aoacmag_median'],
                                       trak_frac=(100 * ((row['n_samples']
                                                          - row['not_tracking_samples'])
                                                         ) / row['n_samples']),
                                       obc_bad_frac=(100 * ((row['obc_bad_status_samples']
                                                             / row['n_samples']))),
                                       mag_source='trak'))
    if vv is not None:
        for sc_row in table:
            if sc_row['type'] != 'ACQ':
                slot = sc_row['slot']
                if str(slot) in vv['slots']:
                    vvslot = vv['slots'][str(slot)]
                    sc_row.update(dict(dr_rms=vvslot['dr_rms'],
                                       dz_rms=vvslot['dz_rms'],
                                       dy_rms=vvslot['dy_rms'],
                                       dy_mean=vvslot['dy_mean'],
                                       dz_mean=vvslot['dz_mean'],
                                       id_status=vvslot.get('id_status'),
                                       cel_loc_flag=vvslot.get('cel_loc_flag')))


    # let's explicitly reformat everything
    format_spec = dict(
        idx="{:d}",
        slot="{:d}",
        id="{:d}",
        type="{}",
        sz="{}",
        mag="{:6.3f}",
        yang="{:d}",
        zang="{:d}",
        pass_notes="&nbsp;{}",
        obc_id="{}",
        mag_obs="{:6.3f}",
        trak_frac="{:.0f}",
        obc_bad_frac="{:.0f}",
        mag_source="{}",
        id_status="{}",
        cel_loc_flag="{}",
        dr_rms="{:6.3f}",
        dz_rms="{:6.3f}",
        dy_rms="{:6.3f}",
        dy_mean="{:6.3f}",
        dz_mean="{:6.3f}",
        )

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
    agasc_info = agasc.get_star(id)
    agasc_list = [(key, agasc_info[key]) for key in agasc_info.dtype.names]
    return agasc_list


def star_info(id):
    agasc_info = get_star(id)
    acqs, agg_acq = get_star_acq_stats(id)
    traks, agg_trak = get_star_trak_stats(id)
    return dict(agasc_info=agasc_info,
                acqs=acqs,
                agg_acq=agg_acq,
                traks=traks,
                agg_trak=agg_trak)

def get_aiprops(obsid):
    aiprops = aca_db.fetchall(
        "select * from aiprops where obsid = {} order by tstart".format(
            obsid))
    return aiprops

def main(obsid, report_root=DEFAULT_REPORT_ROOT):

    strobs = "%05d" % obsid
    chunk_dir = strobs[0:2]
    topdir = os.path.join(report_root, chunk_dir)
    outdir = os.path.join(topdir, strobs)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    jinja_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader('templates'))
    jinja_env.line_comment_prefix = '##'
    jinja_env.line_statement_prefix = '#'

    logger.info("Making report for {}".format(obsid))
    logger.info("Getting target info from axafapstat")
    summary = target_summary(obsid)
    # Links
    logger.info("Looking up obsid links")

    all_progress = dict(science=['ocat',
                             'long_term', 'short_term',
                             'starcheck', 'observed',
                             'engineering', 'aspect_1',
                             'cxcds_vv', 'released'],
                    er=['starcheck', 'observed',
                         'engineering', 'cxcds_vv'])
    report_status = dict()

    er = summary is None and obsid > 40000
    progress = all_progress['er' if er else 'science']
    if er:
        links = obs_links(obsid)
    else:
        if summary is None:
            raise ValueError("Obsid not found in target table")
        report_status['ocat'] = summary['status']
        links = obs_links(obsid, summary['seq_nbr'])

    if not er and (summary['status'] in
                   ['canceled', 'unobserved', 'untriggered']):
        logger.info(
            "Obsid {obsid} has status {status}".format(
                obsid=obsid, status=summary['status']))

    if summary is not None:
        if summary['lts_lt_plan'] is not None:
            report_status['long_term'] = summary['lts_lt_plan']
            plan = summary['lts_lt_plan']
            plan_date = DateTime("{:4d}:{:03d}:{:02d}:{:02d}:{:2d}.000".format(
                    plan.year, plan.dayofyear, plan.hour, plan.minute, plan.second))
            now = DateTime()
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
        obs_sc, mp_dir, status = get_starcheck(obsid)
        logger.info("Plotting starcheck catalog to {}".format(os.path.join(outdir, 'starcheck.png')))
        if obs_sc['obs'][0]['point_ra'] is None:
            raise LookupError("Observation has no pointing.")
        fig, cat, obs = catalog.plot(obsid, mp_dir)
        fig.savefig(os.path.join(outdir, 'starcheck.png'))
        plt.close('all')
    except LookupError:
        logger.info("No starcheck.  Writing out OCAT info only")
        template = jinja_env.get_template('report.html')
        page = template.render(obsid=obsid,
                               target=summary,
                               links=links,
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

    if acqs or trak:
        last_sched = "eng. data available"

    er_status = None
    if er:
        stat_map = dict(ran='ran on',
                        approved='approved for',
                        ran_pretimelines='ran on',
                        planned='planned for')
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
        vv['has_errors'] = ('errors' in vv) or None
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
                cat_row['idlink'] = '<A HREF="star_{id}.html" STYLE="text-decoration: none;" ONMOUSEOVER="return overlib (\'ACQ total:{n_acq} noid:{n_noid} <BR /> GUI total:{n_gui} bad:{n_bad} fail:{n_fail} obc_bad:{n_obc_bad} <BR /> Avg Mag {avg_mag:4.2f}\', WIDTH, 220);", ONMOUSEOUT="return nd();"> {id}</A>'.format(id=int(row['id']), n_acq=s['agg_acq']['n_acqs'], n_noid=s['agg_acq']['n_acq_noid'], n_gui=s['agg_trak']['n_guis'], n_bad=s['agg_trak']['n_bad'], n_fail=s['agg_trak']['n_fail'], n_obc_bad=s['agg_trak']['n_obc_bad'], avg_mag=s['agg_trak']['avg_mag'] or s['agg_acq']['avg_mag'] or 13.94)
            else:
                cat_row['idlink'] = "&nbsp;"
        else:
            cat_row['idlink'] = int(row['id'])

    template = jinja_env.get_template('report.html')
    page = template.render(cat_table=cat_table,
                           obs=obs,
                           sc=obs_sc,
                           vv=vv,
                           links=links,
                           target=summary,
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



if __name__ == '__main__':
    opt = get_options()
    main(opt.obsid, opt.report_root)
