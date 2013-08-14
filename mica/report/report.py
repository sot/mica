from __future__ import division
import os
import sys
import re
import logging
import jinja2
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.io import fits

from Ska.Shell import getenv, tcsh
import agasc
import Ska.DBI

from mica.archive import obspar
from mica.catalog import catalog
from mica.starcheck import starcheck
from mica.archive import asp_l1
from mica.vv import get_vv, get_vv_files, get_arch_vv

plt.rcParams['lines.markeredgewidth'] = 0

logger = logging.getLogger('mica.report')
logger.setLevel(logging.INFO)
if not len(logger.handlers):
    logger.addHandler(logging.StreamHandler())

aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')


def get_starcheck(obsid):
    mp_dir, status = starcheck.get_mp_dir(obsid)
    sc = starcheck.obsid(obsid, mp_dir)
    return (sc, mp_dir, status)


def starcheck_orig_link(obsid):
    """
    Return a link to the original starcheck products for an obsid
    """
    mp_dir, status = starcheck.get_mp_dir(obsid)
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
        return None, None
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
        return None, None

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
    vv_url = None
    if vv is not None:
        vv_url = "https://icxc.cfa.harvard.edu/cgi-bin/vv/vv_report.cgi?vvid={}".format(vv['vvid'])
    del vv_db
    return vv_url

def obs_links(obsid, sequence=None):
    mp_dir, status = starcheck.get_mp_dir(obsid)
    links = dict(
        obscat="https://icxc.cfa.harvard.edu/cgi-bin/mp/target_param.cgi?{}".format(obsid),
        mp_dir=None,
        shortterm=None,
        fot_dir=None,
        starcheck_html=None,
        vv=None)
    if sequence is not None:
        links.update(dict(
                seq_sum="https://icxc.harvard.edu/cgi-bin/mp/target.cgi?{}".format(sequence)))
    if mp_dir is not None:
        links.update(dict(
                shortterm=guess_shortterm(mp_dir),
                fot_dir=guess_fot_summary(mp_dir),
                starcheck_html="{top}{mp_dir}starcheck.html#obsid{obsid}".format(
                    top="https://icxc.harvard.edu/mp/mplogs",
                    obsid=obsid,
                    mp_dir=mp_dir),
                vv=official_vv(obsid)))
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
                sc_row.update(dict(dr_rms=vv[str(slot)]['dr_rms'],
                                   dz_rms=vv[str(slot)]['dz_rms'],
                                   dy_rms=vv[str(slot)]['dy_rms'],
                                   dy_mean=vv[str(slot)]['dy_mean'],
                                   dz_mean=vv[str(slot)]['dz_mean'],
                                   id_status=vv[str(slot)].get('id_status'),
                                   cel_loc_flag=vv[str(slot)].get('cel_loc_flag')))


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
            if k in sc_row:
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


def main():

    outdir = 'test_report_new'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    jinja_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader('templates'))
    jinja_env.line_comment_prefix = '##'
    jinja_env.line_statement_prefix = '#'


    #obsid = 14932
    #obsid = 10980
    #obsid = 2306
    #obsid = 15315
    #obsid = 15293
    #obsid = 433
    #obsid = 15552
    #obsid = 14520
    #obsid = 53521
    #obsid = 15042
    #obsid = 14895
    obsid = 15129

    logger.info("Making report for {}".format(obsid))
    logger.info("Getting target info from axafapstat")
    summary = target_summary(obsid)
    # Links
    logger.info("Looking up obsid links")

    er = summary is None and obsid > 40000
    if er:
        links = obs_links(obsid)
    else:
        links = obs_links(obsid, summary['seq_nbr'])

    if not er and (summary['status'] in
                   ['canceled', 'unobserved', 'untriggered']):
        logger.info(
            "Obsid {obsid} has status {status}".format(
                obsid=obsid, status=summary['status']))

    ## Starcheck
    logger.info("Fetching starcheck catalog")
    try:
        obs_sc, mp_dir, status = get_starcheck(obsid)
        logger.info("Plotting starcheck catalog")
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
                               obs=None)
        full_report_file = os.path.join(outdir, 'report.html')
        logger.info("Writing out full report to {}".format(full_report_file))
        f = open(full_report_file, 'w')
        f.write(page)
        f.close()
        return


    
    # engineering data available
    logger.info("Getting acq and trak stats")
    acqs = get_obs_acq_stats(obsid)
    trak = get_obs_trak_stats(obsid)

    if er:
        run_obspar = None
        vv = None
        logger.info("Processing ER; no V&V available")
    else:
        # obspar ingested
        try:
            run_obspar = obspar.get_obspar(obsid)
        except:
            run_obspar = None

        # v&v available
        try:
            vv_obi = get_arch_vv(obsid, version='last')
            vv = dict(slots=vv_obi.slot_report)['slots']
            #vv = get_vv(obsid, version='last')['slots']
            #vv_files = get_vv_files(obsid, version='last')
        except LookupError:
            logger.info("No V&V available")
            vv = None

    cat_table = catalog_info(obs_sc['catalog'], acqs, trak, vv)

    for row in cat_table:
        if row['type'] != 'FID':
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
            row['idlink'] = '<A HREF="star_{id}.html">{id}</A>'.format(id=int(row['id']))
        else:
            row['idlink'] = int(row['id'])


    template = jinja_env.get_template('report.html')
    page = template.render(cat_table=cat_table,
                           obs=obs,
                           links=links,
                           target=summary,
                           obsid=obsid)
    full_report_file = os.path.join(outdir, 'report.html')
    logger.info("Writing out full report to {}".format(full_report_file))
    f = open(full_report_file, 'w')
    f.write(page)
    f.close()
    raise ValueError

if __name__ == '__main__':
    main()
