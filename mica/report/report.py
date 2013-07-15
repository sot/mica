import os
import jinja2
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits

import Ska.DBI

from mica.archive import obspar
from mica.catalog import catalog
from mica.starcheck import starcheck
from mica.archive import asp_l1

plt.rcParams['lines.markeredgewidth'] = 0

def get_starcheck(obsid):
    mp_dir, status = starcheck.get_mp_dir(obsid)
    sc = starcheck.obsid(obsid, mp_dir)
    return (sc, mp_dir, status)


def get_acq_stats(obsid):
    acqs = aca_db.fetchall("select * from acq_stats_data where obsid = %d" % obsid)
    if len(acqs):
        return Table(acqs)

def get_trak_stats(obsid):
    guis = aca_db.fetchall("select * from trak_stats_data where obsid = %d" % obsid)

    if len(guis):
        return Table(guis)

aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')

#obsid = 14932
obsid = 10980
#obsid = 2306
#obsid = 15315
#obsid = 15293
outdir = 'test_report_out'
if not os.path.exists(outdir):
    os.makedirs(outdir)

try:
    obspar = obspar.get_obspar(obsid)
except:
    obspar = dict()


jinja_env = jinja2.Environment(
    loader=jinja2.FileSystemLoader('templates'))
jinja_env.line_comment_prefix = '##'
jinja_env.line_statement_prefix = '#'

# Starcheck
obs_sc, mp_dir, status = get_starcheck(obsid)
fig, cat, obs = catalog.plot(obsid, mp_dir)
fig.savefig(os.path.join(outdir, 'starcheck.png'))


sc_template = jinja_env.get_template('starcheck.html')
sc_catalog = [dict(zip(cat.dtype.names, row))
              for row in cat]
for row in sc_catalog:
    row['strid'] = (row['idnote'] or '') + str(row['id'])
    row['passnotes'] = (row['pass'] or '') + (row['notes'] or '') or ''
#page = sc_template.render(obs=obs,
#                          catalog=dcat)
#f = open(os.path.join(outdir, 'starcheck.html'), 'w')
#f.write(page)
#f.close()

# acq stats
if status is not None and (status == 'ran' or status == 'ran_pretimelines') :
    acqs = get_acq_stats(obsid)
    acq_patches = []
    for acq_star in acqs[acqs['obc_id'] == 'NOID']:
        acq_box = plt.Rectangle(
            (acq_star['yang'] - (2 * acq_star['halfw']),
             acq_star['zang'] - (2 * acq_star['halfw'])),
            width=acq_star['halfw'] * 4,
            height=acq_star['halfw'] * 4,
            color='red',
            alpha='.5',
            fill=True)
        acq_patches.append(fig.gca().add_patch(acq_box))
    fig.savefig(os.path.join(outdir, 'acqstat_catalog.png'))
    [p.set_visible(False) for p in acq_patches]
    sc_template = jinja_env.get_template('acq_stats.html')
    #dcat = [dict(zip(cat.dtype.names, row))
    #       for row in cat]
    acq_catalog = acqs
    #page = sc_template.render(obs=obs,
    #                          catalog=acqs)
    #f = open(os.path.join(outdir, 'acq_stats.html'), 'w')
    #f.write(page)
    #f.close()

import mica.vv
try:
    vv, files = mica.vv.get_vv(obsid)
except ValueError:
    vv, files = mica.vv.get_vv(obsid, version='last')


# stars used in kalman filter and trak stats
trak_stats = get_trak_stats(obsid)
gspr = asp_l1.get_files(obsid, content=['GSPROPS'])
first_gspr = fits.open(gspr[0])[1].data
fidpr = asp_l1.get_files(obsid, content=['FIDPROPS'])
first_fidpr = fits.open(fidpr[0])[1].data

trak_patches = []
dtrak_stats = [dict(zip(trak_stats.dtype.names, row))
               for row in trak_stats]
fid_stats = []
star_stats = []
for tracked in dtrak_stats:
    for stat in ['def_pix', 'ion_rad', 'mult_star', 'quad_bound', 'common_col', 'sat_pix',
                 'obc_bad_status', 'bad_status', 'not_tracking']:
        stat_key = '{}_frac'.format(stat)
        if tracked['n_samples'] == 0:
            tracked[stat_key] = 0
        else:
            tracked[stat_key] = (tracked["{}_samples".format(stat)] * 1.0 
                                 / tracked['n_samples'])
    if tracked['type'] == 'FID':
        fidpr_rec = first_fidpr[first_fidpr['slot'] == tracked['slot']]
        if not len(fidpr_rec) or fidpr_rec['id_status'][0] != 'GOOD':
            tracked['id_status'] = 'ERR'
            fid_circle = plt.Circle(
                (tracked['cyan_exp'],
                 tracked['czan_exp']),
                radius=250,
                color='red',
                alpha='.5',
                fill=True)
            trak_patches.append(fig.gca().add_patch(fid_circle))
        else:
            tracked['id_status'] = fidpr_rec['id_status'][0]
        tracked['cel_loc_flag'] = 0
        if str(tracked['slot']) in vv['slots']:
            vv['slots'][str(tracked['slot'])].update(dict(cel_loc_flag=0))
            vv['slots'][str(tracked['slot'])].update(dict(id_status=tracked['id_status']))
        fid_stats.append(tracked)
        continue
    gspr_rec = first_gspr[first_gspr['slot'] == tracked['slot']]
    tracked['cel_loc_flag'] = 0 # this should be the right answer if there is no line for the slot
    tracked['id_status'] = 'ERR'
    if len(gspr_rec):
        tracked['cel_loc_flag'] = gspr_rec['cel_loc_flag'][0]
        tracked['id_status'] = gspr_rec['id_status'][0]
    if str(tracked['slot']) in vv['slots']:
        vv['slots'][str(tracked['slot'])].update(dict(cel_loc_flag=tracked['cel_loc_flag']))
        vv['slots'][str(tracked['slot'])].update(dict(id_status=tracked['id_status']))
    if not len(gspr_rec) or gspr_rec['cel_loc_flag'][0] == 0:
        gui_circle = plt.Circle(
            (tracked['cyan_exp'],
             tracked['czan_exp']),
            radius=250,
            color='red',
            alpha='.5',
            fill=True)
        trak_patches.append(fig.gca().add_patch(gui_circle))
    star_stats.append(tracked)
fig.savefig(os.path.join(outdir, 'trakinfo_catalog.png'))
[p.set_visible(False) for p in trak_patches]

#sc_template = jinja_env.get_template('trak_stats.html')
#page = sc_template.render(obs=obs,
#                          fid_stats=fid_stats,
#                          star_stats=star_stats,
#                          )
#f = open(os.path.join(outdir, 'trak_stats.html'), 'w')
#f.write(page)
#f.close()



#fidpr = asp_l1.get_files(obsid, content=['FIDPROPS'])
#first_fidpr = fits.open(fidpr[0])[1].data
#fid_patches = []
#for tracked in trak_stats:
#    if tracked['type'] != 'FID':
#        continue
#    fidpr_rec = first_fidpr[first_fidpr['slot'] == tracked['slot']]
#    if not len(fidpr_rec) or fidpr_rec['id_status'] != 'GOOD':
#        fid_circle = plt.Circle(
#            (tracked['cyan_exp'],
#             tracked['czan_exp']),
#            radius=250,
#            color='red',
#            alpha='.5',
#            fill=True)
#        fid_patches.append(fig.gca().add_patch(fid_circle))
#
#
#

big_template = jinja_env.get_template('one_big.html')
page = big_template.render(obs=obs,
                           obspar=obspar,
                       sc_catalog=sc_catalog,
                       acq_catalog=acq_catalog,
                       star_stats=star_stats,
                       fid_stats=fid_stats,
                       slots=sorted(vv['slots']),
                       vv_slots=vv['slots'])
f = open(os.path.join(outdir, 'obsid.html'), 'w')
f.write(page)
f.close()

#template = jinja_env.get_template('vv_slots.html')
#page = template.render(slots=sorted(vv['slots']),
#                       vv_slots=vv['slots'] )
#f = open(os.path.join(outdir, 'vv_slots.html'), 'w')
#f.write(page)
#f.close()

for slot in vv['slots']:
    template = jinja_env.get_template('vv_slots_single.html')
    page = template.render(obsdir=os.path.dirname(files[0]),
                           slot=slot)
    f = open(os.path.join(outdir, 'slot_{}.html'.format(slot)), 'w')
    f.write(page)
    f.close()
    

#
#

#kal.plot()



#
#asols =  asp_l1.get_files(obsid, content=['ASPSOL'])

#trak_data = get_trak_data(tstart, tstop)











#asp_db = Ska.DBI.DBI(dbi='sqlite',
#                     server='/data/aca/archive/asp1/processing_asp_l1.db3')



#asol = asp_l1.get_files(obsid, content=['ASPSOL'])


#aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')
#
## put in mica.starcheck
#sc_db = Ska.DBI.DBI(dbi='sqlite',
#                           server='/data/aca/archive/starcheck/starcheck.db3')
#in_starcheck = sc_db.fetchall("select * from starcheck_obs where obsid = %d"
#                              % obsid)
#in_cmd_states = aca_db.fetchall("select * from cmd_states "
#                                + "where obsid = %d " % obsid
#                                + "and pcad_mode = 'NPNT'")
#

