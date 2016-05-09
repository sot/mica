from astropy.table import Table
from mica.stats import acq_stats
from Ska.DBI import DBI


def get_acq_data(agasc_id):
    acq = Table(acq_stats.get_stats())
    acq_star = acq[acq['agasc_id'] == agasc_id]
    acq_table = []
    acq_star.sort('guide_start')
    for s in acq_star:
        srec = {}
        for col in ['type', 'obsid', 'obi', 'slot', 'mag', 'mag_obs', 'star_tracked']:
            srec[col] = s[col]
        srec['date'] = s['guide_start']
        srec['acq_dy'] = s['cdy']
        srec['acq_dz'] = s['cdz']
        srec['id'] = s['acqid']
        acq_table.append(srec)
    return acq_table


def get_gui_data(agasc_id):
    db = DBI(dbi='sybase', server='sybase', user='aca_read')
    gui = db.fetchall('select * from trak_stats_data where id = {}'.format(
            agasc_id))
    gui = Table(gui)
    gui.sort('kalman_datestart')
    gui_table = []
    for s in gui:
        srec = {}
        for col in ['type', 'obsid', 'obi', 'slot']:
            srec[col] = s[col]
        srec['date'] = s['kalman_datestart']
        srec['mag'] = s['mag_exp']
        srec['mag_obs'] = s['aoacmag_mean']
        srec['perc_not_track'] = s['not_tracking_samples'] * 100.0 / s['n_samples']
        gui_table.append(srec)
    return gui_table


def get_star_stats(agasc_id):
    acq_table = get_acq_data(agasc_id)
    gui_table = get_gui_data(agasc_id)
    return acq_table, gui_table
