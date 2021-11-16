# Licensed under a 3-clause BSD style license - see LICENSE.rst
from functools import partial

from astropy.table import Table
import numpy as np
import pytest

from mica.archive.cda import (get_archive_file_list, get_ocat_details_local,
                              get_ocat_details_web, get_ocat_summary_web,
                              get_proposal_abstract)
from mica.archive.cda.services import update_ocat_local

from pathlib import Path


DETAIL_COLNAMES = [
    'seq_num', 'status', 'obsid', 'pr_num', 'target_name', 'grid_name', 'instr', 'grat', 'type',
    'obs_cycle', 'prop_cycle', 'charge_cycle', 'start_date', 'public_avail', 'readout_detector',
    'datamode', 'joint', 'hst', 'noao', 'nrao', 'rxte', 'spitzer', 'suzaku', 'xmm', 'swift',
    'nustar', 'category', 'seg_max_num', 'prop_title', 'pi_name', 'observer', 'app_exp',
    'exp_time', 'ra', 'dec', 'soe_roll', 'time_crit', 'y_off', 'z_off', 'x_sim', 'z_sim',
    'raster', 'obj_type', 'obj', 'photo', 'vmag', 'est_cnt_rate', 'forder_cnt_rate',
    'count_rate', 'event_count', 'dither', 'y_amp', 'y_freq', 'y_phase', 'z_amp', 'z_freq',
    'z_phase', 'roll', 'window', 'unint', 'pointing_update', 'monitor', 'pre_id', 'mon_min',
    'mon_max', 'group_id', 'constr', 'epoch', 'period', 'pstart', 'ps_marg', 'pend', 'pe_marg',
    'multitel', 'multitel_obs', 'multitel_int', 'constr_rmk', 'too_type', 'too_start',
    'too_stop', 'alt_group', 'alt_trig', 'simode', 'hrc', 'spect_mode', 'blank_en', 'u_hi',
    'v_hi', 'u_lo', 'v_lo', 'timing', 'z_blk', 'acis', 'mode', 'bep_pack', 'dropped_chip_cnt',
    'i0', 'i1', 'i2', 'i3', 's0', 's1', 's2', 's3', 's4', 's5', 'spectra_max_count',
    'multiple_spectral_lines', 'subary', 'strt_row', 'row_cnt', 'd_cyc', 'sec_cnt',
    'pr_time', 'sec_time', 'f_time', 'oc_sum', 'oc_row', 'oc_col', 'evfil', 'evfil_lo',
    'evfil_ra', 'efficient', 'spwin'
]

# Ensure tests continue to work in future by always using a fixed date range
DATE_RANGE = '1999-01-01/2021-11-01'
get_ocat_details_web = partial(get_ocat_details_web, startDate=DATE_RANGE)
get_ocat_summary_web = partial(get_ocat_summary_web, startDate=DATE_RANGE)


@pytest.fixture(scope="session")
def datafile(tmp_path_factory):
    """Make a temp HDF5 Ocat details file within 60 arcmin of 3c273 for obsids
    before 2021-Nov that persists for the testing session."""
    datafile = str(tmp_path_factory.mktemp('ocat') / 'target_table.h5')
    update_ocat_local(datafile, target_name='3c273', resolve_name=True, radius=60,
                      startDate=DATE_RANGE)
    return datafile


def test_ocat_details_local_all(datafile):
    dat = get_ocat_details_local(datafile=datafile)
    assert dat.colnames == DETAIL_COLNAMES
    assert len(dat) == 38
    ok = dat['obsid'] == 7781
    assert np.all(dat['target_name'][ok] == "SDSS J123132.37+013814.1")
    assert set(dat['target_name']) == {
        ' ',
        '3C 273',
        '3C273',
        '3C273-JET',
        'SDSS J123132.37+013814.1',
        'SDSSJ123215.81+020610.0'}


def test_ocat_details_local_where(datafile):
    dat = get_ocat_details_local(datafile=datafile, target_name='jet')
    assert len(dat) == 4
    assert all('JET' in row['target_name'] for row in dat)


def test_ocat_details_local_filtering(datafile):
    dat = get_ocat_details_local(datafile=datafile, target_name='jet',
                                 public_avail='2004-12-05 10:16:12',
                                 where='strt_row=115',
                                 obsid=4876)
    assert isinstance(dat, dict)
    assert dat['obsid'] == 4876


def test_ocat_details_local_filtering_no_match1(datafile):
    dat = get_ocat_details_local(datafile=datafile, target_name='jet',
                                 public_avail='2004-12-05 10:16:12',
                                 where='strt_row==115000',  # no match
                                 obsid=4876, return_type='table')
    assert len(dat) == 0


def test_ocat_details_local_filtering_no_match2(datafile):
    dat = get_ocat_details_local(datafile=datafile, target_name='jet',
                                 public_avail='2004-12-05 10:16:13',  # no match
                                 where='strt_row==115',
                                 obsid=4876, return_type='table')
    assert len(dat) == 0


def test_ocat_details_local_filtering_no_match3(datafile):
    dat = get_ocat_details_local(datafile=datafile, target_name='jetttttt',  # no match
                                 public_avail='2004-12-05 10:16:12',
                                 where='strt_row==115',
                                 obsid=4876, return_type='table')
    assert len(dat) == 0


def test_ocat_details_local_filtering_no_match4(datafile):
    dat = get_ocat_details_local(datafile=datafile, target_name='jet',
                                 public_avail='2004-12-05 10:16:12',
                                 where='strt_row==115',
                                 obsid=4876000, return_type='table')  # no match
    assert len(dat) == 0


def test_ocat_details_local_obsid_auto(datafile):
    """Test obsid with 'auto' return"""
    dat = get_ocat_details_local(datafile=datafile, obsid=4911)
    assert isinstance(dat, dict)
    assert dat['target_name'] == 'SDSSJ123215.81+020610.0'
    assert dat['pre_id'] is np.ma.masked


def test_ocat_details_local_obsid_table(datafile):
    dat = get_ocat_details_local(datafile=datafile, obsid=4911, return_type='table')
    assert isinstance(dat, Table)
    row = dat[0]
    assert row['target_name'] == 'SDSSJ123215.81+020610.0'
    assert row['pre_id'] is np.ma.masked


def test_ocat_details_local_no_match_auto(datafile):
    with pytest.raises(ValueError, match=r'failed to find obsid 999999'):
        get_ocat_details_local(datafile=datafile, obsid=999999)


def test_ocat_details_local_no_match_table(datafile):
    dat = get_ocat_details_local(datafile=datafile, obsid=999999, return_type='table')
    assert len(dat) == 0
    assert dat.colnames == DETAIL_COLNAMES


def test_get_archive_file_list():
    dat = get_archive_file_list(obsid=2365, detector='pcad',
                                subdetector='aca', level=1, filetype='aspsol')
    assert len(dat) == 2
    assert dat['Filename'][0].startswith('pcadf1200418')
    assert dat.colnames == ['Filename', 'Filesize', 'Timestamp']


def test_get_proposal_abstract():
    exp = {
        'abstract': ('We propose the Chandra-COSMOS survey which will provide an '
                     'unprecedented combination of contiguous area, depth and '
                     'resolution. 36 densely tiled observations will cover the central '
                     '0.7 sq.deg. COSMOS field to a uniform 200ksec depth. COSMOS '
                     'explores the coupled evolution of galaxies, dark matter halos '
                     'and AGNs (massive black holes) largely free of cosmic variance. '
                     'COSMOS is a comprehensive survey including: HST, Spitzer, '
                     'Subaru, VLT, Magellan, VLA, MAMBO, GALEX, & potentially EVLA & '
                     'ALMA. Chandra resolution & sensitivity enables the study of '
                     'large scale phenomena: (1) influence of the surrounding '
                     'environment; (2) interaction between galaxies; (3) influence of '
                     'groups and clusters'),
        'principal_investigator': 'Martin Elvis',
        'proposal_number': '08900073',
        'proposal_title': 'THE CHANDRA-COSMOS SURVEY'
    }

    dat = get_proposal_abstract(obsid=8000)
    assert dat == exp

    dat = get_proposal_abstract(propnum='08900073')
    assert dat == exp


def test_get_proposal_abstract_fail():
    with pytest.raises(ValueError, match=r'must provide obsid or propnum'):
        get_proposal_abstract()
