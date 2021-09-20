# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Update a mica Ska target table
"""
from pathlib import Path
import argparse
from astropy.table import Table, join

from cxotime import CxoTime
import Ska.DBI


def get_options():
    parser = argparse.ArgumentParser(
        description="Update target table")
    parser.add_argument("--datafile",
                        default='target_table.ecsv')
    opt = parser.parse_args()
    return opt


def main():
    opt = get_options()
    table = get_table()
    table.write(opt.datafile, overwrite=True)


def get_table():
    """
    Fetch target, acisparam, hrcparam, and dither Sybase tables to make a target table
    data product.
    """

    # These all look to have 1=1 mapping on Chaser output or map to
    # another public table (acisparam dither hrcparam etc)
    targ_cols = [('obsid', 'i4'),
                 ('targid', 'f4'),
                 ('seq_nbr', 'U6'),
                 ('targname', 'U50'),
                 ('obj_flag', 'U4'),
                 ('object', 'U20'),
                 ('si_mode', 'U10'),
                 ('photometry_flag', 'U4'),
                 ('vmagnitude', 'f8'),
                 ('ra', 'f8'),
                 ('dec', 'f8'),
                 ('est_cnt_rate', 'f4'),
                 ('forder_cnt_rate', 'f4'),
                 ('y_det_offset', 'f4'),
                 ('z_det_offset', 'f4'),
                 ('raster_scan', 'U4'),
                 ('dither_flag', 'U4'),
                 ('uninterrupt', 'U4'),
                 ('approved_exposure_time', 'f8'),
                 ('pre_min_lead', 'f8'),
                 ('pre_max_lead', 'f8'),
                 ('pre_id', 'f4'),
                 ('phase_constraint_flag', 'U4'),
                 ('acisid', 'f4'),
                 ('hrcid', 'f4'),
                 ('grating', 'U5'),
                 ('instrument', 'U7'),
                 ('type', 'U30'),
                 ('lts_lt_plan', 'U24'),
                 ('status', 'U20'),
                 ('observatories', 'U60'),
                 ('constr_in_remarks', 'U4'),
                 ('obs_ao_str', 'U10'),
                 ('window_flag', 'U4'),
                 ('roll_flag', 'U4'),
                 ('public_avail', 'U24'),
                 ('spwindow_flag', 'U4'),
                 ('charge_ao_str', 'U10'),
                 ('pointing_constraint', 'U4'),
                 ('group_id', 'U50'),
                 ('multitelescope_interval', 'f8'),
                 ('remarks', 'U200'),
                 ('mp_remarks', 'U200')]

    acis_cols = [('acisid', 'i4'),
                 ('exp_mode', 'U2'),
                 ('ccdi0_on', 'U2'),
                 ('ccdi1_on', 'U2'),
                 ('ccdi2_on', 'U2'),
                 ('ccdi3_on', 'U2'),
                 ('ccds0_on', 'U2'),
                 ('ccds1_on', 'U2'),
                 ('ccds2_on', 'U2'),
                 ('ccds3_on', 'U2'),
                 ('ccds4_on', 'U2'),
                 ('ccds5_on', 'U2'),
                 ('bep_pack', 'U3'),
                 ('onchip_sum', 'U4'),
                 ('onchip_row_count', 'f8'),
                 ('onchip_column_count', 'f8'),
                 ('frame_time', 'f8'),
                 ('subarray', 'U8'),
                 ('subarray_start_row', 'f8'),
                 ('subarray_row_count', 'f8'),
                 ('duty_cycle', 'U4'),
                 ('secondary_exp_count', 'f8'),
                 ('primary_exp_time', 'f8'),
                 ('secondary_exp_time', 'f8'),
                 ('eventfilter', 'U4'),
                 ('eventfilter_lower', 'f8'),
                 ('eventfilter_higher', 'f8'),
                 ('most_efficient', 'U4'),
                 ('dropped_chip_count', 'f4'),
                 ('multiple_spectral_lines', 'U4'),
                 ('spectra_max_count', 'f8')]

    acis_dither = {'y_amp': 0.002222, 'y_freq': 0.360000, 'y_phase': 0.000000,
                   'z_amp': 0.002222, 'z_freq': 0.509100, 'z_phase': 0.000000}
    hrc_dither = {'y_amp': 0.005556, 'y_freq': 0.331200, 'y_phase': 0.000000,
                  'z_amp': 0.005556, 'z_freq': 0.468400, 'z_phase': 0.000000}
    dither_cols = ['y_amp', 'y_freq', 'y_phase', 'z_amp', 'z_freq', 'z_phase']

    select_cols = [col[0] for col in targ_cols]
    # Remove remarks as these are LOB types
    select_cols.remove('mp_remarks')
    select_cols.remove('remarks')

    # Get the target table but just get remarks as sybase varchars
    targ_query = (f"select {','.join(select_cols)}, "
                  "convert(varchar(200), mp_remarks) as mp_remarks, "
                  "convert(varchar(200), remarks) as remarks from target")

    with Ska.DBI.DBI(dbi='sybase', server='sqlsao', user='aca_ops', database='axafocat',
                     numpy=False) as db:
        targets = db.fetchall(targ_query)

    # Convert the datetime columns to strings on the Python side
    for t in targets:
        for col in ['public_avail', 'lts_lt_plan']:
            try:
                t[col] = CxoTime(t[col]).iso
            except Exception:
                pass

    tables = {}
    # Fetch dither and hrcparams
    with Ska.DBI.DBI(dbi='sybase', server='sqlsao', user='aca_ops', database='axafocat',
                     numpy=False) as db:
        for name in ['dither', 'hrcparam']:
            tables[name] = Table(db.fetchall(f"select * from {name}"))

    # Stick in default dithers.  Assume all but dither_flag = 'N' get dither.
    # This could also be done on the sybase side.
    for t in targets:
        if t['obsid'] in tables['dither']['obsid']:
            ok = t['obsid'] == tables['dither']['obsid']
            t.update(zip(dither_cols, tables['dither'][dither_cols][ok][0]))
        else:
            if t['dither_flag'] != 'N':
                if t['instrument'] in ['HRC-S', 'HRC-I']:
                    t.update(hrc_dither)
                else:
                    t.update(acis_dither)

    # Make an astropy table with types.
    # Many are floats because they handle nans from None sybase types.
    names = [col[0] for col in targ_cols] + dither_cols
    dtypes = [col[1] for col in targ_cols] + ['f4', 'f4', 'f4', 'f4', 'f4', 'f4']
    target_table = Table(targets, names=names, dtype=dtypes)

    # The acisparam table has None-types as well so fetch and make a defined-type table of it
    with Ska.DBI.DBI(dbi='sybase', server='sqlsao', user='aca_ops', database='axafocat',
                     numpy=False) as db:
        acisparam = db.fetchall("select * from acisparam")
        tables['acisparam'] = Table(acisparam,
                                    names=[col[0] for col in acis_cols],
                                    dtype=[col[1] for col in acis_cols])

    # Then join that table on the target table
    target_table = join(target_table, tables['acisparam'], keys=['acisid'], join_type='left')

    # Not sure about how to assign si mode conditionally from acis or hrcparam, so do it explicitly
    hrc_cols = tables['hrcparam'].colnames.copy()
    hrc_cols.remove('si_mode')
    target_table = join(target_table, tables['hrcparam'][hrc_cols], keys=['hrcid'],
                        join_type='left')
    for row in tables['hrcparam']:
        ok = target_table['hrcid'] == row['hrcid']
        target_table['si_mode'][ok] = row['si_mode']

    return target_table


if __name__ == "__main__":
    main()
