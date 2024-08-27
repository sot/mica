# Licensed under a 3-clause BSD style license - see LICENSE.rst
import tables
import numpy as np

VV_DTYPE = np.dtype(
    [
        ('obsid', '<i4'),
        ('revision', '<i4'),
        ('isdefault', '<i4'),
        ('aspect_1_id', '<i4'),
        ('used', '<i4'),
        ('vv_version', '<i4'),
        ('ap_date', '|S21'),
        ('tstart', '<f8'),
        ('tstop', '<f8'),
        ('sim_z', '<f8'),
        ('sim_z_offset', '<f8'),
        ('instrument', '|S10'),
        ('ra_pnt', '<f8'),
        ('dec_pnt', '<f8'),
        ('roll_pnt', '<f8'),
        ('slot', '<i4'),
        ('type', '|S10'),
        ('n_pts', '<i4'),
        ('rad_off', '<f8'),
        ('frac_dy_big', '<f8'),
        ('frac_dz_big', '<f8'),
        ('frac_mag_big', '<f8'),
        ('mean_y', '<f8'),
        ('mean_z', '<f8'),
        ('dy_mean', '<f8'),
        ('dy_med', '<f8'),
        ('dy_rms', '<f8'),
        ('dz_mean', '<f8'),
        ('dz_med', '<f8'),
        ('dz_rms', '<f8'),
        ('dr_mean', '<f8'),
        ('dr_med', '<f8'),
        ('dr_rms', '<f8'),
        ('mag_mean', '<f8'),
        ('mag_med', '<f8'),
        ('mag_rms', '<f8'),
        ('mean_aacccdpt', '<f8'),
    ]
)

vv_desc, byteorder = tables.descr_from_dtype(VV_DTYPE)
filters = tables.Filters(complevel=5, complib='zlib')
h5f = tables.openFile('vv.h5', 'a')
tbl = h5f.createTable('/', 'vv', vv_desc, filters=filters, expectedrows=1e6)
h5f.close()
