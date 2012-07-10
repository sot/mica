import tables

vv_desc = dict(
obsid=tables.IntCol(pos=0),
revision=tables.IntCol(pos=1),
most_recent=tables.IntCol(pos=2),
slot=tables.IntCol(pos=3),
type=tables.StringCol(10,pos=4),
n_pts=tables.IntCol(pos=5),
rad_off=tables.FloatCol(pos=6),
frac_dy_big=tables.FloatCol(pos=7),
frac_dz_big=tables.FloatCol(pos=8),
frac_mag_big=tables.FloatCol(pos=9),
mean_y =tables.FloatCol(pos=10),
mean_z =tables.FloatCol(pos=11),
dy_mean=tables.FloatCol(pos=12),
dy_med =tables.FloatCol(pos=13),
dy_rms =tables.FloatCol(pos=14),
dz_mean=tables.FloatCol(pos=15),
dz_med =tables.FloatCol(pos=16),
dz_rms =tables.FloatCol(pos=17),
mag_mean=tables.FloatCol(pos=18),
mag_med =tables.FloatCol(pos=19),
mag_rms =tables.FloatCol(pos=20),
)

h5f = tables.openFile('/data/aca/archive/vv/vv.h5', 'a')
tbl = h5f.createTable('/', 'vv', vv_desc)
tbl.cols.obsid.createIndex()
h5f.close()
