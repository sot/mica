import tables

vv_desc = dict(
obsid=tables.IntCol(pos=0),
revision=tables.IntCol(pos=1),
most_recent=tables.IntCol(pos=2),
tstart=tables.FloatCol(pos=3),
tstop=tables.FloatCol(pos=4),
sim_z=tables.FloatCol(pos=5),
sim_z_offset=tables.FloatCol(pos=6),
instrument=tables.StringCol(10,pos=7),
ra_pnt=tables.FloatCol(pos=8),
dec_pnt=tables.FloatCol(pos=9),
roll_pnt=tables.FloatCol(pos=10),
slot=tables.IntCol(pos=11),
type=tables.StringCol(10,pos=12),
n_pts=tables.IntCol(pos=13),
rad_off=tables.FloatCol(pos=14),
frac_dy_big=tables.FloatCol(pos=15),
frac_dz_big=tables.FloatCol(pos=16),
frac_mag_big=tables.FloatCol(pos=17),
mean_y =tables.FloatCol(pos=18),
mean_z =tables.FloatCol(pos=19),
dy_mean=tables.FloatCol(pos=20),
dy_med =tables.FloatCol(pos=21),
dy_rms =tables.FloatCol(pos=22),
dz_mean=tables.FloatCol(pos=23),
dz_med =tables.FloatCol(pos=24),
dz_rms =tables.FloatCol(pos=25),
dr_mean=tables.FloatCol(pos=26),
dr_med =tables.FloatCol(pos=27),
dr_rms =tables.FloatCol(pos=28),
mag_mean=tables.FloatCol(pos=29),
mag_med =tables.FloatCol(pos=30),
mag_rms =tables.FloatCol(pos=31),
mean_aacccdpt=tables.FloatCol(pos=32),
)

h5f = tables.openFile('vv.h5', 'a')
tbl = h5f.createTable('/', 'vv', vv_desc)
tbl.cols.obsid.createIndex()
h5f.close()
