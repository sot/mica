CREATE TABLE archfiles (
  filename        text not null,
  obsid int,        

title   text,
observer        text,
ao      text,
object  text,
ss_object       text,
obs_id  text,
obi_num int,
seq_num text,
instrume        text,
grating text,
detector        text,
detnam  text,
si_mode text,
optical_monitor text,
raster_scan     text,
dither  text,
dither_y_amp float, 
dither_y_freq float,
dither_y_phase float,
dither_z_amp float,
dither_z_freq float,
dither_z_phase float,
ra_pnt  float,
dec_pnt float,
roll_pnt        float,
ra_targ float,
dec_targ        float,
y_det_offset    float,
z_det_offset    float,
radial_offset   float,
defocus float,
sim_z_offset    float,
pre_id  text,
uninterrupted   text,
seg_max_num     int,
ra_nom  float,
dec_nom float,
roll_nom        float,
date_obs        text,
date_end        text,
tstart  float,
tstop   float,
sched_start     float,
sched_stop      float,
sched_exp_time  float,
obs_mode        text,
maneuver        text,
maneuver_v1     float,
maneuver_v2     float,
maneuver_v3     float,
maneuver_angle  float,
maneuver_ref    text,
mjdref  float,
timezero        float,
timeunit        text,
timesys text,
timversn        text,
datamode        text,
readmode        text,
ccdi0_on        text,
ccdi1_on        text,
ccdi2_on        text,
ccdi3_on        text,
ccds0_on        text,
ccds1_on        text,
ccds2_on        text,
ccds3_on        text,
ccds4_on        text,
ccds5_on        text,
dropped_chip_count      int,
onchip_sum      text,
sumrow  int,
sumcol  int,
subarray        text,
startrow        int,
rowcnt  int,
subarray_frame_time     float,
duty_cycle      text,
dtycycle        int,
exptimea        float,
exptimeb        float,
most_efficient  text,
eventfilter     text,
phamin  float,
pharange        float,
bias_request    text,
mode    text,
mission text,
telescop        text,
sim_x   float,
sim_y   float,
sim_z   float,
foc_len float,
py_shutter_position int,
range_switch_level int,
antico_enable text,
upper_level_disc int,
timing_mode text,
trigger_level int,
u_blank_hi int,
u_blank_low int,
v_blank_low int,
v_blank_hi int,
uld_enable tetxt,
zero_block text,
blank_enable text,
width_enable text,
spect_mode text,
my_shutter_position int,
width_threshold int,
mt_a float,
mt_aop float,
mt_e float,
mt_epoch float,
mt_i float,
mt_ma float,
mt_raan float,

timeref text,
tassign text,
origin  text,
ascdsver        text,
obi_data        text,
revision        int,
obspar_ver      int,
obspar_type     text,
obspar_stat     text,
content         text default 'OBSPAR',
isdefault       int,

  CONSTRAINT pk_archfiles PRIMARY KEY (filename)
);

CREATE INDEX idx_archfiles_obsid ON archfiles (obsid);
CREATE INDEX idx_archfiles_tstart ON archfiles (tstart);
