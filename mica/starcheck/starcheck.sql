create table starcheck_id
(id             int             not null,
dir             varchar(20)     not null,
mtime           float           not null,
primary key (id)
);


create table starcheck_obs
(sc_id int not null,
obsid		int		not null,
obs_idx         int             not null,
point_ra	float,
point_dec	float,
point_roll	float,
target_id	varchar(30),
sci_instr	varchar(10),
sim_z_offset_steps	int,
sim_z_offset_mm float,
grating		varchar(10),
dither_state	varchar(4),
dither_y_amp	float,
dither_y_period	float,
dither_z_amp	float,
dither_z_period float,
mp_starcat_time		varchar(25),
mp_starcat_vcdu_cnt int,
CONSTRAINT fk_so_id FOREIGN KEY (sc_id) REFERENCES starcheck_id (id)
);

create table starcheck_manvr
(sc_id int not null,
obsid  int not null,
obs_idx int not null,
instance int default 0 not null,
duration_sec  float,
angle_deg     float,
slew_err_arcsec   float,
target_Q1       float,
target_Q2       float,
target_Q3       float,
target_Q4       float,
mp_targquat_time	varchar(25),
mp_targquat_vcdu_cnt int,
CONSTRAINT fk_sm_id FOREIGN KEY (sc_id) REFERENCES starcheck_id (id)
);


create table starcheck_catalog
(sc_id int not null,
obsid 		int		not null,
obs_idx         int             not null,
mp_starcat_time		varchar(25),
idx		int		not null,
slot		int 		not null,
id              int,
idnote		varchar(4),
type		varchar(4) 	not null,
sz		varchar(10),
p_acq           float,
minmag		float,
mag		float,
maxmag		float,
yang		int,
zang		int,
dim		int,
res		int,
halfw		int,
pass            varchar(10),
notes		varchar(10),
CONSTRAINT fk_sc_id FOREIGN KEY (sc_id) REFERENCES starcheck_id (id)
);



create table starcheck_warnings
(sc_id int not null,
obsid		int 		not null,
obs_idx         int             not null,
warning		varchar(100)	not null,
idx		int,
warning_type    varchar(15),
CONSTRAINT fk_sw_id FOREIGN KEY (sc_id) REFERENCES starcheck_id (id)
);


create table starcheck_processing
(sc_id int not null,
mp_path	varchar(30) not null,
type		varchar(30) not null,
value		varchar(100),
section_line_num	int,
last_ap_date	datetime,
CONSTRAINT fk_sp_id FOREIGN KEY (sc_id) REFERENCES starcheck_id (id)
);

create table starcheck_pred_temp
(sc_id int not null,
obsid int not null,
pred_ccd_temp float,
CONSTRAINT fk_spt_id FOREIGN KEY (sc_id) REFERENCES starcheck_id (id)
);

CREATE INDEX mtime_idx on starcheck_id(mtime);
CREATE INDEX cat_id_idx on starcheck_catalog(sc_id);
CREATE INDEX cat_time_idx on starcheck_catalog(mp_starcat_time);
CREATE INDEX id_dir_idx on starcheck_id(dir);
CREATE INDEX manvr_id_idx on starcheck_manvr(sc_id);
CREATE INDEX obs_obs_idx on starcheck_obs(obsid);
CREATE INDEX obs_time_idx on starcheck_obs(mp_starcat_time);
