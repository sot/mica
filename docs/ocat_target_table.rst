OCAT Target Table
---------------------

The :mod:`mica.archive.ocat_target_table` module includes code
to fetch a mica version of the content from the OCAT details
aka Chaser aka target table.

This is the type of content that can be seen directly at:

https://cda.harvard.edu/srservices/ocatDetails.do?obsid=2121

Using this module to get the data for all science observations:

   >>> from mica.archive.ocat_target_table import get_ocat_target_table
   >>> dat = get_ocat_target_table()
   >>> dat[dat['obsid'] == 5438][0]['target_name']
   'R Aqr'

Some of these fields may be described in the chaser help at:

https://cda.harvard.edu/chaser/chaserFieldHelp.html


Data table fields
^^^^^^^^^^^^^^^^^

============================= ================================================================
 Column                       Description
============================= ================================================================
seq_num                       Sequence Number
status                        Status (unobserved, archived, untriggered, etc)
obsid                         Obsid
pr_num                        Proposal Number
target_name                   Target name
grid_name                     database id of grid name if grid observation
instr                         Instrument
grat                          Grating (HETG, LETG, NONE)
type                          Observation type (TOO, GO, GTO)
obs_cycle                     Observation cycle
prop_cycle                    Proposal cycle
charge_cycle                  Charge cycle
start_date                    Observation start date
public_avail                  Date publicly available
readout_detector              Detector (which HRC detector or string of actual ACIS ccds)
datamode                      Instrument data mode
joint                         Joint observatories (string)
hst                           HST time (orbits?)
noao                          NOAO time
nrao                          NRAO time
rxte                          RXTE time
spitzer                       SPITZER time
suzaku                        SUZAKU time
xmm                           XMM time
swift                         SWIFT time
nustar                        NUSTAR time
category                      Science category
seg_max_num
prop_title                    Proposal Title
pi_name                       Principal investigator name
observer                      Observer name
app_exp                       Approved exposure time (ks)
exp_time                      Actual exposure time (ks)
ra                            Target Right Ascension
dec                           Target Declination
soe_roll                      Roll from SOE
time_crit                     Time critical flag
y_off                         Y offset
z_off                         Z offset
x_sim                         X SIM
z_sim                         Z SIM
raster                        Raster flag
obj_type                      Object type
obj                           Solar system object name
photo                         Photometry flag
vmag                          V Mag of object
est_cnt_rate                  Estimated count rate
forder_cnt_rate               First order count rage
count_rate
event_count
dither                        Dither flag
y_amp                         Dither Y amplitude if custom dither
y_freq                        Dither Y frequency if custom dither
y_phase                       Dither Y phase if custom dither
z_amp                         Dither Z amplitude if custom dither
z_freq                        Dither Z frequency if custom dither
z_phase                       Dither Z phase if custom dither
roll                          Roll constraint flag
window                        Window constraint flag
unint                         Uninterrupt constraint flag
pointing_update               Pointing update constraint flag
monitor                       Monitor series flag
pre_id                        Obsid of previous observation in monitor series
mon_min                       Min days from pre_id for monitor observation
mon_max                       Max days from pre_id for monitor observation
group_id                      Database group id
constr
epoch
period
pstart
ps_marg
pend
pe_marg
multitel
multitel_obs
multitel_int
constr_rmk                    Constraint in remarks flag
too_type
too_start
too_stop
alt_group
alt_trig
simode                        Science Instrument (SI) mode
hrc
spect_mode
blank_en
u_hi
v_hi
u_lo
v_lo
timing
z_blk
acis
mode                          ACIS mode (CC or TE)
bep_pack                      ACIS BEP PACK (F, G, VF, F+B)
dropped_chip_cnt              Dropped chip count
i0                            ACIS I0 ccd status (Y, N, optional with number, or D if dropped)
i1                            ACIS I1 ccd status (Y, N, optional with number, or D if dropped)
i2                            ACIS I2 ccd status (Y, N, optional with number, or D if dropped)
i3                            ACIS I3 ccd status (Y, N, optional with number, or D if dropped)
s0                            ACIS S0 ccd status (Y, N, optional with number, or D if dropped)
s1                            ACIS S1 ccd status (Y, N, optional with number, or D if dropped)
s2                            ACIS S2 ccd status (Y, N, optional with number, or D if dropped)
s3                            ACIS S3 ccd status (Y, N, optional with number, or D if dropped)
s4                            ACIS S4 ccd status (Y, N, optional with number, or D if dropped)
s5                            ACIS S5 ccd status (Y, N, optional with number, or D if dropped)
spectra_max_count             Spectra Max Count
multiple_spectral_lines       Multiple spectral lines expected (Y, N)
subary                        ACIS subarray (CUSTOM, NONE)
strt_row                      Start row of ACIS subarray
row_cnt                       Number of rows of ACIS subarray
d_cyc
sec_cnt
pr_time
sec_time
f_time
oc_sum
oc_row
oc_col
evfil
evfil_lo
evfil_ra
efficient                     ACIS use most efficient (Y, N)
spwin                         Spatial window (Y, N)
============================= ================================================================

The HDF5 in-kernel searches may be faster working with the table directly for some
operations.
