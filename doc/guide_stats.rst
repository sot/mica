Guide star statistics
---------------------

The :mod:`mica.stats.guide_stats` module
includes code to gather acquisition and guide statistics data for each observation and
return those data to the user.

To get the whole guide stats data table:

   >>> from mica.stats.guide_stats import get_stats
   >>> stats = get_stats()
   >>> stats[(stats['obsid'] == 5438) & (stats['slot'] == 3)][0]['dy_std']
   0.18406899294435713

The hdf5 in-kernel searches may be faster working with the table directly for some
operations.

Data table fields
^^^^^^^^^^^^^^^^^^
The guide statistics data table contains the following columns:

======================= ====================================================================
 Column                 Description
======================= ====================================================================
obsid                   obsid
obi                     observation interval number
kalman_tstart           Kalman transition tstart from kadi event
kalman_datestart        Kalman transition datestart from kadi event
kalman_datestop         Kalman transition datestop from kadi event
revision                revision string representing the software version
slot                    ACA readout slot id
idx                     starcheck star catalog index id
type                    star catalog type (BOT or GUI)
yang                    commanded y-angle position for center of readout window (arcsec)
zang                    commanded z-angle position for center of readout window (arcsec)
sz                      image readout window size
mag                     catalog MAG_ACA of acquisition star
n_samples               number of telemetry samples in Kalman (only samples with good telem quality)
n_track                 number of samples with AOACFCT == 'TRAK'
f_track                 fraction of samples with AOACFCT == 'TRAK'
f_racq                  fraction samples AOACFCT == 'RACQ'
f_srch                  fraction samples AOACFCT == 'SRCH'
f_none                  fraction samples AOACFCT == 'NONE'
n_kalman                number of tracked samples with no (modern) bad status flags set (no IR, SP)
no_track                fraction of no track samples (n_samples - n_track / n_samples)
f_within_0.3            fraction of tracked samples with 0.3 arcsecs of expected position
f_within_1              fraction of tracked samples with 1 arcsecs of expected position
f_within_5              fraction of tracked samples with 5 arcsecs of expected position
f_outside_5             fraction of tracked samples outside 5 arcsecs of expected position
f_obc_bad               fraction of tracked samples with modern bad status flags set (IR, SP)
f_common_col            fraction of tracked samples with AOACICC flag set to ERR
f_quad_bound            fraction of tracked samples with AOACIQB flag set to ERR
f_sat_pix               fraction of tracked samples with AOACISP flag set to ERR
f_def_pix               fraction of tracked samples with AOACIDP flag set to ERR
f_ion_rad               fraction of tracked samples with AOACIIR flag set to ERR
f_mult_star             fraction of tracked samples with AOACIMS flag set to ERR
aoacmag_min             min magnitude over interval
aoacmag_mean            mean magnitude over interval
aoacmag_max             max magnitude over interval
aoacmag_std             std of magnitude over interval
aoacyan_mean            mean AOACYAN
aoaczan_mean            mean AOACZAN
dy_min                  dy min (see notes on residuals)
dy_mean                 dy max (see notes on residuals)
dy_std                  dy std (see notes on residuals)
dy_max                  dy max (see notes on residuals)
dz_min                  dz min (see notes on residuals)
dz_mean                 dz max (see notes on residuals)
dz_std                  dz std (see notes on residuals)
dz_max                  dz max (see notes on residuals)
dr_min                  dr min (see notes on residuals)
dr_mean                 dr max (see notes on residuals)
dr_std                  dr std (see notes on residuals)
dr_5th                  dr 5th percentile value
dr_95th                 dr 95th percentile value
dr_max                  dr max (see notes on residuals)
n_track_interv          number of intervals of continuous track
n_long_track_interv     number of intervals of continuous track with > 60 samples
n_long_no_track_interv  number of intervals of continuous loss of track with > 60 samples
n_racq_interv           number of intervals during Kalman when RACQ was set for the slot
n_srch_interv           number of intervals during Kalman when SRCH was set
agasc_id                AGASC (catalog) id
color1                  AGASC COLOR1, estimated B-V color
ra                      AGASC right ascension (degrees)
dec                     AGASC declination (degrees)
epoch                   AGASC EPOCH
pm_ra                   AGASC PM_RA, proper motion in ra (milli-arcsec/year)
pm_dec                  AGASC PM_DEC, proper motion in dec (milli-arsec/year)
var                     AGASC VAR, known or suspected variable star
pos_err                 AGASC POS_ERR, position error (milli-arcsec)
mag_aca                 AGASC MAG_ACA, ACA mag
mag_err                 AGASC MAG_ERR, (to be fixed, this should have been MAG_ACA_ERR)
mag_band                AGASC MAG_BAND, integer code for spectral band
pos_catid               AGASC POS_CATID, integer code for position catalog
aspq1                   AGASC ASPQ1, integer spoiler code
aspq2                   AGASC ASPQ2, integer proper motion flag
aspq3                   AGASC ASPQ3, integer dist in 100m-arcsec to nearest Tycho2 star
acqq1                   AGASC ACQQ1, mag diff to brightest star within 53.3" (0.01mags)
acqq2                   AGASC ACQQ2, mag diff to brightest star within 107" (0.01mags)
acqq4                   AGASC ACQQ4, mag diff to brightest star within 267.5" (0.01mags)
n100_warm_frac          estimated n100 fraction of CCD pixels for this observation
tccd_mean               mean CCD temperature over Kalman interval
tccd_max                max CCD temperature over Kalman interval
known_bad               ignore this star in standard processing (boolean)
bad_comment             reason to ignore a "known_bad" star
======================= ====================================================================

Notes on residuals:

* The dy, dz, dr values are residuals calculated by subtracting the
  expected star position from the telemetered star position (AOACYAN AOACZAN).
* The expected star position has been calculated using the onboard estimated attitude (AOATTQT*) and the
  AGASC RA/Dec for the commanded star.
* ``dr`` is defined as ``sqrt(dy**2 + dz**2)``

Processing
^^^^^^^^^^

For each observation, after the observation has run and telemetry is available, the guide
stats process does the following:

* Fetches the AGASC information for each star in the catalog
* Fetches the PCAD data for the Kalman interval
* For each guide star in the Kalman interval star calculates statistics on metrics for
  that star over the interval.

