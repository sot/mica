Acquisition Data
================

Processing
------------------------------------

For each observation, after the observation has run and telemetry is available:

The acquisition and guide stats process

* fetches the AGASC information for each star in the catalog
* fetches the PCAD data at the end of the acquisition interval

and for each acquisition star determines

* if that star was "successfully" acquired
* what the observed magnitude and position of the star were in the last PCAD telemetry
  readout before the guide transition.

Acquisition data products
-------------------------

The acquisition database table may be retrieved with::

  from mica.stats.acq_stats import get_stats
  acq_data = get_stats()

Alternatively, the raw hdf5 may be read directly.  It includes the following columns:


=============== ====================================================================
 Column         Description
=============== ====================================================================
obsid           the obsid
obi             observation interval number
acq_start       acquisition datestart from kadi event
guide_start     guide transition datestart from kadi event
guide_tstart    guide transition start from kadi event in Chandra secs
one_shot_length one shot delta quaternion arc length in arcsecs for this obsid manvr
revision        revision string representing the software version
slot            ACA readout slot id
idx             starcheck star catalog index id
type            star catalog type (BOT or ACQ)
yang            commanded y-angle position for center of readout window (arcsec)
zang            commanded z-angle position for center of readout window (arcsec)
halfw           acquisition search box half width (arcsec)
mag             catalog MAG_ACA of acquisition star
acqid           acquisition success indicator (boolean)
star_tracked    (boolean) something tracked within 5 arcsecs during acq interval
spoiler_tracked (boolean) something tracked outside 5 arcsecs during acq interval
img_func        image status at transition: can be None, star, spoiler, RACQ, SRCH
n_trak_interv   number of intervals during acq when something was tracked
max_trak_cdy    max "corrected" delta-y during acquisition (arcsec)
mean_trak_cdy   mean "corrected" delta-y during acquisition (arcsec)
min_trak_cdy    min "corrected" delta-y during acquisition (arcsec)
max_trak_cdz    max "corrected" delta-z during acquisition (arcsec)
mean_trak_cdy   mean "corrected" delta-z during acquisition (arcsec)
min_trak_cdz    min "corrected" delta-z during acquisition (arcsec)
max_trak_mag    max observed magnitude during acquisition
mean_trak_mag   mean observed magnitude during acquisition
min_trak_mag    min observed magnitude during acquisition
cdy             "corrected" delta-y at guide transition (arcsec)
cdy             "corrected" delta-z at guide transition (arcsec)
dy              delta-y at guide transition (arcsec)
dz              delta-z at guide transition (arcsec)
ion_rad         ionizing radiation flag set at transition (boolean)
def_pixel       defective pixel flag set at transition
mult_star       multiple star flag set at transition
sat_pix         saturated pixel flag set at transition
mag_obs         observed magnitude at transition
yang_obs        observed y-angle at transition (arcsec)
zang_obs        observed z-angle at transition (arcsec)
agasc_id        AGASC (catalog) id
color1          AGASC COLOR1, estimated B-V color
ra              AGASC right ascension in degrees
dec             AGASC declination in degrees
epoch           AGASC EPOCH
pm_ra           AGASC PM_RA, proper motion in ra in milli-arcsec/year
pm_dec          AGASC PM_DEC, proper motion in dec in milli-arsec/year
var             AGASC VAR, known or suspected varliable star
pos_err         AGASC POS_ERR, position error in milli-arcsec
mag_aca         AGASC MAG_ACA, ACA mag
mag_err         AGASC MAG_ERR, (to be fixed, this should have been MAG_ACA_ERR)
mag_band        AGASC MAG_BAND, integer code for spectral band
pos_catid       AGASC POS_CATID, integer code for position catalog
aspq1           AGASC ASPQ1, integer spoiler code
aspq2           AGASC ASPQ2, integer proper motion flag
aspq3           AGASC ASPQ3, integer dist in 100m-arcsec to nearest Tycho2 star
acqq1           AGASC ACQQ1, mag diff to brightest star within 53.3" (unit 0.01mags)
acqq2           AGASC ACQQ2, mag diff to brightest star within 107" (unit 0.01mags)
acqq4           AGASC ACQQ4, mag diff to brightest star within 267.5" (unit 0.01mags)
n100_warm_frac  estimated n100 fraction of CCD pixels for this observation
ccd_temp        mean CCD temperature over 500 seconds surrounding guide transition
known_bad       (boolean) ignore this star in standard processing
bad_comment     reason to ignore a "known_bad" star
=============== ====================================================================

For the columns that reference "corrected" y-angle and z-angle, these have been
corrected by the one-shot quaternion update used during the acquisition sequence.


