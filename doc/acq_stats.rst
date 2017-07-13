Acquisition star statistics
---------------------------
The :mod:`mica.stats.acq_stats` module
includes code to gather acquisition star statistics data for each observation and
return those data to the user.

To get the whole acquisition table data:

   >>> from mica.stats.acq_stats import get_stats
   >>> stats = get_stats()
   >>> stats[(stats['obsid'] == 5438) & (stats['slot'] == 1)][0]['dy']
   8.4296239297976854

The hdf5 in-kernel searches may be faster working with the table directly for some
operations.


Data table fields
^^^^^^^^^^^^^^^^^
The acquisition statistics data table contains the following columns:

=============== ====================================================================
 Column         Description
=============== ====================================================================
obsid           obsid
obi             observation interval number
acq_start       acquisition datestart from kadi event
guide_start     guide transition datestart from kadi event
guide_tstart    guide transition start from kadi event in Chandra secs
one_shot_length one shot delta quaternion arc length for this obsid manvr (arcsec)
revision        revision string representing the software version
slot            ACA readout slot id
idx             starcheck star catalog index id
type            star catalog type (BOT or ACQ)
yang            commanded y-angle position for center of readout window (arcsec)
zang            commanded z-angle position for center of readout window (arcsec)
halfw           acquisition search box half width (arcsec)
mag             catalog MAG_ACA of acquisition star
acqid           acquisition success indicator (boolean)
star_tracked    something tracked within 5 arcsecs during acq interval (boolean)
spoiler_tracked something tracked outside 5 arcsecs during acq interval (boolean)
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
cdy             "corrected" delta-y at guide transition. (arcsec)
cdy             "corrected" delta-z at guide transition (arcsec)
dy              delta-y at guide transition. expected position determined using onboard est. attitude (arcsec)
dz              delta-z at guide transition. expected position determined using onboard est. attitude (arcsec)
ion_rad         ionizing radiation flag set at transition (boolean)
def_pixel       defective pixel flag set at transition (boolean)
mult_star       multiple star flag set at transition (boolean)
sat_pix         saturated pixel flag set at transition (boolean)
mag_obs         observed magnitude at transition
yang_obs        observed y-angle at transition (arcsec)
zang_obs        observed z-angle at transition (arcsec)
agasc_id        AGASC (catalog) id
color1          AGASC COLOR1, estimated B-V color
ra              AGASC right ascension (degrees)
dec             AGASC declination (degrees)
epoch           AGASC EPOCH
pm_ra           AGASC PM_RA, proper motion in ra (milli-arcsec/year)
pm_dec          AGASC PM_DEC, proper motion in dec (milli-arsec/year)
var             AGASC VAR, known or suspected variable star
pos_err         AGASC POS_ERR, position error (milli-arcsec)
mag_aca         AGASC MAG_ACA, ACA mag
mag_err         AGASC MAG_ERR, (to be fixed, this should have been MAG_ACA_ERR)
mag_band        AGASC MAG_BAND, integer code for spectral band
pos_catid       AGASC POS_CATID, integer code for position catalog
aspq1           AGASC ASPQ1, integer spoiler code
aspq2           AGASC ASPQ2, integer proper motion flag
aspq3           AGASC ASPQ3, integer dist in 100m-arcsec to nearest Tycho2 star
acqq1           AGASC ACQQ1, mag diff to brightest star within 53.3" (0.01mags)
acqq2           AGASC ACQQ2, mag diff to brightest star within 107" (0.01mags)
acqq4           AGASC ACQQ4, mag diff to brightest star within 267.5" (0.01mags)
n100_warm_frac  estimated n100 fraction of CCD pixels for this observation
ccd_temp        mean CCD temperature over 500 seconds surrounding guide transition
known_bad       ignore this star in standard processing (boolean)
bad_comment     reason to ignore a "known_bad" star
=============== ====================================================================

For the columns that reference "corrected" y-angle and z-angle, these have been
corrected by the one-shot quaternion update used during the acquisition sequence.

Processing
^^^^^^^^^^

For each observation, after the observation has run and telemetry is available, the
acquisition and guide stats process does the following:

* Fetches the AGASC information for each star in the catalog
* Fetches the PCAD data at the end of the acquisition interval
* For each acquisition star determines:

  * If that star was "successfully" acquired
  * What the observed magnitude and position of the star were in the last PCAD telemetry
    readout before the guide transition.


