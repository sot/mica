Acquisition Data
================

Processing (still done outside mica)
------------------------------------

For each observation, after the observation has run and telemetry is available:

* an independent process is run to ingest the star catalog for the as-run observation

The acquisition and guide stats process

* fetches the AGASC information for each star in the catalog
* fetches the PCAD data at the end of the acquisition interval

and for each acquisition star determines

* if that star was "successfully" acquired
* what the observed magnitude and position of the star were in the last PCAD telemetry
  readout before the guide transition.

Acquisition data products
-------------------------

The acquisition database includes two tables

* acq_stats_data
* acq_stats_warnings

These tables are on the sybase server and there is presently no mica API to access them.

The *acq_stats_data* table includes

========== ========================================================
 Column    Description
========== ========================================================
obsid      the obsid
obi        observation interval number
tstart     observation (obspar) tstart (Chandra secs)
tstop      observation tstart (Chandra secs)
slot       ACA readout slot id
idx        star catalog index id
cat_pos    position in uploaded star catalog
type       star catalog type (BOT or ACQ)
agasc_id   AGASC (catalog) id
obc_id     acquisition success indicator ('ID', or 'NOID')
yang       commanded y-angle position for center of readout window
zang       commanded z-angle position for center of readout window
mag        catalog MAG_ACA of acquisition star
color      catalog COLOR1 of acquisition star
halfw      acquisition search box half width
mag_obs    observed magnitude of star
yang_obs   observed y-angle position of centroid of star
zang_obs   observed z-angle position of centroid of star
y_offset   mean y offset of all of the other acquired stars
z_offset   mean z offset of all of the other acquired stars
d_mag      observed mag - catalog mag
d_yang     observed y-angle - (catalog y + y_offset)
d_zang     observed z_angle - (catalog z + z_offset)
ap_date    last automatic processing date
revision   acq stats processing software version
========== ========================================================


The *acq_stats_warnings* table includes

========== ========================================================
 Column    Description
========== ========================================================
obsid      the obsid
obi        observation interval number
slot       ACA readout slot id
details    text of starcheck warning for slot
========== ========================================================


Acquisition Success
-------------------

The probability of acquisition success can be a function of:

* Star magnitude
* Dark current distribution (probably need to parametrize this, function of time and
* temperature)
* Starcheck warnings like spoiler, common column, bad pixel etc
* Search box size (for acq stars)
* Orbital position / radiation environment
* Others?

Several of these parameters are available within the acquisition statistics database
(magnitude, warnings, search box size) as seen above.

A value for the dark current distribution can presently be obtained using the current dark
current model within mica.archive.aca_dark.dark_model.

A histogram of the dark current distribution maybe obtained with::

  x, xbins, y = mica.archive.aca_dark.dark_model.get_dark_hist(date='2014:001', T_ccd='-14')

The immediately useful value of the estimate of "warm" pixels from this model may be
obtained with::

  warm_frac = mica.archive.aca_dark.dark_model.get_warm_fracs(warm_threshold=100,
                                                            date='2014:001', T_ccd='-14')


Handy Data for Acquisition Success Fitting
------------------------------------------

Code to make a table with acquisition data plus a warm pixel estimate and columns to hold
starcheck warning status::

  import Ska.DBI
  import Ska.Numpy
  from Ska.engarchive import fetch_sci
  import mica.archive.aca_dark.dark_model
  import numpy as np
  import re
  
  db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')
  # Get useful columns from the acq stats table, but don't get
  # ap_date because it is awkward to pickle the sybase datetime later
  cols = ['obsid', 'obi', 'tstart', 'tstop', 'slot', 'idx', 'cat_pos',
          'type', 'agasc_id', 'obc_id', 'yang', 'zang', 'mag', 'color',
          'halfw', 'mag_obs', 'yang_obs', 'zang_obs', 'y_offset', 'z_offset',
          'd_mag', 'd_yang', 'd_zang', 'revision']
  acq_stars = db.fetchall("select {} from acq_stats_data".format(",".join(cols)))
  warnings = db.fetchall("select * from acq_stats_warnings")
  db.conn.close()
  
  # Make a time-map/lookup table of the warm pixel values for the times
  # of the acquisition stars
  warm_pix_tmap = {}
  warm_threshold = 100.
  for time in np.unique(acq_stars['tstart']):
      temp = np.mean(fetch_sci.MSID('AACCCDPT',
                                    time-250,
                                    time+250).vals)
      warm_frac = mica.archive.aca_dark.dark_model.get_warm_fracs(
          warm_threshold, time, temp)
      warm_pix_tmap[time] = warm_frac
  
  # Get a warm pixel value for every star
  warm_pix = []
  for star in acq_stars:
      warm_pix.append(warm_pix_tmap[star['tstart']])
  acq = Ska.Numpy.add_column(acq_stars, 'warm_pix', warm_pix)
  
  # Add columns for some starcheck warnings
  warn_types = ['red_spoiler', 'yellow_spoiler', 'bad_pixel',
                'common_column', 'known_bad_star']
  for k in warn_types:
      acq = Ska.Numpy.add_column(acq, k, np.zeros(len(acq), dtype=bool))
  # Read through warnings and put matching ones in the table
  for warn in warnings:
      acq_match = ((acq['obsid'] == warn['obsid'])
                   & (acq['slot'] == warn['slot'])
                   & (acq['obi'] == warn['obi']))
      if re.search('Search spoiler', warn.details):
          dmag = float(re.search('.*spoil.*\s(\S+)$', warn.details).group(1))
          if (dmag > -0.2):
              acq['red_spoiler'][acq_match] = True
          else:
              acq['yellow_spoiler'][acq_match] = True
      if re.search('Bad Acquisition Star', warn.details):
          acq['known_bad_star'][acq_match] = True
      if re.search('Common Column', warn.details):
          acq['common_column'][acq_match] = True
      if re.search('bad pixel', warn.details):
          acq['bad_pixel'][acq_match] = True
  
  np.save('acq_table.npy', acq)

The resulting numpy save file of the recarray/table would include these columns

============== =======================================================
 Column        Description
============== =======================================================
warm_pix       estimated CCD warm fraction from dark model
red_spoiler    starcheck spoiler with mag difference > -0.2
yellow_spoiler starcheck spoiler with mag difference > -1.0 and < -0.2
bad_pixel      known bad pixel in search box
common_column  star is marked as in common column
known_bad_star acquisition star has already had multiple failures
============== =======================================================

in addition to the columns copied in from the acq_stats_data table:

========== ========================================================
 Column    Description
========== ========================================================
obsid      the obsid
obi        observation interval number
tstart     observation (obspar) tstart (Chandra secs)
tstop      observation tstart (Chandra secs)
slot       ACA readout slot id
idx        star catalog index id
cat_pos    position in uploaded star catalog
type       star catalog type (BOT or ACQ)
agasc_id   AGASC (catalog) id
obc_id     acquisition success indicator ('ID', or 'NOID')
yang       commanded y-angle position for center of readout window
zang       commanded z-angle position for center of readout window
mag        catalog MAG_ACA of acquisition star
color      catalog COLOR1 of acquisition star
halfw      acquisition search box half width
mag_obs    observed magnitude of star
yang_obs   observed y-angle position of centroid of star
zang_obs   observed z-angle position of centroid of star
y_offset   mean y offset of all of the other acquired stars
z_offset   mean z offset of all of the other acquired stars
d_mag      observed mag - catalog mag
d_yang     observed y-angle - (catalog y + y_offset)
d_zang     observed z_angle - (catalog z + z_offset)
revision   acq stats processing software version
========== ========================================================
