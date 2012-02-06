
Example: Acquisition Report
----------------------------
::

  acq = mica.acq(obsid=?????)

The acq function would run the acq tests on the obsid acquisition
statistics data, general acquisitions statistics, and starcheck catalog
data and return a container object with everything.  Everything would
include:

* Annotated data tables from each type of included data
* Annotated test table
* A starcheck plot
* A test result plot
* The data used to make those plots

Data
====

Acq Stats Table Data includes, for each obsid / obi / slot:

::

 [102] .aca.1> sp_cols acq_stats_data;
 name         type         length     
 ------------ ------------ -----------
 obsid        int                    4
 obi          int                    4
 tstart       float                  8
 tstop        float                  8
 slot         int                    4
 idx          int                    4
 cat_pos      int                    4
 type         varchar                3
 agasc_id     int                    4
 obc_id       varchar                6
 yang         float                  8
 zang         float                  8
 mag          float                  8
 color        float                  8
 halfw        float                  8
 mag_obs      float                  8
 yang_obs     float                  8
 zang_obs     float                  8
 y_offset     float                  8
 z_offset     float                  8
 d_mag        float                  8
 d_yang       float                  8
 d_zang       float                  8
 ap_date      datetime               8
 revision     varchar               10

:cat_pos:
  Position of entry in catalog

:idx:
  Index of entry in starcheck

:type:
  ACQ BOT GUI FID or MON

:agasc_id:
  AGASC id or fid id

:obc_id:
  ID or NOID

:yang, zang, mag, color:
  Planned y-angle, z-angle, and catalog mag and color

:mag_obs, yang_obs, zang_obs:
  Observed mag, y-angle, and z_angle

:y_offset, z_offset:
  mean y and z offset of the observation over the acquired stars with valid magnitudes

:d_mag:
  Observed mag - catalog mag

:d_yang:
  Observed y-angle - (planned y-angle + y_offset)

:d_zang:
  Observed z-angle - (planned z-angle + z_offset)

:ap_date:
  Last signed-off obidet_0_5 processing date

:revision:
  Last version of acq_stat processing tool


Tests
=====
  

Tests could include:

* on a per-star basis:

  * was this star acquired
  * was this star part of an acquisition anomaly
  * did we acquire this star in an earlier observation but now failed
  * is the star delta magnitude outside the 95% percentile

* on a per-obsid basis:

  * does this obsid have an acquisition anomaly
  * does this obsid have an acquisition failure
  * how many

The test directive would also mark any features of interest on a
new star catalog plot extended from the starcheck plot.

The acq test template includes directives for how these data may be
organized into an html report and how they may be placed in a file
archive/cache.  Each test/entry also includes a directive for how it
should generalize itself if called for a collection of obsids.

Example: Collective Acquisition Report
---------------------------------------
::

  acq = mica.collective_acq(obsid_list)

The collective_acq template instructs collection of acq objects for the
obsids in obsid_list, gathers those data, and provides methods to at least
recreate the current acq_stats_reports for an arbitrary set of obsids.

There should exist helper tools in the set to extract a set of obsids
based on criteria, as in:

* display complete reports for all observations that use this star

