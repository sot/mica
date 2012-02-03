
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

