Mica is a suite of tools that is designed to provide convenient access and
display of all Aspect-relevant data for an obsid or collection of obsids.

As such, the suite provides interfaces to:

* aggregate data on a per-obsid basis
* aggregate per-obsid data based on obsid search criteria?
* return these data as meaningful structures to other tools
* generate customizable, human-readable report or reports
* archive or cache obsid data for quick re-query-ing

The suite aggregates all available aspect planning and trending data and these
data will also be extended as needed to answer questions to best support operations.

The current aspect trending tools consist of a variety of report-based
projects.  Most of these tools are designed to gather data on a per-obsid
basis, but the tools used to visualize the data are often just over
monthly (or longer) time ranges.  To evaluate an individual obsid, one
must manually examine several database tables, possibly retrieve a V&V
report using the web interface, and then retrieve and examine telemetry.
Mica does all of this (and MORE!).

Data Sources
==============

Planning data
--------------

Starcheck outputs.  
  The starcheck outputs are already stored in sybase tables (starcheck_catalog,
  starcheck_manvr, starcheck_obs, starcheck_warnings).  The timelines table
  also maps the obsid to its SOT MP directory to get a link to the original
  products (including the plots).

Maneuver sequence data  
  The maneuvers to and from the obsid target
  should be presented in context.  The ms*sum file in the planned products
  directory should be parsed for these data.  A new tool should archive
  these data in sybase or using a file archive.

cmd_states outputs
  The predicted command states for the obsid will be available in the
  cmd_states table.

From a processing perspective, the aggregation tool will need to be able
to grab planning data when loads are approved, and should provide access
and display of these pieces immediately.

Misc data
----------

Processing Status

While not perfectly accurate, the queries used to make the processing
status web page will be used to make a table of current status.  This
status is dynamic and thus won't be archived / cached.

L0 data
---------

ACA img0 data
 *Define structure and content of ACA img0 archive*.  (If it is needed
  at all... how big would it be?)

L0 derived products
--------------------

observations_table / obspar data
  The observations table and obspar tables are available in the aca sybase
  database. These are flat tables.

acq and tracking stats
  The acq_stats_data table and trak_stats_data tables contain acquisition,
  guide star tracking, and fid light tracking statistics.  These are
  archived on a per-obsid basis and are also flat tables.  These data
  should be integrated with the catalog data for report construction.

aonstar aokalstr
  These data are are available in telemetry but not presently gathered on a
  per-obsid basis.  For aggregation by Mica these data may be gathered
  directly from the engineering archive, or archived in sybase tables by a
  to-be-developed pre-processing tool.

temperatures / header 3 and those in regular telemetry
  ACA temperatures such as the CCD temperature AACCCDPT are available in
  telemetry from the engineering archive.  ACA hdr 3 telemetry provides
  higher resolution data and is presently archived in the perigee_health
  processing.  That processing also makes pickle and json files of the data
  for each pass, though the json files may be of less use for per-obsid
  aggregation.

aca_bgd_mon if available
  The aca bgd data has dates for high background events.  For aggregation on
  a per-obsid basis it would be useful to include a yes/no and link to high
  background plots if available.

periscope data
  The obs_periscope_tilt table includes per-obsid periscope change
  statistics

obc_rate_noise data
  ???

L1 data
------------

Raw aspect L1 data products.  *Define structure and content of L1 archive*.

L1 derived products
--------------------

fid_drift data
 ...

astromon data
 ...

vv data
  Existing V&V data will not be accessed by the aggregation tool.  Instead,
  a new v&v tool will be run on new data to confirm established V&V limits
  and also allow the convenient creation of new tests, archiving and
  querying of useful metrics, and interactive plotting for V&V
  investigation.  The new v&v tool will require for source data only access
  to the new local ASP1 archive and an obspar for each obsid query.

Not Per-Obsid Data
-------------------

It may be natural to ask questions of these data sets that compare the
per-obsid data for a single obsid to averages or percentiles of all other
data of that type.  New archiving methods for V&V data etc should consider
this.  This will likely be bootstrapped.

Data Types
===========

Will we need any web scraping or the like to assist with aggregation of
obsid data?  That doesn't appear to be the case.  Necessary data appears
to be available as either:

* entries in an already-available sybase table
* raw fits files
* engineering archive data
* task data from perigee_health_plots and aca_bgd_mon
* combinations of the above
* testing output for tests based on above data


Interface
============

The suite will provide uniform access functions over the varied
data products listed above.  For each predefined topic area or test
collection, a configuration file will exist to assist the tool to collect
the appropriate data, run any tests as needed, and return a data structure
to the aggregation tool.  The user shall be able to request the result of
individual topic areas or tests directly, by instantiating those
collection objects, or by calling a general function that retrieves all
defined items.  New test or topic definitions may also be passed directly
to the aggregator.

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
