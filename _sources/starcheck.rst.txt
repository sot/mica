.. _starcheck-label:

Starcheck products
==================

The :mod:`mica.starcheck` module provides tools to update a database of starcheck products
and methods to retrieve content from that database.  The `mica.starcheck` database
inherits much of the structure of the original Sybase starcheck_database project which
organized the per-obsid data into these tables:

   * starcheck_obs
   * starcheck_manvr
   * starcheck_warnings
   * starcheck_catalog

That Sybase database was updated for new obsids as they became available in the CXCDS
archive and is intended to contain only records for observations in their final/observed
form. When developing the `mica.report` system, we had a new driver to be able to track
observations in the planning stages, so a new database of starcheck products, this
`mica.starcheck` database, was developed which contains entries that include observations
in released schedules that have not yet been observed.  The `mica.starcheck` database
includes all versions of observations that came out in released products and does not
attempt to pre-filter on final/observed configuration on ingest.  These distinctions may
be useful for developers who wish to use the database directly.  For most users, the
included methods filter the records available in the database to return useful values.

The `mica.starcheck` sqlite3 database, as mentioned inherits most of the per-obsid tables,
but adds to each an id `sc_id` that maps the entry back to a new `starcheck_id` table
which has one record for each directory that has been ingested (each directory is assumed
to have a single starcheck.txt file as the source for that entry).  A new table,
`starcheck_pred_temp` has also been added which includes the ACA CCD temperature
prediction that was written directly into starcheck.txt at the time of review (when applicable).


Get catalog data
----------------

The `get_starcheck_catalog` and `get_starcheck_catalog_at_date` methods return
dictionaries of the content of the starcheck database tables.

   >>> from mica.starcheck import get_starcheck_catalog
   >>> obsid_starcheck = get_starcheck_catalog(5438)
   >>> obsid_starcheck.keys()
   ['status', 'warnings', 'cat', 'mp_dir', 'obs', 'manvr']
   >>> obsid_starcheck['cat']
   <Table masked=False length=11>
   sc_id obsid obs_idx    mp_starcat_time    ...  res  halfw    pass    notes 
   int64 int64  int64        unicode672      ... int64 int64 unicode128 object
   ----- ----- ------- --------------------- ... ----- ----- ---------- ------
     941  5438       2 2005:282:18:53:05.394 ...     1    25              None
     941  5438       2 2005:282:18:53:05.394 ...     1    25              None
     941  5438       2 2005:282:18:53:05.394 ...     1    25              None
     941  5438       2 2005:282:18:53:05.394 ...     1   120       a2g2   None
     941  5438       2 2005:282:18:53:05.394 ...     1   120       a2g2   None
     941  5438       2 2005:282:18:53:05.394 ...     1   120         a3   None
     941  5438       2 2005:282:18:53:05.394 ...     1   120       a2g3   None
     941  5438       2 2005:282:18:53:05.394 ...     1   120              None
     941  5438       2 2005:282:18:53:05.394 ...     1   120         a2   None
     941  5438       2 2005:282:18:53:05.394 ...     1   120         a2   None
     941  5438       2 2005:282:18:53:05.394 ...     1   120         a3   None
   >>> obsid_starcheck['cat'].colnames
   ['sc_id', 'obsid', 'obs_idx', 'mp_starcat_time', 'idx', 'slot', 'id', 'idnote', 'type',
   'sz', 'minmag', 'mag', 'maxmag', 'yang', 'zang', 'dim', 'res', 'halfw', 'pass',
   'notes']

The `get_starcheck_catalog_at_date` routine is designed to retrieve the catalog that would
apply at the given date.  "Apply" in this context isn't quite the "onboard catalog at the
time" but instead, if the given date is during NMM is the catalog that was commanded
during that NMM, or if the date is in NPNT, the catalog is the catalog that was commanded
in the previous NMM (that should be tracked during the NPNT interval).

`get_starcheck_catalog_at_date` is intended both to help with larger scale processing (to
get the "right" catalog for a whole bunch of kadi dwells) and for some special cases. For
example, it will get the commanded catalog during a vehicle-only period of time on the
spacecraft.

   >>> from mica.starcheck import get_starcheck_catalog_at_date
   >>> catalog_by_date = get_starcheck_catalog_by_date('2015:174:01:25:30.068')
   >>> catalog_by_date['obs']['obsid']
   16689

   >>> from kadi import events
   >>> dwells = events.dwells.filter('2015:174:01:25:30.068', '2015:174:01:25:30.068')
   >>> dwells[0].get_obsid()
   17696

This is an example of a time in which the commanded obsid (COBSQID) was not updated due to
SCS107, so a fetch of the catalog for obsid 17696 will not give the right star positions
if one desires to calculate residuals during the dwell.

