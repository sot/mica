.. mica documentation master file, created by
   sphinx-quickstart on Sun Aug  5 16:28:18 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Mica Documentation
================================

Contents:

.. toctree::
   :maxdepth: 2

Mica Archive
============

The Archive components of the mica suite provide tools to:

 * retrieve telemetry and processed products from the CXCDS archive
 * store these products in a Ska file archive for speed and convenience
 * quickly retrieve these products 

Mica provides interfaces to these products:

 * ACA L0 telemetry
    * ACA diagnostic telemetry (header three)
 * Aspect L1 products
 * Obspars 

ACA L0 module
----------------

The :mod:`mica.archive.aca_l0` module provides tools to build and fetch from
a file archive of ACA L0 telemetry.  This telemetry is stored in
directories by year and day-of-year, and ingested filenames are stored
in a lookup table.  

Methods are provided to retrieve files and read those data files into
data structures.

   >>> from mica.archive import aca_l0

   >>> obsid_files = aca_l0.get_files(obsid=5438)
   >>> time_files = aca_l0.get_files(start='2012:001', stop='2012:002')
   >>> time_8x8 = aca_l0.get_files(start='2011:001', stop='2011:010',
   ...                             imgsize=[8])


   >>> slot_data = aca_l0.get_slot_data('2012:001', '2012:002', slot=7)
   >>> temp_ccd_8x8 = aca_l0.get_slot_data('2005:001', '2005:010',
   ...                                     slot=6, imgsize=[8],
   ...                                     columns=['TIME', 'TEMPCCD'])

The :func:`~mica.archive.aca_l0.get_slot_data()` method will retrieve all columns by default and
the resulting data structure will have masked columns where those
values are not available (i.e. HD3TLM64 in 6x6 or 4x4 image data).


ACA Diagnostic module
---------------------

The :mod:`mica.archive.aca_hdr3` module works with Header 3 data from 8x8 ACA
L0 image data.  An MSID class and MSIDset class are provided to fetch 
pseudo-MSIDs and return masked array data structures to work with them.

   >>> from mica.archive import aca_hdr3
   >>> ccd_temp = aca_hdr3.MSID('ccd_temp', '2012:001', '2012:020')
   >>> type(ccd_temp.vals)
   'numpy.ma.core.MaskedArray'
   >>> perigee_data = aca_hdr3.MSIDset(['ccd_temp', 'aca_temp', 'dac'],
   ...                                 '2012:001', '2012:030')


Aspect L1 module
------------------

The :mod:`mica.archive.asp_l1` module provides tools to build and fetch from
a file archive of Aspect level 1 products.

Methods are provided to find the archive directory:

   >>> from mica.archive import asp_l1
   >>> asp_l1.get_dir(2121)
   '/data/aca/archive/asp1/02/02121'
   >>> obsdirs = asp_l1.get_obs_dirs(6000)

The obsdirs dictionary should look something like::

  {'default': '/data/aca/archive/asp1/06/06000',
  2: '/data/aca/archive/asp1/06/06000_v02',
  3: '/data/aca/archive/asp1/06/06000_v03',
  'last': '/data/aca/archive/asp1/06/06000',
  'revisions': [2, 3]}

Methods are also provided to retrieve a list files by obsid and time range.

   >>> obs_files = asp_l1.get_files(6000)
   >>> obs_gspr = asp_l1.get_files(6000, content=['GSPROPS'])
   >>> range_fidpr = asp_l1.get_files(start='2012:001',
   ...                                stop='2012:030',
   ...                                content=['FIDPROPS'])




Obspar module
------------------

The :mod:`mica.archive.obspar` module provides tools to build and fetch from
a file archive of obspars.

Methods are provided to find the archive directory and obspar files:


   >>> obspar.get_obspar_file(7000)
   '/data/aca/archive/obspar/07/07000/axaff07000_000N002_obs0a.par.gz'
   >>> obspar.get_dir(2121)
   '/data/aca/archive/obspar/02/02121'
   >>> obsdirs = obspar.get_obs_dirs(6000)

The obsdirs dictionary should look something like::

   {'default': '/data/aca/archive/obspar/06/06000',
   2: '/data/aca/archive/obspar/06/06000_v02',
   3: '/data/aca/archive/obspar/06/06000_v03',
   'last': '/data/aca/archive/obspar/06/06000',
   'revisions': [2, 3]}

A method is provided to read the obspar into a dictionary:

   >>> from mica.archive import obspar
   >>> obspar.get_obspar(7001)['detector']
   'ACIS-I'

Methods are also provided to retrieve a list of files by obsid and time range.

   >>> obs_files = obspar.get_files(6000)
   >>> range = obspar.get_files(start='2012:001',
   ...                          stop='2012:030')




API docs
========

.. toctree::
   :maxdepth: 1

   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

