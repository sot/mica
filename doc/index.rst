.. mica documentation master file, created by
   sphinx-quickstart on Sun Aug  5 16:28:18 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to mica's documentation!
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


ACA L0 telemetry
----------------

The mica.archive.aca_l0 module provides tools to build and fetch from
a file archive of ACA L0 telemetry.  This telemetry is stored in
directories by year and day-of-year, and ingested filenames are stored
in a lookup table.  

ACA L0 module
~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   mica.archive.aca_l0

Retrieving Data from the ACA L0 Archive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Module examples::

  import mica.archive.aca_l0
  interval_files = mica.archive.aca_l0.get_interval_files('2012:001', '2012:002')
  slot_data = mica.archive.aca_l0.get_slot_data('2012:001', '2012:002', slot=7)

ACA Diagnostic module
~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   
   mica.archive.aca_hdr3

Retrieving Data with the ACA diagnostic module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Module examples::

  import mica.archive.aca_hdr3
  ccd_temp = mica.archive.aca_hdr3.MSID('ccd_temp', '2012:001', '2012:020')
  perigee_data = mica.archive.aca_hdr3.MSIDset(['ccd_temp','aca_temp', 'dac'],
                                                     '2012:001', '2012:030')

Aspect L1 products
------------------

ASP L1 module
~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   mica.archive.asp_l1

Retrieving Data with the ASP L1 module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example::

  import mica.archive.asp_l1
  obsdir = mica.archive.asp_l1.get_dir(2121)

Obspar products
------------------

Obspar module
~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   mica.archive.obspar
  
Retrieving Data with the obspar module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example::

  import mica.archive.obspar
  obspar = mica.archive.obspar.get_obspar(2121)




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

