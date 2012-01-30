-----------
Mica Layout
-----------

ObsReport
---------

For a list of topics, collect Topic and App objects.

 - Initialized with obsid
 - Options
   - database-insert-enabled
   - version
   - directory

App
---

 - Plugin with configurable raw HTML (aca_movie CGI?)
   provides menu item and a page frame of content

Topic
-----

object is a collection of

 - Template
 - Table
 - Plot
 - Builder

Template
--------

includes methods to

 - set position of Table and Plot objects 
 - returns HTML or other layout
 - returns menu item

Plot
----   

(given needs, these may just end up being matplotlib figures) includes
methods to

 - receive and return plots, raw data, and annotations
 - write out plots

Table
-----

Table style data with description information and data that can be
used for color/status includes methods to

 - receive and return tables, rows, and annotations
 - insert in database
 - write out as json

Builder
-------

From a set of requests or tests builds Tables and Plots for the
Topic. Structure of Request or Test TBD


Topics that inherit from Topic
------------------------------

Catalog
-------
 - Builder builds catalog Table object 
   - from starcheck_* tables
   - from reference to original starcheck
   - fetches plain plot (rebuild plot as Plot?)  

Guide
-----
 - Builder builds stats Table
   - from trak_stats table

GuideCatalog
------------
 - combines Catalog and Guide objects
   (new tests based on combination?)

Acq
---
 - Builder builds stats Table
   - from acq_stats_data table

AcqCatalog
----------
 - combines Catalog and Acq objects

Thermal
-------
 - Build thermal data from ACA0 Header 3 records
 - Fetch data from AOACCCDPT

Maneuver
--------
 - Build from manuevers from planning ms*sum etc

OneShot
-------
 - retrieve one shot data
 - integrate with Catalog

V&V
---
 - Run task to fetch new ASP1 to archive as needed
 - Builder fetches ASP1 from ASP1 archive
 - Builder fetches Obspar
 - Build Tables for slots overall, Slots, obsid tests
 - Plots per Slot and AspSol
 - Run multirepro for image reconstruction
 - Nest V&V for multirepro

