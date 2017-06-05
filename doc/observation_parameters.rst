Observation parameters
----------------------

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


And a method is provided to fetch obsids that have obspars:

   >>> obspar.get_obsids('2016:001', '2016:002')
   [51365, 18736, 18036, 18203]


.. _stats-label:

