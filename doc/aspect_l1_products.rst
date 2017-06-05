Aspect L1 products
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


