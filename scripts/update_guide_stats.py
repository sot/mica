#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import mica.stats.guide_stats
mica.stats.guide_stats.main()

import os
table_file = mica.stats.guide_stats.TABLE_FILE
file_stat = os.stat(table_file)
if file_stat.st_size > 200e6:
    print """
Warning: {tfile} is larger than 200MB and may need
Warning: to be manually repacked (i.e.):
Warning:
Warning: ptrepack --chunkshape=auto --propindexes --keep-source-filters {tfile} compressed.h5
Warning: cp compressed.h5 {tfile}

""".format(tfile=table_file)
