#!/usr/bin/env python

import mica.stats.acq_stats
mica.stats.acq_stats.update()

import os
table_file = mica.stats.acq_stats.table_file
file_stat = os.stat(table_file)
if file_stat.st_size > 50e6:
    print """
Warning: {tfile} is larger than 50MB and may need
Warning: to be manually repacked (i.e.):
Warning:
Warning: ptrepack --chunkshape=auto --propindexes --keep-source-filters {tfile} compressed.h5
Warning: cp compressed.h5 {tfile}

""".format(tfile=table_file)
