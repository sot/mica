#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import argparse

from mica.stats import update_guide_stats
import mica.stats.guide_stats

# Cheat and pass options directly.  Needs entrypoint scripts
opt = argparse.Namespace(
    datafile=mica.stats.guide_stats.TABLE_FILE,
    obsid=None,
    check_missing=False,
    start=None,
    stop=None,
)
update_guide_stats.update(opt)


table_file = mica.stats.guide_stats.TABLE_FILE
file_stat = os.stat(table_file)
if file_stat.st_size > 200e6:
    print(
        """
Warning: {tfile} is larger than 200MB and may need
Warning: to be manually repacked (i.e.):
Warning:
Warning: ptrepack --chunkshape=auto --propindexes --keep-source-filters {tfile} compressed.h5
Warning: cp compressed.h5 {tfile}

""".format(tfile=table_file)
    )
