import os
import Ska.arc5gl
import tempfile
from glob import glob
import re
import logging
import shutil
import Ska.DBI
import numpy as np


archive_dir = "."
aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')
apstat = Ska.DBI.DBI(dbi='sybase', server='sqlocc', database='axafapstat')

def get_asp(obsid, version='last'):

    tempdir = tempfile.mkdtemp()
    arc5 = Ska.arc5gl.Arc5gl()
    arc5.sendline("cd %s" % tempdir)
    arc5.sendline("obsid=%d" % obsid)
    # if multi-obi, we'll need to be specific
    obis = aca_db.fetchall("select * from obspar where obsid = %d" % obsid)
    if len(np.unique(obspar['obi'])) > 1:
        minobi = np.min(obspar['obi'])
        arc5.sendline("obi=%d" % minobi)
    arc5.sendline("version=%s" % version)
    # just fetch the aspect solution to start
    arc5.sendline("get asp1{aspsol}")
    asols = glob(os.path.join(tempdir, "*asol*fits*"))
    if not asols:
        raise ValueError("No asol files in %s" % tempdir)
    versions = {}
    for asol in asols:
        fmatch = re.search('pcadf\d+N(\d{3})_asol1', asol)
        if fmatch:
            versions[int(fmatch.group(1))] = 1
    if len(versions) > 1:
        raise ValueError("Different version asol files in %s" % tempdir)
    got_vers = versions.keys()[0]
    if version != 'last' and version != 'default' \
            and int(got_vers) != int(version):
        raise ValueError("Did not receive requested version got %d vs %d"
                         % (got_vers, version))
    strobs = "%05d_v%02d" % (obsid, got_vers)
    chunk_dir = strobs[0:2]
    obs_dir = os.path.join(archive_dir, chunk_dir, strobs)
    if not os.path.exists(obs_dir):
        os.makedirs(obs_dir)
    else:
        for afile in glob(os.path.join(tempdir, "*")):
            os.remove(afile)
        os.removedirs(tempdir)
        return
    # if this is new, fetch everything
    arc5.sendline("get asp1")
    for afile in glob(os.path.join(tempdir, "*")):
        shutil.copy(afile, obs_dir)
        os.remove(afile)
    os.removedirs(tempdir)


last_id = 22011
todo = apstat.fetchall("select * from aspect_1 where aspect_1_id > %d"
                       % last_id)
for obs in todo:
    get_asp(obs['obsid'])

    
