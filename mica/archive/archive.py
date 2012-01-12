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
last_id_file = os.path.join(archive_dir, 'last_asp1_id.txt')
aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')
apstat = Ska.DBI.DBI(dbi='sybase', server='sqlocc', database='axafapstat')
arc5 = Ska.arc5gl.Arc5gl()
logger = logging.getLogger('asp1 fetch')
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())


def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--obsid",
                      type='int')
    parser.add_option("--version",
                      default='last')
    opt, args = parser.parse_args()
    return opt, args


def get_asp(obsid, version='last'):

    tempdir = tempfile.mkdtemp()
    arc5.sendline("reset")
    arc5.sendline("cd %s" % tempdir)
    logger.info("retrieving data for %d in %s" % (obsid, tempdir))
    arc5.sendline("obsid=%d" % obsid)
    # if multi-obi, we'll need to be specific
    obis = apstat.fetchall(
        "select distinct obi from obidet_0_5 where obsid = %d" % obsid)
    if len(obis) > 1:
        minobi = np.min(obis['obi'])
        logger.info("limiting arc5gl to obi %d" % minobi)
        arc5.sendline("obi=%d" % minobi)
    arc5.sendline("version=%s" % version)
    if version == 'last' or version == 'default':
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
        # update version to number
        version = versions.keys()[0]
    #if version != 'last' and version != 'default' \
    #        and int(got_vers) != int(version):
    #    raise ValueError("Did not receive requested version got %d vs %d"
    #                     % (got_vers, version))
    strobs = "%05d_v%02d" % (obsid, version)
    chunk_dir = strobs[0:2]
    obs_dir = os.path.join(archive_dir, chunk_dir, strobs)
    if not os.path.exists(obs_dir):
        logger.info("making directory %s" % obs_dir)
        os.makedirs(obs_dir)
    else:
        logger.info("obsid dir %s already exists" % obs_dir)
        for afile in glob(os.path.join(tempdir, "*")):
            logger.info("deleting %s" % afile)
            os.remove(afile)
        os.removedirs(tempdir)
        return
    last_obs_ln = os.path.join(archive_dir, chunk_dir, "%05d" % obsid)
    from Ska.Shell import bash
    obs_dirs = glob(os.path.join(archive_dir, chunk_dir, '%05d_v*' % obsid))
    chunk_dir_path = os.path.join(archive_dir, chunk_dir)
    max_obs = sorted(obs_dirs)[-1]
    logger.info("ln -sf %s %s" % (os.path.relpath(max_obs, chunk_dir_path),
                               last_obs_ln))
    bash("ln -sf %s %s" % (os.path.relpath(max_obs, chunk_dir_path),
                           last_obs_ln))
    # if this is new, fetch everything
    arc5.sendline("get asp1")
    logger.info("copying asp1 from %s to %s" % (tempdir, obs_dir))
    for afile in glob(os.path.join(tempdir, "*")):
        shutil.copy(afile, obs_dir)
        os.remove(afile)
    os.removedirs(tempdir)


def main(opt):
    # if an obsid is requested, just do that
    if opt.obsid:
        get_asp(opt.obsid, version=opt.version)
        return
    # if no obsid specified, try to retrieve all asp_1 runs
    # since tool last run
    # use an aspect_1_id in a file
    last_id_fh = open(last_id_file)
    last_id = int(last_id_fh.read().rstrip())
    last_id_fh.close()
    todo = apstat.fetchall("""select * from aspect_1
                              where aspect_1_id > %d
                              order by aspect_1_id"""
                           % last_id)
    for obs in todo:
        logger.info("running get_asp for obsid %d aspect run on %s"
                    % (obs['obsid'], obs['ap_date']))
        get_asp(obs['obsid'], obs['revision'])
    if not len(todo):
        logger.info("No new aspect_1 data")
    else:
        last_id_fh = open(last_id_file, 'w')
        last_id_fh.write("%d" % todo[-1]['aspect_1_id'])
        last_id_fh.close()

if __name__ == '__main__':
    opt, args = get_options()
    main(opt)
