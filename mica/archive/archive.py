import os
import Ska.arc5gl
import tempfile
from glob import glob
import re
import logging
import shutil
import Ska.DBI
from Ska.Shell import bash
import numpy as np


archive_dir = "."
last_id_file = os.path.join(archive_dir, 'last_asp1_id.txt')
aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')
apstat = Ska.DBI.DBI(dbi='sybase', server='sqlocc', database='axafapstat')
arc5 = Ska.arc5gl.Arc5gl()
#arc5.echo = True
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


def get_file_ver(tempdir, fileglob="*fidpr1*fits*"):
    files = glob(os.path.join(tempdir, fileglob))
    if not files:
        return None
    versions = {}
    for f in files:
        fmatch = re.search('pcadf\d+N(\d{3})_', f)
        if fmatch:
            versions[int(fmatch.group(1))] = 1
        if len(versions) > 1:
            raise ValueError("Different version files in %s" % tempdir)
    # update version to number
    version = versions.keys()[0]
    return version


def get_ver_num(obsid, version='default'):
    # this is obi agnostic, all obis should be same
    # version anyway...
    tempdir = tempfile.mkdtemp()
    print tempdir
    arc5.sendline("reset")
    arc5.sendline("cd %s" % tempdir)
    arc5.sendline("obsid=%d" % obsid)
    if version != 'default':
        arc5.sendline("version=%s" % version)
    # just fetch the aspect solution to start
    arc5.sendline("get asp1{fidprops}")
    version = get_file_ver(tempdir)
    for afile in glob(os.path.join(tempdir, "*")):
        os.remove(afile)
    os.removedirs(tempdir)
    return version


def get_asp(obsid, version='last'):
    n_version = get_ver_num(obsid, version=version)
    if n_version is None:
        raise ValueError("No ASP1 for ver %s" % version)
    # use numeric version instead of 'last' or 'default'
    version = n_version

    # set up directory for data
    strobs = "%05d_v%02d" % (obsid, version)
    chunk_dir = strobs[0:2]
    obs_dir = os.path.join(archive_dir, chunk_dir, strobs)
    if not os.path.exists(obs_dir):
        logger.info("making directory %s" % obs_dir)
        os.makedirs(obs_dir)
    else:
        logger.info("obsid dir %s already exists" % obs_dir)
        return

    # get data
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
    arc5.sendline("get asp1")
    logger.info("copying asp1 from %s to %s" % (tempdir, obs_dir))
    for afile in glob(os.path.join(tempdir, "*")):
        shutil.copy(afile, obs_dir)
        os.remove(afile)
    os.removedirs(tempdir)
    return obs_dir, chunk_dir


# this needs help
def update_link(obsid):
    # for the obsid get the default and last version numbers
    default_ver = get_ver_num(obsid, 'default')
    last_ver = get_ver_num(obsid, 'last')

    # set up directory for data
    chunk_dir = ("%05d" % obsid)[0:2]

    obs_ln = os.path.join(archive_dir, chunk_dir, "%05d" % obsid)

    obs_ln_last = os.path.join(archive_dir, chunk_dir,
                               "%05d_last" % obsid)

    if default_ver is not None:
        def_ver_dir = os.path.join(archive_dir, chunk_dir,
                                   '%05d_v%02d' % (obsid, default_ver))
        # link the default version
        chunk_dir_path = os.path.join(archive_dir, chunk_dir)
        logger.info("ln -sf %s %s" % (
                os.path.relpath(def_ver_dir, chunk_dir_path),
                obs_ln))
        bash("ln -sf %s %s" % (
                os.path.relpath(def_ver_dir, chunk_dir_path),
                obs_ln))

    if last_ver is not None:
        last_ver_dir = os.path.join(archive_dir, chunk_dir,
                                    '%05d_v%02d' % (obsid, last_ver))
        # if last and default are different, and we have last
        # make a link to it
        if last_ver != default_ver:
            if os.path.exists(last_ver_dir):
                obs_ln_last = os.path.join(archive_dir, chunk_dir,
                                           "%05d_last" % obsid)
                logger.info("ln -sf %s %s" % (
                        os.path.relpath(last_ver_dir, chunk_dir_path),
                        obs_ln_last))
                bash("ln -sf %s %s" % (
                        os.path.relpath(last_ver_dir, chunk_dir_path),
                        obs_ln_last))
        else:
            # if default and last are the same and we have a last link
            # delete it
            if os.path.exists(obs_ln_last):
                logger.info("removing outdated link %s"
                            % obs_ln_last)
                os.remove(obs_ln_last)


def get_prov_data():
    # find obsids that have a "_last" link and check to see if
    # they have been updated to have new
    # flight/approved/released/default data
    chunk_dirs = glob(os.path.join(archive_dir, "??"))
    todo_obs = []
    for cdir in chunk_dirs:
        last_links = glob(os.path.join(cdir, "*_last"))
        for link in last_links:
            lmatch = re.search('(\d{5})_last', link)
            if lmatch:
                obs = dict(obsid=int(lmatch.group(1)),
                           revision='default')
                todo_obs.append(obs)
    return todo_obs


def main(opt):
    # if an obsid is requested, just do that
    if opt.obsid:
        get_asp(opt.obsid, version=opt.version)
        update_link(opt.obsid)
        return

    # look for all previous "provisional" aspect solutions
    # and update if possible
    prov_data = get_prov_data()
    for obs in prov_data:
        logger.info("running get_asp for prov obsid %d ver %s, "
                    % (obs['obsid'], obs['revision']))
        get_asp(obs['obsid'], obs['revision'], newonly=True)
        update_link(obs['obsid'])

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
        # let's not pass the revision and thus just get
        # released data when run in update mode
        get_asp(obs['obsid'], version='default')
        update_link(obs['obsid'])
    if not len(todo):
        logger.info("No new aspect_1 data")
    else:
        last_id_fh = open(last_id_file, 'w')
        last_id_fh.write("%d" % todo[-1]['aspect_1_id'])
        last_id_fh.close()


if __name__ == '__main__':
    opt, args = get_options()
    main(opt)
