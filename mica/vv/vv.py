from __future__ import division
import os
import json
import tables
import shutil
import logging
from glob import glob
import numpy as np
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

import Ska.DBI

import mica.archive.asp_l1 as asp_l1_arch
import mica.archive.obspar as obspar_arch
from mica.archive import obsid_archive
from .core import Obi

# get rid of the black edges on the plot markers
plt.rcParams['lines.markeredgewidth'] = 0

mica_archive = os.environ.get('MICA_ARCHIVE') or '/data/aca/archive'
FILES = dict(
    data_root=os.path.join(mica_archive, 'vv'),
    temp_root=os.path.join(mica_archive, 'tempvv'),
    shelf_file=os.path.join(mica_archive, 'vv', 'vv_shelf.db'),
    h5_file=os.path.join(mica_archive, 'vv', 'vv.h5'),
    last_file=os.path.join(mica_archive, 'vv', 'last_id.txt'),
    asp1_proc_table=os.path.join(mica_archive, 'asp1', 'processing_asp_l1.db3'),
    bad_obsid_list=os.path.join(mica_archive, 'vv', 'bad_obsids.json'))

KNOWN_BAD_OBSIDS = json.loads(open(FILES['bad_obsid_list']).read())
#    save=True,
#    known_bad_obsids=[1283, 5542, 722, 721],

logger = logging.getLogger('vv')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def get_vv_dir(obsid, version="default"):
    num_version = None
    if version == 'last' or version == 'default':
        asp_l1_proc = Ska.DBI.DBI(dbi="sqlite", server=FILES['asp1_proc_table'])
        if version == 'default':
            obs = asp_l1_proc.fetchall("""select * from aspect_1_proc
                                          where obsid = {} and isdefault = 1
                                       """.format(obsid))
            if not len(obs):
                raise ValueError("Version {} not found for obsid {}".format(
                    version, obsid))
            num_version = obs['revision'][0]
        if version == 'last':
            obs = asp_l1_proc.fetchall("""select * from aspect_1_proc
                                          where obsid = {}
                                       """.format(obsid))
            if not len(obs):
                raise ValueError("Version {} not found for obsid {}".format(
                    obsid, version))
            num_version = np.max(obs['revision'])
    else:
        num_version = version
    strobs = "%05d_v%02d" % (obsid, num_version)
    chunk_dir = strobs[0:2]
    chunk_dir_path = os.path.join(FILES['data_root'], chunk_dir)
    obs_dir = os.path.join(chunk_dir_path, strobs)
    return obs_dir

def get_vv_files(obsid, version="default"):
    vv_dir = get_vv_dir(obsid, version)
    return glob(os.path.join(vv_dir, "*"))

def get_vv(obsid, version="default"):
    vv_dir = get_vv_dir(obsid, version)
    json_file = glob(os.path.join(vv_dir, "*.json"))[0]
    return json.loads(open(json_file).read())

def get_rms_data():
    h5f = tables.openFile(FILES['h5_file'], 'r')
    tbl = h5f.getNode('/', 'vv')
    data = tbl[:]
    h5f.close()
    return data

def file_vv(obi):
    obsid = int(obi.info()['obsid'])
    version = int(obi.info()['revision'])
    # set up directory for data
    strobs = "%05d_v%02d" % (obsid, version)
    chunk_dir = strobs[0:2]
    chunk_dir_path = os.path.join(FILES['data_root'], chunk_dir)
    obs_dir = os.path.join(chunk_dir_path, strobs)
    if not os.path.exists(obs_dir):
        logger.info("making directory %s" % obs_dir)
        os.makedirs(obs_dir)
    else:
        logger.info("obsid dir %s already exists" % obs_dir)
    for f in glob(os.path.join(obi.tempdir, "*")):
        os.chmod(f, 0775)
        shutil.copy(f, obs_dir)
        os.remove(f)
    logger.info("moved VV files to {}".format(obs_dir))
    os.removedirs(obi.tempdir)
    logger.info("removed directory {}".format(obi.tempdir))
    # make any desired link
    obs_ln = os.path.join(FILES['data_root'], chunk_dir, "%05d" % obsid)
    obs_ln_last = os.path.join(
        FILES['data_root'], chunk_dir, "%05d_last" % obsid)
    obsdirs = asp_l1_arch.get_obs_dirs(obsid)
    isdefault = 0
    if 'default' in obsdirs:
        if (os.path.realpath(obsdirs[version])
                == os.path.realpath(obsdirs['default'])):
            if os.path.islink(obs_ln):
                os.unlink(obs_ln)
            os.symlink(os.path.relpath(obs_dir, chunk_dir_path), obs_ln)
            isdefault = 1
    if 'last' in obsdirs:
        if ('default' in obsdirs
                and (os.path.realpath(obsdirs['last'])
                     != os.path.realpath(obsdirs['default']))
                or 'default' not in obsdirs):
            if (os.path.realpath(obsdirs[version])
                    == os.path.realpath(obsdirs['last'])):
                if os.path.islink(obs_ln_last):
                    os.unlink(obs_ln_last)
                os.symlink(os.path.relpath(obs_dir, chunk_dir_path),
                           obs_ln_last)
        if ('default' in obsdirs
            and (os.path.realpath(obsdirs['last'])
                 == os.path.realpath(obsdirs['default']))):
            if os.path.exists(obs_ln_last):
                os.unlink(obs_ln_last)
    obi.isdefault = isdefault

def update(obsids=[]):
    asp_l1_proc = Ska.DBI.DBI(dbi="sqlite", server=FILES['asp1_proc_table'])
    # if an obsid is requested, just do that
    # there's a little bit of duplication in this block
    if len(obsids):
        for obsid in obsids:
            todo = asp_l1_proc.fetchall(
                "SELECT * FROM aspect_1_proc where obsid = {}".format(
                    obsid))
            for obs in todo:
                if obs['obsid'] in KNOWN_BAD_OBSIDS:
                    logger.info("Skipping known bad obsid {}".format(obs['obsid']))
                    continue
                logger.info("running VV for obsid {} run on {}".format(
                    obs['obsid'], obs['ap_date']))
                process(obs['obsid'], version=obs['revision'])
                asp_l1_proc.execute("""UPDATE aspect_1_proc set vv_complete = 1
                                       where obsid = {} and revision = {}
                                    """.format(obs['obsid'], obs['revision']))
    else:
        # if no obsid specified, try to retrieve all data without vv
        todo = asp_l1_proc.fetchall(
            """SELECT * FROM aspect_1_proc
               where vv_complete != 1 order by aspect_1_id""")
        for obs in todo:
            if obs['obsid'] in KNOWN_BAD_OBSIDS:
                logger.info("Skipping known bad obsid {}".format(obs['obsid']))
                continue
            logger.info("running VV for obsid {} run on {}".format(
                obs['obsid'], obs['ap_date']))
            process(obs['obsid'], version=obs['revision'])
            asp_l1_proc.execute("""UPDATE aspect_1_proc set vv_complete = 1
                                   where obsid = {} and revision = {}
                                """.format(obs['obsid'], obs['revision']))


def get_arch_vv(obsid, version='last'):
    """
    Get obsid paths from mica.archive and create mica.vv.Obi object
    """
    asp_l1_dirs = asp_l1_arch.get_obs_dirs(obsid)
    if version not in asp_l1_dirs:
        logger.error(
            "Requested version {} not in asp_l1 archive".format(version))
        return None
    l1_dir = asp_l1_dirs[version]
    # find the obspar that matches the requested aspect_1 products
    # this is in the aspect processing table
    asp_l1_proc = Ska.DBI.DBI(dbi="sqlite", server=FILES['asp1_proc_table'])
    asp_obs = asp_l1_proc.fetchall(
        "SELECT * FROM aspect_1_proc where obsid = {}".format(
            obsid))
    asp_proc = None
    if version == 'last':
        asp_proc = asp_obs[asp_obs['aspect_1_id']
                           == np.max(asp_obs['aspect_1_id'])][0]
    if version == 'default':
        asp_proc = asp_obs[asp_obs['isdefault'] == 1][0]
    if asp_proc is None:
        asp_proc = asp_obs[asp_obs['revision'] == version][0]
    obspar_dirs = obspar_arch.get_obs_dirs(obsid)
    if asp_proc['obspar_version'] not in obspar_dirs:
        # try to update the obspar archive with the missing version
        config = obspar_arch.CONFIG.copy()
        config.update(dict(obsid=obsid, version=asp_proc['obspar_version']))
        oa = obsid_archive.ObsArchive(config)
        oa.logger.setLevel(logging.INFO)
        oa.logger.addHandler(logging.StreamHandler())
        oa.update()
        obspar_dirs = obspar_arch.get_obs_dirs(obsid)
    obspar_file = glob(os.path.join(obspar_dirs[asp_proc['obspar_version']],
                                    'axaf*par*'))[0]
    return Obi(obspar_file, l1_dir, temproot=FILES['temp_root'])


def process(obsid, version='last'):
    obi = get_arch_vv(obsid, version)
    obi.save_plots_and_resid()
    if obi is None:
        return None
    shelf_file = os.path.join(FILES['data_root'],
                              FILES['shelf_file'])
    obi.shelve_info(shelf_file)
    file_vv(obi)
    h5 = tables.openFile(FILES['h5_file'], 'a')
    tbl = h5.getNode('/', 'vv')
    obi.set_tbl(tbl)
    obi.slots_to_table()
    tbl.flush()
    h5.flush()
    h5.close()
    return obi


