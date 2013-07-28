from __future__ import division
import os
import json
import tables
import shutil
import logging
from glob import glob
import numpy as np

import Ska.DBI

import mica.archive.asp_l1 as asp_l1_arch
import mica.archive.obspar as obspar_arch
from mica.archive import obsid_archive
from .core import Obi
from .vv import FILES

KNOWN_BAD_OBSIDS = json.loads(open(FILES['bad_obsid_list']).read())

logger = logging.getLogger('vv')


def _file_vv(obi):
    """
    Save processed V&V data to per-obsid archive
    """
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
    """
    For a list of obsids or for all new obsids, run V&V processing
    and save V&V info to archive.

    :param obsids: optional list of obsids
    """
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
        # set to only use aspect_1_id greater than those currently
        # in the archive at time of vv install.
        # this is to avoid processing old no-longer-default data for now
        todo = asp_l1_proc.fetchall(
            """SELECT * FROM aspect_1_proc
               where vv_complete != 1 and aspect_1_id > 34878
               order by aspect_1_id""")
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
    Given obsid and version, find archived ASP1 and obspar products and
    run V&V.  Effort is made to find the obspar that was actually used during
    creation of the ASP1 products.

    :param obsid: obsid
    :param version: 'last', 'default', or revision number of ASP1 products
    :returns: mica.vv.Obi V&V object
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
    """
    For requested obsid/version, run V&V, make plots,
    save plots and JSON, save info to shelve file, and
    update RMS HDF5.

    :param obsid: obsid
    :param version: 'last', 'default' or revision number of ASP1 products
    :returns: mica.vv.Obi V&V object
    """
    obi = get_arch_vv(obsid, version)
    obi.save_plots_and_resid()
    if obi is None:
        return None
    shelf_file = os.path.join(FILES['data_root'],
                              FILES['shelf_file'])
    obi.shelve_info(shelf_file)
    _file_vv(obi)
    h5 = tables.openFile(FILES['h5_file'], 'a')
    tbl = h5.getNode('/', 'vv')
    obi.set_tbl(tbl)
    obi.slots_to_table()
    tbl.flush()
    h5.flush()
    h5.close()
    return obi
