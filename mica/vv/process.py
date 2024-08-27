# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division

import json
import logging
import os
import shutil
import tempfile
from glob import glob

import numpy as np
import ska_dbi
import tables

import mica.archive.asp_l1 as asp_l1_arch
import mica.archive.obspar as obspar_arch
from mica.archive import obsid_archive

from .core import VV_VERSION, Obi
from .vv import FILES

VV_DTYPE = np.dtype(
    [
        ("obsid", "<i4"),
        ("revision", "<i4"),
        ("isdefault", "<i4"),
        ("aspect_1_id", "<i4"),
        ("used", "<i4"),
        ("vv_version", "<i4"),
        ("ap_date", "|S21"),
        ("tstart", "<f8"),
        ("tstop", "<f8"),
        ("sim_z", "<f8"),
        ("sim_z_offset", "<f8"),
        ("instrument", "|S10"),
        ("ra_pnt", "<f8"),
        ("dec_pnt", "<f8"),
        ("roll_pnt", "<f8"),
        ("slot", "<i4"),
        ("type", "|S10"),
        ("n_pts", "<i4"),
        ("rad_off", "<f8"),
        ("frac_dy_big", "<f8"),
        ("frac_dz_big", "<f8"),
        ("frac_mag_big", "<f8"),
        ("mean_y", "<f8"),
        ("mean_z", "<f8"),
        ("dy_mean", "<f8"),
        ("dy_med", "<f8"),
        ("dy_rms", "<f8"),
        ("dz_mean", "<f8"),
        ("dz_med", "<f8"),
        ("dz_rms", "<f8"),
        ("dr_mean", "<f8"),
        ("dr_med", "<f8"),
        ("dr_rms", "<f8"),
        ("mag_mean", "<f8"),
        ("mag_med", "<f8"),
        ("mag_rms", "<f8"),
        ("mean_aacccdpt", "<f8"),
    ]
)


KNOWN_BAD_OBSIDS = []
if os.path.exists(FILES["bad_obsid_list"]):
    KNOWN_BAD_OBSIDS = json.loads(open(FILES["bad_obsid_list"]).read())

logger = logging.getLogger("vv")


def _file_vv(obi):
    """
    Save processed V&V data to per-obsid archive
    """
    obsid = int(obi.info()["obsid"])
    version = int(obi.info()["revision"])
    # set up directory for data
    strobs = "%05d_v%02d" % (obsid, version)
    chunk_dir = strobs[0:2]
    chunk_dir_path = os.path.join(FILES["data_root"], chunk_dir)
    obs_dir = os.path.join(chunk_dir_path, strobs)
    if not os.path.exists(obs_dir):
        logger.info("making directory %s" % obs_dir)
        os.makedirs(obs_dir)
    else:
        logger.info("obsid dir %s already exists" % obs_dir)
    for f in glob(os.path.join(obi.tempdir, "*")):
        os.chmod(f, 0o775)
        shutil.copy(f, obs_dir)
        os.remove(f)
    logger.info("moved VV files to {}".format(obs_dir))
    os.removedirs(obi.tempdir)
    logger.info("removed directory {}".format(obi.tempdir))
    # make any desired link
    obs_ln = os.path.join(FILES["data_root"], chunk_dir, "%05d" % obsid)
    obs_ln_last = os.path.join(FILES["data_root"], chunk_dir, "%05d_last" % obsid)
    obsdirs = asp_l1_arch.get_obs_dirs(obsid)
    isdefault = 0
    if "default" in obsdirs:
        if os.path.realpath(obsdirs[version]) == os.path.realpath(obsdirs["default"]):
            if os.path.islink(obs_ln):
                os.unlink(obs_ln)
            os.symlink(os.path.relpath(obs_dir, chunk_dir_path), obs_ln)
            isdefault = 1
    if "last" in obsdirs:
        if (
            "default" in obsdirs
            and (
                os.path.realpath(obsdirs["last"])
                != os.path.realpath(obsdirs["default"])
            )
            or "default" not in obsdirs
        ):
            if os.path.realpath(obsdirs[version]) == os.path.realpath(obsdirs["last"]):
                if os.path.islink(obs_ln_last):
                    os.unlink(obs_ln_last)
                os.symlink(os.path.relpath(obs_dir, chunk_dir_path), obs_ln_last)
        if "default" in obsdirs and (
            os.path.realpath(obsdirs["last"]) == os.path.realpath(obsdirs["default"])
        ):
            if os.path.exists(obs_ln_last):
                os.unlink(obs_ln_last)
    obi.isdefault = isdefault


def update(obsids=None):
    """
    Run V&V process over obsids.

    For a list of obsids or for all new obsids, run V&V processing
    and save V&V info to archive.

    :param obsids: optional list of obsids
    """
    if obsids is None:
        # If no obsid specified, run on all with vv_complete = 0 in
        # the local processing database.
        with ska_dbi.DBI(dbi="sqlite", server=FILES["asp1_proc_table"]) as db:
            obsids = db.fetchall(
                """SELECT * FROM aspect_1_proc
                where vv_complete = 0
                order by aspect_1_id"""
            )["obsid"]
    for obsid in obsids:
        with ska_dbi.DBI(dbi="sqlite", server=FILES["asp1_proc_table"]) as db:
            proc = db.fetchall(
                f"SELECT obsid, revision, ap_date FROM aspect_1_proc where obsid = {obsid}"
            )
        for obs in proc:
            if obsid in KNOWN_BAD_OBSIDS:
                logger.info(f"Skipping known bad obsid {obsid}")
                continue
            logger.info(f"running VV for obsid {obsid} run on {obs['ap_date']}")
            try:
                process(obsid, version=obs["revision"])
            except LookupError:
                logger.warn(f"Skipping obs:ver {obsid}:{obs['revision']}. Missing data")
                continue
            update_str = f"""UPDATE aspect_1_proc set vv_complete = {VV_VERSION}
                              where obsid = {obsid} and revision = {obs['revision']}"""

            logger.info(update_str)
            with ska_dbi.DBI(dbi="sqlite", server=FILES["asp1_proc_table"]) as db:
                db.execute(update_str)


def get_arch_vv(obsid, version="last"):
    """
    Retrieve and load archived V&V.

    Given obsid and version, find archived ASP1 and obspar products and
    run V&V.  Effort is made to find the obspar that was actually used during
    creation of the ASP1 products.

    :param obsid: obsid
    :param version: 'last', 'default', or revision number of ASP1 products
    :returns: mica.vv.Obi V&V object
    """
    logger.info("Generating V&V for obsid {}".format(obsid))
    asp_l1_dirs = asp_l1_arch.get_obs_dirs(obsid)
    if asp_l1_dirs is None or version not in asp_l1_dirs:
        raise LookupError("Requested version {} not in asp_l1 archive".format(version))
    l1_dir = asp_l1_dirs[version]
    # find the obspar that matches the requested aspect_1 products
    # this is in the aspect processing table
    asp_l1_proc = ska_dbi.DBI(dbi="sqlite", server=FILES["asp1_proc_table"])
    asp_obs = asp_l1_proc.fetchall(
        "SELECT * FROM aspect_1_proc where obsid = {}".format(obsid)
    )
    asp_proc = None
    if len(asp_obs) == 0:
        return None
    if version == "last":
        asp_proc = asp_obs[asp_obs["aspect_1_id"] == np.max(asp_obs["aspect_1_id"])][0]
    if version == "default":
        asp_proc = asp_obs[asp_obs["isdefault"] == 1][0]
    if asp_proc is None:
        asp_proc = asp_obs[asp_obs["revision"] == version][0]
    obspar_dirs = obspar_arch.get_obs_dirs(obsid)
    if obspar_dirs is None or asp_proc["obspar_version"] not in obspar_dirs:
        # try to update the obspar archive with the missing version
        config = obspar_arch.CONFIG.copy()
        config.update(dict(obsid=obsid, version=asp_proc["obspar_version"]))
        oa = obsid_archive.ObsArchive(config)
        oa.logger.setLevel(logging.INFO)
        oa.logger.addHandler(logging.StreamHandler())
        oa.update()
        obspar_dirs = obspar_arch.get_obs_dirs(obsid)
    try:
        obspar_file = glob(
            os.path.join(obspar_dirs[asp_proc["obspar_version"]], "axaf*par*")
        )[0]
    except IndexError as err:
        raise LookupError(f"Requested version {version} not in obspar archive") from err
    return Obi(obspar_file, l1_dir, temproot=FILES["temp_root"])


def process(obsid, version="last"):
    """
    Run V&V process for an obsid.

    For requested obsid/version, run V&V, make plots,
    save plots and JSON, save info to shelve file, and
    update RMS HDF5.

    :param obsid: obsid
    :param version: 'last', 'default' or revision number of ASP1 products
    :returns: mica.vv.Obi V&V object
    """
    obi = get_arch_vv(obsid, version)
    if not os.path.exists(FILES["temp_root"]):
        os.makedirs(FILES["temp_root"])
    if obi is None:
        return None
    obi.tempdir = tempfile.mkdtemp(dir=FILES["temp_root"])
    obi.save_plots_and_resid()
    _file_vv(obi)
    if not os.path.exists(FILES["h5_file"]):
        vv_desc, byteorder = tables.descr_from_dtype(VV_DTYPE)
        filters = tables.Filters(complevel=5, complib="zlib")
        h5f = tables.open_file(FILES["h5_file"], "a")
        tbl = h5f.create_table("/", "vv", vv_desc, filters=filters, expectedrows=1e6)
        h5f.close()
    h5 = tables.open_file(FILES["h5_file"], "a")
    tbl = h5.get_node("/", "vv")
    obi.set_tbl(tbl)
    obi.slots_to_table()
    tbl.flush()
    h5.flush()
    h5.close()
    return obi
