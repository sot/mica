# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division

import json
import logging
import os
from glob import glob

import numpy as np
import ska_dbi
import tables

from mica.common import MICA_ARCHIVE

FILES = dict(
    data_root=os.path.join(MICA_ARCHIVE, "vv"),
    temp_root=os.path.join(MICA_ARCHIVE, "tempvv"),
    shelf_file=os.path.join(MICA_ARCHIVE, "vv", "vv_shelf"),
    h5_file=os.path.join(MICA_ARCHIVE, "vv", "vv.h5"),
    last_file=os.path.join(MICA_ARCHIVE, "vv", "last_id.txt"),
    asp1_proc_table=os.path.join(MICA_ARCHIVE, "asp1", "processing_asp_l1.db3"),
    bad_obsid_list=os.path.join(MICA_ARCHIVE, "vv", "bad_obsids.json"),
)


logger = logging.getLogger("vv")
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def get_vv_dir(obsid, version="default"):
    """
    Get directory containing V&V products.

    Get directory containing V&V products for a requested obsid/version,
    including plots and json.

    :param obsid: obsid
    :param version: 'last', 'default' or version number
    :returns: directory name for obsid/version
    """
    num_version = None
    if version in ["last", "default"]:
        asp_l1_proc = ska_dbi.DBI(dbi="sqlite", server=FILES["asp1_proc_table"])
        if version == "default":
            obs = asp_l1_proc.fetchall(
                """select * from aspect_1_proc
                                          where obsid = {} and isdefault = 1
                                       """.format(obsid)
            )
            if not len(obs):
                raise LookupError(
                    "Version {} not found for obsid {}".format(version, obsid)
                )
            num_version = obs["revision"][0]
        if version == "last":
            obs = asp_l1_proc.fetchall(
                """select * from aspect_1_proc
                                          where obsid = {}
                                       """.format(obsid)
            )
            if not len(obs):
                raise LookupError("No entries found for obsid {}".format(obsid))
            num_version = np.max(obs["revision"])
    else:
        num_version = version
    strobs = "%05d_v%02d" % (obsid, num_version)
    chunk_dir = strobs[0:2]
    chunk_dir_path = os.path.join(FILES["data_root"], chunk_dir)
    obs_dir = os.path.join(chunk_dir_path, strobs)
    if not os.path.exists(obs_dir):
        raise LookupError("Expected vv archive dir {} not found".format(obs_dir))
    return obs_dir


def get_vv_files(obsid, version="default"):
    """
    Get list of V&V files available for a requested obsid/version.

    :param obsid: obsid
    :param version: 'default', 'last' or version number
    :returns: list of files
    """
    vv_dir = get_vv_dir(obsid, version)
    return glob(os.path.join(vv_dir, "*"))


def get_vv(obsid, version="default"):
    """
    Retrieve V&V data for an obsid/version.

    This reads the saved JSON and returns the previously-
    calculated V&V data.

    :param obsid: obsid
    :param version: 'last', 'default', or version number
    :returns: dict of V&V data
    """
    vv_dir = get_vv_dir(obsid, version)
    json_file = glob(os.path.join(vv_dir, "*.json"))[0]
    return json.loads(open(json_file).read())


def get_rms_data():
    """
    Retrieve/return all data from RMS trending H5 archive

    :returns: numpy array of RMS data for each star/obsid/version
    """
    tables_open_file = getattr(tables, "open_file", None) or tables.openFile
    with tables_open_file(FILES["h5_file"], "r") as h5f:
        data = h5f.root.vv[:]
    return data
