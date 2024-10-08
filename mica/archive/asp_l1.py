#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Update aspect archive L1 products.

Script to update Ska file archive aspect L1 products.  Module
also provides methods to retrieve the directory (or directories)
for an obsid.

This uses the obsid_archive module with a configuration specific
to the aspect L1 products.

"""

import logging
import os

import numpy as np
from astropy.table import Table
from Quaternion import Quat

from mica.archive import asp_l1_proc, obsid_archive  # noqa: F401 imported but unused
from mica.common import MICA_ARCHIVE

# these columns are available in the headers of the fetched telemetry
# for this product (ASP L1) and will be included in the file lookup table
ARCHFILES_HDR_COLS = (
    "tstart",
    "tstop",
    "caldbver",
    "content",
    "ascdsver",
    "revision",
    "date",
)

# config = ConfigObj('asp1.conf')
CONFIG = dict(
    data_root=os.path.join(MICA_ARCHIVE, "asp1"),
    temp_root=os.path.join(MICA_ARCHIVE, "temp"),
    bad_obsids=os.path.join(MICA_ARCHIVE, "asp1", "asp_l1_bad_obsids.dat"),
    sql_def="archfiles_asp_l1_def.sql",
    apstat_table="aspect_1",
    apstat_id="aspect_1_id",
    label="asp_l1",
    small="asp1{fidprops}",
    small_glob="*fidpr*",
    small_ver_regex=r"pcadf\d+N(\d{3})_",
    delete=["kalm1", "gdat", "adat"],
    full="asp1",
    filecheck=False,
    cols=ARCHFILES_HDR_COLS,
    content_types=[
        "ASPQUAL",
        "ASPSOL",
        "ACADATA",
        "GSPROPS",
        "GYRODATA",
        "KALMAN",
        "ACACAL",
        "ACACENT",
        "FIDPROPS",
        "GYROCAL",
        "ACA_BADPIX",
    ],
)


def get_options():
    import argparse

    desc = """
Run the update process to get new ASP L1 telemetry, save it in the Ska
file archive, and include it in the file lookup database.  This is intended
to be run as a cron task, and in regular processing, the update will fetch
and ingest all telemetry since the task's last run.  Options also provided
to fetch and ingest specific obsids and versions.

See the ``CONFIG`` in the asp_l1.py file and the config description in
obsid_archive for more information on the asp l1 default config if parameters
without command-line options need to be changed.
"""
    parser = argparse.ArgumentParser(description=desc)
    defaults = dict(CONFIG)
    parser.set_defaults(**defaults)
    parser.add_argument("--obsid", type=int, help="specific obsid to process")
    parser.add_argument(
        "--version", default="last", help="specific processing version to retrieve"
    )
    parser.add_argument(
        "--firstrun",
        action="store_true",
        help="for archive init., ignore rev in aspect_1 table",
    )
    parser.add_argument("--data-root", help="parent directory for all data")
    parser.add_argument("--temp-root", help="parent temp directory")
    parser.add_argument(
        "--filecheck",
        action="store_true",
        help="for provisional data, download files and check"
        " that all are present.  If unset, proceed if dir"
        " exists",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Allow update to rebuild archive from obsid 1",
    )
    opt = parser.parse_args()
    return opt


# set up an archive object with default config for use by the other
# get_* methods
archive = obsid_archive.ObsArchive(CONFIG)


def get_dir(obsid):
    """
    Get ASP L1 directory for default/released products for an obsid.

      >>> from mica.archive import asp_l1
      >>> asp_l1.get_dir(2121)
      '/proj/sot/ska/data/mica/archive/asp1/02/02121'

    :param obsid: obsid
    :returns: directory
    :rtype: string
    """
    return archive.get_dir(obsid)


def get_obs_dirs(obsid):
    """
    Get all ASP L1 directories for an obsid in the Ska file archive.

      >>> from mica.archive import asp_l1
      >>> obsdirs = asp_l1.get_obs_dirs(6000)

    obsdirs will look something like::

      {'default': '/proj/sot/ska/data/mica/archive/asp1/06/06000',
      2: '/proj/sot/ska/data/mica/archive/asp1/06/06000_v02',
      3: '/proj/sot/ska/data/mica/archive/asp1/06/06000_v03',
      'last': '/proj/sot/ska/data/mica/archive/asp1/06/06000',
      'revisions': [2, 3]}

    :param obsid: obsid
    :returns: map of obsid version to directories
    :rtype: dictionary
    """
    return archive.get_obs_dirs(obsid)


def get_files(obsid=None, start=None, stop=None, revision=None, content=None):
    """
    List asp_l1 files for an obsid or a time range.

      >>> from mica.archive import asp_l1
      >>> obs_files = asp_l1.get_files(6000)
      >>> obs_gspr = asp_l1.get_files(6000, content=['GSPROPS'])
      >>> range_fidpr = asp_l1.get_files(start='2012:001',
      ...                                stop='2012:030',
      ...                                content=['FIDPROPS'])


    The available content types are: ASPQUAL, ASPSOL, ASPSOLOBI, ACACAL,
    ACA_BADPIX, FIDPROPS, GYROCAL, GSPROPS, and ACACENT.

    :param obsid: obsid
    :param start: time range start (Chandra.Time compatible)
    :param stop: time range stop (Chandra.Time compatible)
    :param revision: revision integer or 'last'
                     defaults to current released version
    :param content: archive CONTENT type
                    defaults to all available ASP1 types
    :returns: full path of files matching query
    """
    return archive.get_files(
        obsid=obsid, start=start, stop=stop, revision=revision, content=content
    )


def get_atts(obsid=None, start=None, stop=None, revision=None, filter=True):
    """
    Get ground aspect attitudes

    Get the ground aspect solution quaternions and times covering obsid or start to stop,
    in the ACA frame.

    :obsid: obsid
    :start: start time (DateTime compat)
    :stop: stop time (DateTime compat)
    :revision: aspect pipeline processing revision (integer version, None, or 'last')
    :filter: boolean, true means returned values will not include quaternions during times
                      when asp_sol_status is non-zero

    :returns: Nx4 np.array of quaternions, np.array of N times, list of dict with header
              from each asol file.
    """
    if revision == "all":
        raise ValueError("revision 'all' doesn't really make sense for this function")
    # These are in time order by default from get_files
    asol_files = get_files(
        obsid=obsid, start=start, stop=stop, revision=revision, content=["ASPSOL"]
    )
    acal_files = get_files(
        obsid=obsid, start=start, stop=stop, revision=revision, content=["ACACAL"]
    )
    aqual_files = get_files(
        obsid=obsid, start=start, stop=stop, revision=revision, content=["ASPQUAL"]
    )
    return get_atts_from_files(asol_files, acal_files, aqual_files, filter=filter)


def get_atts_from_files(asol_files, acal_files, aqual_files, filter=True):
    """
    Get ground aspect attitudes from ASP1 source files.

    From ASP1 source files (asol, acal, aqual) get the ground aspect solution quaternions
    and times covering the range of asol_files in the ACA frame.  The asol, acl, and aqual files
    are assumed to have one-to-one correspondence (though the asol to acal times are checked).

    :asol_files: list of aspect asol1 files
    :acal_files: list of acal1 files associated with asol_files
    :aqual_files: list of aqual files associated with asol_files
    :filter: boolean, true means returned values will not include quaternions during times when
                      asp_sol_status is non-zero

    :returns: Nx4 np.array of quaternions, np.array of N times, list of dict with header
              from each asol file.
    """
    # There should be one asol and one acal file for each aspect interval in the range
    att_chunks = []
    time_chunks = []
    records = []
    for asol_f, acal_f, aqual_f in zip(asol_files, acal_files, aqual_files):
        asol = Table.read(asol_f, hdu=1)
        acal = Table.read(acal_f, hdu=1)
        aqual = Table.read(aqual_f, hdu=1)
        # Check that the time ranges match from the fits headers (meta in the table)
        if not np.allclose(
            np.array([asol.meta["TSTART"], asol.meta["TSTOP"]]),
            np.array([acal.meta["TSTART"], acal.meta["TSTOP"]]),
            atol=10,
        ):
            raise ValueError("ACAL and ASOL have mismatched time ranges")
        if filter and np.any(aqual["asp_sol_status"] != 0):
            # For each sample with bad status find the overlapping time range in asol
            # and remove from set
            for idx in np.flatnonzero(aqual["asp_sol_status"]):
                nok = (asol["time"] >= (aqual["time"][idx] - 1.025)) & (
                    asol["time"] <= (aqual["time"][idx] + 1.025)
                )
                asol = asol[~nok]

        # Transpose of transform/rotation matrix is the inverse
        aca_mis_inv = acal["aca_misalign"][0].transpose()

        # Repeat the transform to match asol
        aca_mis_inv = np.repeat(aca_mis_inv[np.newaxis, ...], repeats=len(asol), axis=0)
        q_mis_inv = Quat(transform=aca_mis_inv)

        # Quaternion multiply the asol quats with that inv misalign and save
        q_att_name = "q_att_raw" if "q_att_raw" in asol.colnames else "q_att"
        att_chunks.append((Quat(q=asol[q_att_name]) * q_mis_inv).q)
        time_chunks.append(np.array(asol["time"]))
        records.append(asol.meta)
    if len(att_chunks) > 0 and len(time_chunks) > 0:
        return np.vstack(att_chunks), np.hstack(time_chunks), records
    else:
        return np.array([]), np.array([]), []


def main():
    """
    Run the ASP L1 update process.

    Run the update process to get new ASP L1 telemetry, save it in the Ska
    file archive, and include it in the file lookup database.
    """
    opt = get_options()
    config = vars(opt)
    archive = obsid_archive.ObsArchive(config)
    archive.logger.setLevel(logging.INFO)
    archive.logger.addHandler(logging.StreamHandler())
    obsids = archive.update()  # noqa: F841 assigned but not used


if __name__ == "__main__":
    main()
