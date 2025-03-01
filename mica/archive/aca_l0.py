# Licensed under a 3-clause BSD style license - see LICENSE.rst
import argparse
import collections
import gzip
import logging
import os
import re
import shutil
import time
from glob import glob
from itertools import count
from pathlib import Path

import astropy.io.fits as pyfits
import numpy as np
import Ska.arc5gl
import Ska.File
import ska_dbi
import tables
from astropy.table import Table, vstack
from Chandra.Time import DateTime
from chandra_aca.aca_image import ACAImage
from cxotime import CxoTimeLike
from numpy import ma

from mica.common import MICA_ARCHIVE, MissingDataError

# import kadi later in obsid_times


logger = logging.getLogger("aca0 fetch")
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

# borrowed from eng_archive
ARCHFILES_HDR_COLS = (
    "tstart",
    "tstop",
    "startmjf",
    "startmnf",
    "stopmjf",
    "stopmnf",
    "tlmver",
    "ascdsver",
    "revision",
    "date",
    "imgsize",
)

FILETYPE = {
    "level": "L0",
    "instrum": "PCAD",
    "content": "ACADATA",
    "arc5gl_query": "ACA0",
    "fileglob": "aca*fits*",
}

ACA_DTYPE = (
    ("TIME", ">f8"),
    ("QUALITY", ">i4"),
    ("MJF", ">i4"),
    ("MNF", ">i4"),
    ("END_INTEG_TIME", ">f8"),
    ("INTEG", ">f4"),
    ("PIXTLM", "<U4"),
    ("BGDTYP", "<U4"),
    ("GLBSTAT", "|u1"),
    ("COMMCNT", "|u1"),
    ("COMMPROG", "|u1"),
    ("IMGFID1", "|u1"),
    ("IMGNUM1", "|u1"),
    ("IMGFUNC1", "|u1"),
    ("IMGSTAT", "|u1"),
    ("IMGROW0", ">i2"),
    ("IMGCOL0", ">i2"),
    ("IMGSCALE", ">i2"),
    ("BGDAVG", ">i2"),
    ("IMGFID2", "|u1"),
    ("IMGNUM2", "|u1"),
    ("IMGFUNC2", "|u1"),
    ("BGDRMS", ">i2"),
    ("TEMPCCD", ">f4"),
    ("TEMPHOUS", ">f4"),
    ("TEMPPRIM", ">f4"),
    ("TEMPSEC", ">f4"),
    ("BGDSTAT", "|u1"),
    ("IMGFID3", "|u1"),
    ("IMGNUM3", "|u1"),
    ("IMGFUNC3", "|u1"),
    ("IMGFID4", "|u1"),
    ("IMGNUM4", "|u1"),
    ("IMGFUNC4", "|u1"),
    ("IMGRAW", ">f4", (64,)),
    ("HD3TLM62", "|u1"),
    ("HD3TLM63", "|u1"),
    ("HD3TLM64", "|u1"),
    ("HD3TLM65", "|u1"),
    ("HD3TLM66", "|u1"),
    ("HD3TLM67", "|u1"),
    ("HD3TLM72", "|u1"),
    ("HD3TLM73", "|u1"),
    ("HD3TLM74", "|u1"),
    ("HD3TLM75", "|u1"),
    ("HD3TLM76", "|u1"),
    ("HD3TLM77", "|u1"),
    ("IMGSIZE", ">i4"),
    ("FILENAME", "<U128"),
)

ACA_DTYPE_8x8 = tuple(
    (dt[0], dt[1], (8, 8)) if dt[0] == "IMGRAW" else dt for dt in ACA_DTYPE
)

ACA_DTYPE_NAMES = tuple([k[0] for k in ACA_DTYPE])

CONFIG = dict(
    data_root=os.path.join(MICA_ARCHIVE, "aca0"),
    temp_root=os.path.join(MICA_ARCHIVE, "temp"),
    days_at_once=30.0,
    sql_def="archfiles_aca_l0_def.sql",
    cda_table="cda_aca0.h5",
)

# Names of bits in GLBSTAT and IMGSTAT fields, in order from most to least significant
GLBSTAT_NAMES = [
    "HIGH_BGD",
    "RAM_FAIL",
    "ROM_FAIL",
    "POWER_FAIL",
    "CAL_FAIL",
    "COMM_CHECKSUM_FAIL",
    "RESET",
    "SYNTAX_ERROR",
]

IMGSTAT_NAMES = [
    "SAT_PIXEL",
    "DEF_PIXEL",
    "QUAD_BOUND",
    "COMMON_COL",
    "MULTI_STAR",
    "ION_RAD",
]


def get_options():
    parser = argparse.ArgumentParser(
        description="Fetch aca level 0 products and make a file archive"
    )
    defaults = dict(CONFIG)
    parser.set_defaults(**defaults)
    parser.add_argument("--data-root", help="parent directory for all data")
    parser.add_argument("--temp-root", help="parent temp directory")
    parser.add_argument(
        "--start",
        help="start date for retrieve " + "(defaults to max date of archived files)",
    )
    parser.add_argument("--stop", help="stop date for retrieve " + "(defaults to now)")
    parser.add_argument(
        "--days-at-once",
        type=float,
        help="if over this number, " + "bin to chunks of this number of days",
    )
    parser.add_argument(
        "--cda-table", help="file name for h5 list from Chandra Archive"
    )
    opt = parser.parse_args()
    return opt


def _extract_bits_as_uint8(
    arr: np.ndarray, bit_start: int, bit_stop: int
) -> np.ndarray:
    """Extract UINT8's from a 2-D array of uint8 bit values into a new array.

    Input array is Nx8 where N is the number of rows and 8 is the number of bits.
    It must be in little-endian order with the most significant bit in the rightmost
    position of the array.
    """
    # This returns a shape of (N, 1) so we need to squeeze it.
    out = np.packbits(arr[:, bit_start:bit_stop], axis=1, bitorder="little")
    return out.squeeze()


def get_aca_images(
    start: CxoTimeLike,
    stop: CxoTimeLike,
):
    """
    Get ACA images from mica L0 archive in same format as maude_decom.get_aca_images.

    This gets a table of images from the mica L0 file archive that includes the same
    output columns as ``chandra_aca.maude_decom.get_aca_images()``. The images are 8x8
    and are centered in the 8x8 image, with masks indicating empty pixels (for 4x4 or
    6x6 data).

    The units of the image data are DN, so you need to multiply by 5.0 / INTEG to get
    e-/sec.

    This function returns additional columns that are not in the maude images:
    - 'HD3TLM*': uint8, raw 8-bit header-3 telemetry values
    - 'IMGSIZE': int32, 4, 6, or 8

    Parameters
    ----------
    start : CxoTimeLike
        Start time of the time range.
    stop : CxoTimeLike
        Stop time of the time range.

    Returns
    -------
    astropy.table.Table
        ACA images and accompanying data.
    """
    slot_data_list = []
    for slot in range(8):
        slot_data_raw = get_slot_data(
            start,
            stop,
            slot=slot,
            centered_8x8=True,
        )
        # Convert from numpy structured array to astropy table, keeping only good data
        ok = slot_data_raw["QUALITY"] == 0
        slot_data = Table(slot_data_raw[ok])

        # Unmask everything that can be unmasked
        # Remove mask from columns where no values are masked. IMGRAW is excepted because
        # the mask reflects the presence of 4x4 or 6x6 images, not entirely missing data
        # per row.
        for col in slot_data.itercols():
            if not np.any(col.mask) and col.name != "IMGRAW":
                slot_data[col.name] = col.data.data

        # Add slot number if there are any rows (if statement not needed after
        # https://github.com/astropy/astropy/pull/17102).
        # if len(slot_data) > 0:
        #     slot_data["IMGNUM"] = slot

        for name in ["IMGFID", "IMGFUNC", "IMGNUM"]:
            slot_data[name] = slot_data[f"{name}1"]
            for ii in (1, 2, 3, 4):
                name_ii = f"{name}{ii}"
                if name_ii in slot_data.colnames:
                    del slot_data[name_ii]

        slot_data_list.append(slot_data)

    # Combine all slot data into one output table and sort by time and slot
    out = vstack(slot_data_list)
    out.sort(keys=["TIME", "IMGNUM"])

    # Rename and rejigger colums to match chandra_aca.maude_decum.get_aca_images
    is_6x6 = out["IMGSIZE"] == 6
    is_4x4 = out["IMGSIZE"] == 4
    for rc in ["ROW", "COL"]:
        # IMGROW/COL0 are row/col of lower/left of 4x4, 6x6, or 8x8 image
        # Make these new columns:
        # IMGROW/COL_A1: row/col of pixel A1 of 4x4, 6x6, or 8x8 image
        # IMGROW/COL0_8x8: row/col of lower/left of 8x8 image with smaller img centered
        out[f"IMG{rc}_A1"] = out[f"IMG{rc}0"]
        out[f"IMG{rc}_A1"][is_6x6] += 1
        out[f"IMG{rc}0_8X8"] = out[f"IMG{rc}0"]
        out[f"IMG{rc}0_8X8"][is_6x6] -= 1
        out[f"IMG{rc}0_8X8"][is_4x4] -= 2

    out["IMGTYPE"] = np.full(shape=len(out), fill_value=4)  # 8x8
    out["IMGTYPE"][is_6x6] = 1  # 6x6
    out["IMGTYPE"][is_4x4] = 0  # 4x4

    out.rename_column("IMGRAW", "IMG")
    out.rename_column("BGDTYP", "AABGDTYP")
    out.rename_column("PIXTLM", "AAPIXTLM")
    del out["FILENAME"]
    del out["QUALITY"]

    # maude_decom.get_aca_images() provides both VCDUCTR (VCDU for each sub-image) and
    # IMG_VCDUCTR (VCDU of first sub-image). For combined images, these are identical.
    out["VCDUCTR"] = out["MJF"] * 128 + out["MNF"]
    out["IMG_VCDUCTR"] = out["VCDUCTR"]

    # Split uint8 GLBSTAT and IMGSTAT into individual bits. The stat_bit_names are
    # in MSB order so we need to reverse them below.
    for stat_name, stat_bit_names in (
        ("GLBSTAT", GLBSTAT_NAMES),
        ("IMGSTAT", IMGSTAT_NAMES),
    ):
        for bit, name in zip(range(8), reversed(stat_bit_names)):
            out[name] = (out[stat_name] & (1 << bit)) > 0

    # Extract fields from the mica L0 8-bit COMMCNT value which is really three MSIDs.
    # We need to use .data attribute of MaskedColumn to get a numpy masked array. The
    # unpackbits function fails on MaskedColumn because it looks for a _mask attribute.
    bits = np.unpackbits(out["COMMCNT"].data).reshape(-1, 8)
    bits_rev = bits[:, ::-1]  # MSB is on the right (index=7)
    out["COMMCNT"] = _extract_bits_as_uint8(bits_rev, 2, 8)
    out["COMMCNT_CHECKSUM_FAIL"] = bits_rev[:, 0] > 0
    out["COMMCNT_SYNTAX_ERROR"] = bits_rev[:, 1] > 0

    # Extract fields from the mica L0 8-bit COMMPROG value which is really two MSIDs
    bits = np.unpackbits(out["COMMPROG"].data).reshape(-1, 8)
    bits_rev = bits[:, ::-1]
    out["COMMPROG"] = _extract_bits_as_uint8(bits_rev, 2, 8)
    out["COMMPROG_REPEAT"] = _extract_bits_as_uint8(bits_rev, 0, 2)

    return out


def get_slot_data(
    start,
    stop,
    slot,
    imgsize=None,
    db=None,
    data_root=None,
    columns=None,
    centered_8x8=False,
):
    """
    Get slot data.

    For a the given parameters, retrieve telemetry and construct a
    masked array of the MSIDs available in that telemetry.

      >>> from mica.archive import aca_l0
      >>> slot_data = aca_l0.get_slot_data('2012:001', '2012:002', slot=7)
      >>> temp_ccd_8x8 = aca_l0.get_slot_data('2005:001', '2005:010',
      ...                                     slot=6, imgsize=[8],
      ...                                     columns=['TIME', 'TEMPCCD'])

    :param start: start time of requested interval
    :param stop: stop time of requested interval
    :param slot: slot number integer (in the range 0 -> 7)
    :param imgsize: list of integers of desired image sizes
                    (defaults to all -> [4, 6, 8])
    :param db: handle to archive lookup table
    :param data_root: parent directory that contains archfiles.db3
                      (for use when db handle not available)
    :param columns: list of desired columns in the ACA0 telemetry
                    (defaults to all in 8x8 telemetry)
    :param centered_8x8: boolean flag to reshape the IMGRAW field to (-1, 8, 8)
                         (defaults to False)
    :returns: data structure for slot
    :rtype: numpy masked recarray
    """
    if data_root is None:
        data_root = CONFIG["data_root"]
    if columns is None:
        columns = ACA_DTYPE_NAMES
    if imgsize is None:
        imgsize = [4, 6, 8]
    if db is None:
        dbfile = os.path.join(data_root, "archfiles.db3")
        db = dict(dbi="sqlite", server=dbfile)

    data_files = _get_file_records(
        start, stop, slots=[slot], imgsize=imgsize, db=db, data_root=data_root
    )
    aca_dtype = ACA_DTYPE_8x8 if centered_8x8 else ACA_DTYPE
    dtype = [k for k in aca_dtype if k[0] in columns]

    if not len(data_files):
        # return an empty masked array
        return ma.zeros(0, dtype=dtype)
    rows = np.sum(data_files["rows"])
    zero_row = ma.zeros(1, dtype=dtype)
    zero_row.mask = ma.masked
    all_rows = zero_row.repeat(rows)
    rowcount = 0
    for f in data_files:
        fp = os.path.join(
            data_root, str(f["year"]), "{0:03d}".format(f["doy"]), f["filename"]
        )
        hdu = pyfits.open(fp)
        chunk = hdu[1].data
        idx0, idx1 = rowcount, (rowcount + len(chunk))
        f_imgsize = int(np.sqrt(chunk[0]["IMGRAW"].size))
        for fname in all_rows.dtype.names:
            if fname in ["IMGRAW", "IMGSIZE"]:
                continue
            if fname in chunk.dtype.names:
                all_rows[fname][idx0:idx1] = chunk[fname]
            elif fname == "PIXTLM":
                all_rows[fname][idx0:idx1] = "ORIG"
            elif fname == "BGDTYP":
                all_rows[fname][idx0:idx1] = "FLAT"
        if "IMGSIZE" in columns:
            all_rows["IMGSIZE"][idx0:idx1] = f_imgsize
        if "FILENAME" in columns:
            all_rows["FILENAME"][idx0:idx1] = f["filename"]
        if "IMGRAW" in columns:
            if centered_8x8:
                if f_imgsize == 8:
                    all_rows["IMGRAW"][idx0:idx1] = chunk["IMGRAW"]
                else:
                    i = (8 - f_imgsize) // 2
                    all_rows["IMGRAW"][idx0:idx1, i:-i, i:-i] = chunk["IMGRAW"]
            else:
                all_rows["IMGRAW"].reshape(rows, 8, 8)[
                    idx0:idx1, 0:f_imgsize, 0:f_imgsize
                ] = chunk["IMGRAW"].reshape(len(chunk), f_imgsize, f_imgsize)
        rowcount += len(chunk)

    # just include the rows in the requested time range in the returned data
    oktime = (all_rows["TIME"] >= DateTime(start).secs) & (
        all_rows["TIME"] <= DateTime(stop).secs
    )
    return all_rows[oktime]


def get_l0_images(start, stop, slot, imgsize=None, columns=None):
    """
    Get ACA L0 Images.

    Get ACA L0 images for the given  ``start`` and ``stop`` times and
    the given ``slot``.  Optionally filter on image size via ``imgsize``
    or change the default image metadata via ``columns``.

      >>> from mica.archive import aca_l0
      >>> imgs = aca_l0.get_l0_images('2012:001', '2012:002', slot=7)
      >>> imgs = aca_l0.get_l0_images('2005:001', '2005:002', slot=6, imgsize=[8])

    The default columns are:
    ['TIME', 'IMGROW0', 'IMGCOL0', 'BGDAVG', 'IMGSTAT', 'IMGFUNC1', 'IMGSIZE', 'IMGSCALE', 'INTEG']

    The image pixel values are given in units of DN.  One can convert to e-/sec
    by multiplying by (5 / INTEG).

    :param start: start time of requested interval
    :param stop: stop time of requested interval
    :param slot: slot number integer (in the range 0 -> 7)
    :param imgsize: list of integers of desired image sizes (default=[4, 6, 8])
    :param columns: image meta-data columns

    :returns: list of ACAImage objects
    """
    if columns is None:
        columns = [
            "TIME",
            "BGDAVG",
            "IMGSTAT",
            "IMGFUNC1",
            "IMGSIZE",
            "IMGSCALE",
            "INTEG",
        ]
    if "IMGROW0" not in columns:
        columns.append("IMGROW0")
    if "IMGCOL0" not in columns:
        columns.append("IMGCOL0")

    # It doesn't make any sense to include IMGRAW in the meta-data columns
    # of ACAImages so just remove it if the user has requested that.
    if "IMGRAW" in columns:
        columns.remove("IMGRAW")

    slot_columns = list(set(columns + ["QUALITY", "IMGRAW", "IMGSIZE"]))
    dat = get_slot_data(start, stop, slot, imgsize=imgsize, columns=slot_columns)

    ok = dat["QUALITY"] == 0
    if not (np.all(ok)):
        dat = dat[ok]

    # Convert temperatures from K to degC
    for temp in ("TEMPCCD", "TEMPHOUS", "TEMPPRIM", "TEMPSEC"):
        if temp in dat.dtype.names:
            dat[temp] -= 273.15

    dat["IMGRAW"].set_fill_value(-9999)
    dat = dat[~dat["IMGSIZE"].mask]
    imgraws = dat["IMGRAW"].filled()

    imgs = []
    for row, raw_imgraw in zip(dat, imgraws):
        imgraw = raw_imgraw.reshape(8, 8)
        sz = row["IMGSIZE"]
        if sz < 8:
            imgraw = imgraw[:sz, :sz]
        meta = {name: row[name] for name in columns if not np.any(row.mask[name])}
        imgs.append(ACAImage(imgraw, meta=meta))

    return imgs


class MSID(object):
    def __init__(self, msid, slot, start, stop):
        self.tstart = DateTime(start).secs
        self.tstop = DateTime(stop).secs
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        self.slot = slot
        self._check_msid(msid)
        self._get_data()

    def _check_msid(self, req_msid):
        if req_msid.upper() in ACA_DTYPE_NAMES:
            self.msid = req_msid.lower()
        else:
            raise MissingDataError("msid %s not found" % req_msid)

    def _get_data(self):
        slot_data = get_slot_data(
            self.tstart,
            self.tstop,
            self.slot,
            imgsize=[8],
            columns=["TIME", self.msid.upper()],
        )
        self.vals = slot_data[self.msid.upper()]
        self.times = slot_data["TIME"]


class MSIDset(collections.OrderedDict):
    def __init__(self, msids, start, stop):
        super(MSIDset, self).__init__()
        self.tstart = DateTime(start).secs
        self.tstop = DateTime(stop).secs
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        for msid in msids:
            self[msid] = MSID(msid, self.tstart, self.tstop)


def obsid_times(obsid):
    from kadi import events  # Kadi is a big import so defer

    dwells = events.dwells.filter(obsid=obsid)
    n_dwells = len(dwells)
    tstart = dwells[0].tstart
    tstop = dwells[n_dwells - 1].tstop

    return tstart, tstop


def get_files(
    obsid=None, start=None, stop=None, slots=None, imgsize=None, db=None, data_root=None
):
    """
    Get list of ACA0 files.

    Retrieve list of files from ACA0 archive lookup table that
    match arguments.  The database query returns files with

      tstart < stop
      and
      tstop > start

    which returns all files that contain any part of the interval
    between start and stop.  If the obsid argument is provided, the
    archived obspar tstart/tstop (sybase aca.obspar table) are used.

      >>> from mica.archive import aca_l0
      >>> obsid_files = aca_l0.get_files(obsid=5438)
      >>> time_files = aca_l0.get_files(start='2012:001', stop='2012:002')
      >>> time_8x8 = aca_l0.get_files(start='2011:001', stop='2011:010',
      ...                             imgsize=[8])

    :param obsid: obsid
    :param start: start time of requested interval
    :param stop: stop time of requested interval
    :param slots: list of integers of desired image slots to retrieve
                   (defaults to all -> [0, 1, 2, 3, 4, 5, 6, 7, 8])
    :param imgsize: list of integers of desired image sizes
                    (defaults to all -> [4, 6, 8])
    :param db: handle to archive lookup table
    :param data_root: parent directory of Ska aca l0 archive

    :returns: interval files
    :rtype: list
    """
    if slots is None:
        slots = [0, 1, 2, 3, 4, 5, 6, 7]
    if data_root is None:
        data_root = CONFIG["data_root"]
    if imgsize is None:
        imgsize = [4, 6, 8]
    if db is None:
        dbfile = os.path.join(data_root, "archfiles.db3")
        db = dict(dbi="sqlite", server=dbfile)
    if obsid is None:
        if start is None or stop is None:
            raise TypeError("Must supply either obsid or start and stop")
    else:
        start, stop = obsid_times(obsid)

    file_records = _get_file_records(
        start, stop, slots=slots, imgsize=imgsize, db=db, data_root=data_root
    )
    files = [
        os.path.join(
            data_root, "%04d" % f["year"], "%03d" % f["doy"], str(f["filename"])
        )
        for f in file_records
    ]
    return files


def _get_file_records(
    start, stop=None, slots=None, imgsize=None, db=None, data_root=None
):
    """
    Get list of ACA0 files records.

    Retrieve list of files from ACA0 archive lookup table that
    match arguments.  The database query returns files with

      tstart < stop
      and
      tstop > start

    which returns all files that contain any part of the interval
    between start and stop.

    :param start: start time of requested interval
    :param stop: stop time of requested interval
                 (default of None will get DateTime(None) which is
                 equivalent to 'now')
    :param slots: list of integers of desired image slots to retrieve
                   (defaults to all -> [0, 1, 2, 3, 4, 5, 6, 7, 8])
    :param imgsize: list of integers of desired image sizes
                    (defaults to all -> [4, 6, 8])
    :param db: handle to archive lookup table
    :param data_root: parent directory of Ska aca l0 archive

    :returns: interval files
    :rtype: list
    """
    if slots is None:
        slots = [0, 1, 2, 3, 4, 5, 6, 7]
    if data_root is None:
        data_root = CONFIG["data_root"]
    if imgsize is None:
        imgsize = [4, 6, 8]
    if db is None:
        dbfile = os.path.join(data_root, "archfiles.db3")
        db = dict(dbi="sqlite", server=dbfile)

    tstart = DateTime(start).secs
    tstop = DateTime(stop).secs
    imgsize_str = ",".join([str(x) for x in imgsize])
    slot_str = ",".join([str(x) for x in slots])
    # a composite index isn't as fast as just doing a padded search on one
    # index first (tstart).  This gets extra files to make sure we don't
    # miss the case where tstart is in the middle of an interval, but
    # drastically reduces the size of the bsearch on the tstop index
    tstart_pad = 10 * 86400
    db_query = (
        "SELECT * FROM archfiles "
        "WHERE tstart >= %f - %f "
        "AND tstart < %f "
        "AND tstop > %f "
        "AND slot in (%s) "
        "AND imgsize in (%s) "
        "order by filetime asc "
        % (tstart, tstart_pad, tstop, tstart, slot_str, imgsize_str)
    )
    with ska_dbi.DBI(**db) as db:
        files = db.fetchall(db_query)
    return files


class Updater(object):
    def __init__(
        self,
        db=None,
        data_root=None,
        temp_root=None,
        days_at_once=None,
        sql_def=None,
        cda_table=None,
        filetype=None,
        start=None,
        stop=None,
    ):
        for init_opt in [
            "data_root",
            "temp_root",
            "days_at_once",
            "sql_def",
            "cda_table",
        ]:
            setattr(self, init_opt, vars()[init_opt] or CONFIG[init_opt])
        self.data_root = os.path.abspath(self.data_root)
        self.temp_root = os.path.abspath(self.temp_root)
        self.filetype = filetype or FILETYPE
        if db is None:
            dbfile = os.path.join(data_root, "archfiles.db3")
            db = dict(dbi="sqlite", server=dbfile, autocommit=False)
        self.db = db
        self.start = start
        self.stop = stop

    # def _rebuild_database(db=None, db_file=None,
    #                      data_root=config['data_root'],
    #                      sql_def=config['sql_def']):
    #    """
    #    Utility routine to rebuild the file lookup database using the
    #    package defaults and the files in the archive.
    #    """
    #    if db is None and db_file is None:
    #        raise ValueError
    #    if db is None and db_file is not None:
    #        logger.info("creating archfiles db from %s"
    #                    % sql_def)
    #        db_sql = os.path.join(os.environ['SKA_DATA'],
    #                              'mica', sql_def)
    #        db_init_cmds = file(db_sql).read()
    #        db = ska_dbi.DBI(dbi='sqlite', server=db_file,
    #                         autocommit=False)
    #        db.execute(db_init_cmds, commit=True)
    #    year_dirs = sorted(glob(
    #            os.path.join(data_root, '[12][0-9][0-9][0-9]')))
    #    for ydir in year_dirs:
    #        day_dirs = sorted(glob(
    #                os.path.join(ydir, '[0-3][0-9][0-9]')))
    #        for ddir in day_dirs:
    #            archfiles = sorted(glob(
    #                    os.path.join(ddir, '*_img0*')))
    #            db.execute("begin transaction")
    #            for i, f in enumerate(archfiles):
    #                arch_info = read_archfile(i, f, archfiles, db)
    #                if arch_info:
    #                    db.insert(arch_info, 'archfiles')
    #            db.commit()

    def _get_missing_archive_files(self, start, only_new=False):
        ingested_files = self._get_arc_ingested_files()
        startdate = DateTime(start).date
        logger.info("Checking for missing files from %s" % startdate)
        # find the index in the cda archive list that matches
        # the first entry with the "start" date
        for idate, backcnt in zip(ingested_files["ingest_date"][::-1], count(1)):  # noqa: B007 backcnt not used in loop
            if idate < startdate:
                break

        last_ok_date = None
        missing = []
        # for the entries after the start date, see if we have the
        # file or a later version
        with ska_dbi.DBI(**self.db) as db:
            for file, idx in zip(ingested_files[-backcnt:], count(0)):
                filename = file["filename"]
                db_match = db.fetchall(
                    "select * from archfiles where "
                    + "filename = '%s' or filename = '%s.gz'" % (filename, filename)
                )
                if len(db_match):
                    continue
                file_re = re.search(
                    r"acaf(\d+)N(\d{3})_(\d)_img0.fits(\.gz)?", filename
                )
                if not file_re:
                    continue
                slot = int(file_re.group(3))
                filetime = int(file_re.group(1))
                version = int(file_re.group(2))
                # if there is an old file and we're just looking for new ones
                range_match = db.fetchall(
                    """SELECT * from archfiles
                       WHERE filetime = %(filetime)d
                       and slot = %(slot)d"""
                    % dict(filetime=filetime, slot=slot)
                )
                if range_match and only_new:
                    continue
                # if there is a newer file there already
                version_match = db.fetchall(
                    """SELECT * from archfiles
                       WHERE filetime = %(filetime)d
                       and slot = %(slot)d
                       and revision >= %(version)d"""
                    % dict(slot=slot, version=version, filetime=filetime)
                )
                if version_match:
                    continue
                # and if made it this far add the file to the list
                missing.append(file)
                # and update the date through which data is complete to
                # the time of the previous file ingest
                if last_ok_date is None:
                    last_ok_date = ingested_files["ingest_date"][-backcnt:][idx - 1]

        if last_ok_date is None:
            last_ok_date = ingested_files["ingest_date"][-1]
        return missing, last_ok_date

    def _get_arc_ingested_files(self):
        table_file = os.path.join(self.data_root, self.cda_table)
        with tables.open_file(table_file) as h5f:
            tbl = h5f.get_node("/", "data")
            arc_files = tbl[:]
        return Table(arc_files)

    def _get_archive_files(self, start, stop):
        """
        Update FITS file archive with arc5gl and ingest files into file archive

        :param filetype: a dictionary (or dictionary-like) object with keys for
                         'level', 'instrum', 'content', 'arc5gl_query'
                         and 'fileglob' for arc5gl.  For ACA0:
                         {'level': 'L0', 'instrum': 'PCAD',
                         'content': 'ACADATA', 'arc5gl_query': 'ACA0',
                         'fileglob': 'aca*fits*'}
        :param start: start of interval to retrieve (Chandra.Time compatible)
        :param stop: end of interval to retrieve (Chandra.Time compatible)

        :returns: retrieved file names
        :rtype: list
        """

        filetype = self.filetype
        # Retrieve CXC archive files in a temp directory with arc5gl
        arc5 = Ska.arc5gl.Arc5gl()
        arc5.sendline("tstart=%s" % DateTime(start).date)
        arc5.sendline("tstop=%s" % DateTime(stop).date)
        arc5.sendline("get %s" % filetype["arc5gl_query"].lower())
        return sorted(glob(filetype["fileglob"]))

    def _read_archfile(self, i, f, archfiles):
        """
        Read index fits file to get lookup data.

        Read FITS filename ``f`` with index ``i`` (position within list of
        filenames) and get dictionary of values to store in file lookup
        database.  These values include all header items in
        ``ARCHFILES_HDR_COLS`` plus the header checksum, the image slot
        as determined by the filename, the imagesize as determined
        by IMGRAW, the year, day-of-year, and number of data rows in the file.

        :param i: index of file f within list of files archfiles
        :param f: filename
        :param archfiles: list of filenames for this batch
        :param db: database handle for file lookup database (ska_dbi handle)

        :returns: info for a file.
        :rtype: dictionary
        """

        # Check if filename is already in file lookup table
        # If so then delete temporary file and abort further processing.
        filename = os.path.basename(f)
        with ska_dbi.DBI(**self.db) as db:
            if db.fetchall(
                "SELECT filename FROM archfiles WHERE filename=?", (filename,)
            ):
                logger.debug(
                    "File %s already in archfiles - unlinking and skipping" % f
                )
                os.unlink(f)
                return None

        logger.debug("Reading (%d / %d) %s" % (i, len(archfiles), filename))
        hdus = pyfits.open(f)
        hdu = hdus[1]

        # Accumlate relevant info about archfile that will be ingested
        # (this is borrowed from eng-archive but is set to archive to a
        # database table instead of an h5 table in this version)
        archfiles_row = dict((x, hdu.header.get(x.upper())) for x in ARCHFILES_HDR_COLS)
        archfiles_row["checksum"] = hdu.header.get("checksum") or hdu._checksum
        imgsize = hdu.data[0]["IMGRAW"].shape[0]
        archfiles_row["imgsize"] = int(imgsize)
        archfiles_row["slot"] = int(
            re.search(r"acaf\d+N\d{3}_(\d)_img0.fits(\.gz)?", filename).group(1)
        )
        archfiles_row["filename"] = filename
        archfiles_row["filetime"] = int(
            re.search(r"(\d+)", archfiles_row["filename"]).group(1)
        )
        filedate = DateTime(archfiles_row["filetime"]).date
        year, doy = (
            int(x) for x in re.search(r"(\d\d\d\d):(\d\d\d)", filedate).groups()
        )
        archfiles_row["year"] = year
        archfiles_row["doy"] = doy
        archfiles_row["rows"] = len(hdu.data)
        hdus.close()

        with ska_dbi.DBI(**self.db) as db:
            # remove old versions of this file
            oldmatches = db.fetchall(
                """SELECT * from archfiles
                   WHERE filetime = %(filetime)d
                   and slot = %(slot)d
                   and startmjf = %(startmjf)d and startmnf = %(startmnf)d
                   and stopmjf = %(stopmjf)d and stopmnf = %(stopmnf)d """
                % archfiles_row
            )
            if len(oldmatches):
                self._arch_remove(oldmatches)

        interval_matches = _get_file_records(
            archfiles_row["tstart"],
            archfiles_row["tstop"],
            slots=[archfiles_row["slot"]],
        )
        if len(interval_matches):
            # if there are files there that still overlap the new file
            # and they are all older, remove the old files.
            if np.all(interval_matches["revision"] < archfiles_row["revision"]):
                logger.info("removing overlapping files at older revision(s)")
                logger.info(interval_matches)
                self._arch_remove(interval_matches)
            # if the overlapping files are all from the same revision
            # just ingest them and hope for the best
            elif np.all(interval_matches["revision"] == archfiles_row["revision"]):
                logger.info("ignoring overlap for same process revision")
            # if the files that overlap are all newer, let's not ingest
            # the "missing" file
            elif np.all(interval_matches["revision"] > archfiles_row["revision"]):
                return None
            else:
                logger.error(archfiles_row)
                logger.error(interval_matches)
                # throw an error if there is still overlap
                raise ValueError("Cannot ingest %s, overlaps existing files" % filename)

        return archfiles_row

    def _arch_remove(self, defunct_matches):
        with ska_dbi.DBI(**self.db) as db:
            for file_record in defunct_matches:
                query = (
                    """delete from archfiles
                              WHERE filetime = %(filetime)d
                              and slot = %(slot)d
                              and startmjf = %(startmjf)d
                              and startmnf = %(startmnf)d
                              and stopmjf = %(stopmjf)d
                              and stopmnf = %(stopmnf)d """
                    % file_record
                )
                logger.info(query)
                db.execute(query)
                db.commit()
                archdir = os.path.abspath(
                    os.path.join(
                        self.data_root,
                        str(file_record["year"]),
                        "{0:03d}".format(file_record["doy"]),
                    )
                )
                logger.info(
                    "deleting %s" % os.path.join(archdir, file_record["filename"])
                )
                real_file = os.path.join(archdir, file_record["filename"])
                if os.path.exists(real_file):
                    os.unlink(real_file)

    def _move_archive_files(self, archfiles):
        """
        Move archive files where they should go.

        Move ACA L0 files into the file archive into directories
        by YYYY/DOY under the specified data_root
        """

        data_root = self.data_root
        if not os.path.exists(data_root):
            os.makedirs(data_root)
        for f in archfiles:
            if not os.path.exists(f):
                continue
            basename = os.path.basename(f)
            # use the timestamp right from the ACA0 filename
            tstart = re.search(r"(\d+)", str(basename)).group(1)
            datestart = DateTime(tstart).date
            year, doy = re.search(r"(\d\d\d\d):(\d\d\d)", datestart).groups()
            archdir = os.path.abspath(os.path.join(data_root, year, doy))
            # construct the destination filepath/name
            archfile = os.path.abspath(os.path.join(archdir, basename))
            if not os.path.exists(archdir):
                os.makedirs(archdir)
            if not os.path.exists(archfile):
                logger.debug("mv %s %s" % (os.path.abspath(f), archfile))
                os.chmod(f, 0o775)
                shutil.move(f, archfile)
            if os.path.exists(f):
                logger.info("Unlinking %s" % os.path.abspath(f))
                os.unlink(f)

    def _fetch_by_time(self, range_tstart, range_tstop):
        logger.info(
            "Fetching %s from %s to %s"
            % ("ACA L0 Data", DateTime(range_tstart).date, DateTime(range_tstop).date)
        )
        archfiles = self._get_archive_files(
            DateTime(range_tstart), DateTime(range_tstop)
        )
        return archfiles

    def _fetch_individual_files(self, files):
        arc5 = Ska.arc5gl.Arc5gl(echo=True)
        logger.info(
            "********** %s %s **********" % (self.filetype["content"], time.ctime())
        )
        fetched_files = []
        ingest_dates = []
        # get the files, store in file archive, and record in database
        for file in files:
            # Retrieve CXC archive files in a temp directory with arc5gl
            missed_file = file["filename"]
            arc5.sendline("dataset=flight")
            arc5.sendline("detector=pcad")
            arc5.sendline("subdetector=aca")
            arc5.sendline("level=0")
            arc5.sendline("filename=%s" % missed_file)
            arc5.sendline("version=last")
            arc5.sendline("operation=retrieve")
            arc5.sendline("go")
            have_files = sorted(glob(f"{missed_file}*"))
            filename = have_files[0]
            # if it isn't gzipped, just gzip it
            if re.match(r".*\.fits$", filename):
                f_in = open(file, "rb")
                f_out = gzip.open("%s.gz" % filename, "wb")
                f_out.writelines(f_in)
                f_out.close()
                f_in.close()
                filename = "%s.gz" % filename
            fetched_files.append(filename)
            ingest_dates.append(file["ingest_date"])
        return fetched_files, ingest_dates

    def _insert_files(self, files):
        count_inserted = 0
        for i, f in enumerate(files):
            arch_info = self._read_archfile(i, f, files)
            if arch_info:
                self._move_archive_files([f])
                with ska_dbi.DBI(**self.db) as db:
                    db.insert(arch_info, "archfiles")
                    db.commit()
                count_inserted += 1
        logger.info("Ingested %d files" % count_inserted)

    def update(self):
        """
        Run ACA0 update process.

        Retrieve ACA0 telemetry files from the CXC archive, store in the
        Ska/ACA archive, and update database of files.
        """
        contentdir = self.data_root
        if not os.path.exists(contentdir):
            os.makedirs(contentdir)
        if not os.path.exists(self.temp_root):
            os.makedirs(self.temp_root)
        archdb = os.path.join(contentdir, "archfiles.db3")
        # if the database of the archived files does not exist,
        # or is empty, make it
        if not os.path.exists(archdb) or os.stat(archdb).st_size == 0:
            logger.info("creating archfiles db from %s" % self.sql_def)
            db_sql = Path(__file__).parent / self.sql_def
            db_init_cmds = open(db_sql).read()
            with ska_dbi.DBI(**self.db) as db:
                db.execute(db_init_cmds, commit=True)
        if self.start:
            datestart = DateTime(self.start)
        else:
            # Get datestart as the most-recent file time from archfiles table
            # will need min-of-max-slot-datestart
            with ska_dbi.DBI(**self.db) as db:
                last_time = min(
                    [
                        db.fetchone(
                            "select max(filetime) from archfiles where slot = %d" % s
                        )["max(filetime)"]
                        for s in range(0, 8)
                    ]
                )
                if last_time is None:
                    raise ValueError(
                        "No files in archive to do update-since-last-run mode.\n"
                        "Please specify a time with --start"
                    )
                datestart = DateTime(last_time)
        datestop = DateTime(self.stop)
        padding_seconds = 10000
        # loop over the specified time range in chunks of
        # days_at_once in seconds with some padding
        for tstart in np.arange(
            datestart.day_start().secs,
            datestop.day_end().secs,
            self.days_at_once * 86400,
        ):
            # set times for a chunk
            range_tstart = tstart - padding_seconds
            range_tstop = tstart + self.days_at_once * 86400
            if range_tstop > datestop.day_end().secs:
                range_tstop = datestop.day_end().secs
            range_tstop += padding_seconds
            # make a temporary directory
            tmpdir = Ska.File.TempDir(dir=self.temp_root)
            dirname = tmpdir.name
            logger.debug("Files save to temp dir %s" % dirname)
            # get the files, store in file archive, and record in database
            with Ska.File.chdir(dirname):
                fetched_files = self._fetch_by_time(range_tstart, range_tstop)
                self._insert_files(fetched_files)

        timestamp_file = os.path.join(self.data_root, "last_timestamp.txt")
        # get list of missing files since the last time the tool ingested
        # files.  If this is first run of the tool, check from the start of
        # the requested time range
        if os.path.exists(timestamp_file) and os.stat(timestamp_file).st_size > 0:
            cda_checked_timestamp = open(timestamp_file).read().rstrip()
        else:
            cda_checked_timestamp = DateTime(self.start).date
        missing_datetime = DateTime(cda_checked_timestamp)
        missing_files, last_ingest_date = self._get_missing_archive_files(
            missing_datetime, only_new=True
        )
        # update the file to have up through the last confirmed good file
        # even before we try to fetch missing ones
        open(timestamp_file, "w").write("%s" % last_ingest_date)

        if len(missing_files):
            logger.info("Found %d missing individual files" % len(missing_files))
            # make a temporary directory
            tmpdir = Ska.File.TempDir(dir=self.temp_root)
            dirname = tmpdir.name
            logger.info("File save to temp dir %s" % dirname)
            with Ska.File.chdir(dirname):
                fetched_files, ingest_times = self._fetch_individual_files(
                    missing_files
                )
                self._insert_files(fetched_files)

            last_ingest_date = missing_files[-1]["ingest_date"]
            # update the file to have up through the last confirmed good file
            # even before we try to fetch missing ones
            open(timestamp_file, "w").write("%s" % last_ingest_date)
        else:
            logger.info("No missing files")


def main():
    """
    ACA L0 cmdline interface.

    Command line interface to fetch ACA L0 telemetry from the CXC Archive
    and store it in the Ska archive.
    """
    opt = get_options()
    kwargs = vars(opt)
    updater = Updater(**kwargs)
    updater.update()


if __name__ == "__main__":
    main()
