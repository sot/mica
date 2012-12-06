import os
from glob import glob
import re
import logging
import shutil
import time
import pyfits
import numpy as np
import numpy.ma as ma
import argparse
import collections
import tables
from itertools import izip, count

import Ska.DBI
import Ska.arc5gl
from Chandra.Time import DateTime
import Ska.File


logger = logging.getLogger('aca0 fetch')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

# borrowed from eng_archive
archfiles_hdr_cols = ('tstart', 'tstop', 'startmjf', 'startmnf',
                      'stopmjf', 'stopmnf',
                      'tlmver', 'ascdsver', 'revision', 'date',
                      'imgsize')


dtype = [('level', '|S4'), ('instrum', '|S7'), ('content', '|S13'),
         ('arc5gl_query', '|S27'), ('fileglob', '|S9')]
filetype = np.rec.fromrecords([('L0', 'PCAD', 'ACADATA', 'ACA0', '*fits.gz')],
                               dtype=dtype)[0]

aca_dtype = (('TIME', '>f8'), ('QUALITY', '>i4'), ('MJF', '>i4'),
             ('MNF', '>i4'),
             ('END_INTEG_TIME', '>f8'), ('INTEG', '>f4'), ('GLBSTAT', '|u1'),
             ('COMMCNT', '|u1'), ('COMMPROG', '|u1'), ('IMGFID1', '|u1'),
             ('IMGNUM1', '|u1'), ('IMGFUNC1', '|u1'), ('IMGSTAT', '|u1'),
             ('IMGROW0', '>i2'), ('IMGCOL0', '>i2'), ('IMGSCALE', '>i2'),
             ('BGDAVG', '>i2'), ('IMGFID2', '|u1'), ('IMGNUM2', '|u1'),
             ('IMGFUNC2', '|u1'), ('BGDRMS', '>i2'), ('TEMPCCD', '>f4'),
             ('TEMPHOUS', '>f4'), ('TEMPPRIM', '>f4'), ('TEMPSEC', '>f4'),
             ('BGDSTAT', '|u1'), ('IMGFID3', '|u1'), ('IMGNUM3', '|u1'),
             ('IMGFUNC3', '|u1'), ('IMGFID4', '|u1'), ('IMGNUM4', '|u1'),
             ('IMGFUNC4', '|u1'), ('IMGRAW', '>f4', (64,)),
             ('HD3TLM62', '|u1'),
             ('HD3TLM63', '|u1'), ('HD3TLM64', '|u1'), ('HD3TLM65', '|u1'),
             ('HD3TLM66', '|u1'), ('HD3TLM67', '|u1'), ('HD3TLM72', '|u1'),
             ('HD3TLM73', '|u1'), ('HD3TLM74', '|u1'), ('HD3TLM75', '|u1'),
             ('HD3TLM76', '|u1'), ('HD3TLM77', '|u1'),
             ('IMGSIZE', '>i4'))

aca_dtype_names = tuple([k[0] for k in aca_dtype])

config = dict(data_root='/data/aca/archive/aca0',
              temp_root='/data/aca/archive/temp',
              days_at_once=30.0,
              sql_def='archfiles_aca_l0_def.sql',
              cda_table='cda_aca0.h5')


def get_options():
    parser = argparse.ArgumentParser(
        description="Fetch aca level 0 products and make a file archive")
    defaults = dict(config)
    parser.set_defaults(**defaults)
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--temp-root",
                        help="parent temp directory")
    parser.add_argument("--start",
                        help="start date for retrieve "
                        + "(defaults to max date of archived files)")
    parser.add_argument("--stop",
                        help="stop date for retrieve "
                        + "(defaults to now)")
    parser.add_argument("--days-at-once",
                        type=float,
                        help="if over this number, "
                        + "bin to chunks of this number of days")
    parser.add_argument("--cda-table",
                        help="file name for h5 list from Chandra Archive")
    opt = parser.parse_args()
    return opt


def get_slot_data(start, stop, slot, imgsize=None,
                  db=None, data_root=None, columns=None,
                  ):
    """
    For a the given parameters, retrieve telemetry and construct a
    masked array of the MSIDs available in that telemetry.

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

    :returns: data structure for slot
    :rtype: numpy masked recarray
    """
    if data_root is None:
        data_root = config['data_root']
    if columns is None:
        columns = aca_dtype_names
    if imgsize is None:
        imgsize = [4, 6, 8]
    if db is None:
        dbfile = os.path.join(data_root, 'archfiles.db3')
        db = Ska.DBI.DBI(dbi='sqlite', server=dbfile)

    data_files = get_interval_files(start, stop, slots=[slot],
                                    imgsize=imgsize, db=db,
                                    data_root=data_root)
    dtype = [k for k in aca_dtype if k[0] in columns]
    if not len(data_files):
        # return an empty masked array
        return ma.zeros(0, dtype=dtype)
    rows = np.sum(data_files['rows'])
    zero_row = ma.zeros(1, dtype=dtype)
    zero_row.mask = ma.masked
    all_rows = zero_row.repeat(rows)
    rowcount = 0
    for f in data_files:
        fp = os.path.join(data_root,
                          str(f['year']),
                          "{0:03d}".format(f['doy']),
                          f['filename'])
        hdu = pyfits.open(fp)
        chunk = hdu[1].data
        for fname in all_rows.dtype.names:
            if fname == 'IMGRAW' or fname == 'IMGSIZE':
                continue
            if fname in chunk.dtype.names:
                all_rows[fname][rowcount:(rowcount + len(chunk))] \
                    = chunk.field(fname)
        if 'IMGSIZE' in columns and 'IMGRAW' in columns:
            f_imgsize = int(np.sqrt(chunk[0].field('IMGRAW').size))
            all_rows['IMGSIZE'][rowcount:(rowcount + len(chunk))] = f_imgsize
            all_rows['IMGRAW'].reshape(rows, 8, 8)[
                rowcount:(rowcount + len(chunk)), 0:f_imgsize, 0:f_imgsize] = (
                chunk.field('IMGRAW').reshape(len(chunk),
                                              f_imgsize, f_imgsize))
        rowcount += len(chunk)

    # just include the rows in the requested time range in the returned data
    oktime = ((all_rows['TIME'] >= DateTime(start).secs)
               & (all_rows['TIME'] <= DateTime(stop).secs))
    return all_rows[oktime]


class MissingDataError(Exception):
    pass


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
        if req_msid.upper() in aca_dtype_names:
            self.msid = req_msid.lower()
        else:
            raise MissingDataError("msid %s not found" % req_msid)

    def _get_data(self):
        slot_data = get_slot_data(
            self.tstart, self.tstop, self.slot, imgsize=[8],
            columns=['TIME', self.msid.upper()])
        self.vals = slot_data[self.msid.upper()]
        self.times = slot_data['TIME']


class MSIDset(collections.OrderedDict):
    def __init__(self, msids, start, stop):
        super(MSIDset, self).__init__()
        self.tstart = DateTime(start).secs
        self.tstop = DateTime(stop).secs
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        for msid in msids:
            self[msid] = MSID(msid, self.tstart, self.tstop)


def get_interval_files(start, stop=None, slots=None,
                       imgsize=None, db=None, data_root=None):
    """
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
        data_root = config['data_root']
    if imgsize is None:
        imgsize = [4, 6, 8]
    if db is None:
        dbfile = os.path.join(data_root, 'archfiles.db3')
        db = Ska.DBI.DBI(dbi='sqlite', server=dbfile)

    tstart = DateTime(start).secs
    tstop = DateTime(stop).secs
    imgsize_str = ','.join([str(x) for x in imgsize])
    slot_str = ','.join([str(x) for x in slots])
    # a composite index isn't as fast as just doing a padded search on one
    # index first (tstart).  This gets extra files to make sure we don't
    # miss the case where tstart is in the middle of an interval, but
    # drastically reduces the size of the bsearch on the tstop index
    tstart_pad = 10 * 86400
    db_query = ('SELECT * FROM archfiles '
                'WHERE tstart >= %f - %f '
                'AND tstart < %f '
                'AND tstop > %f '
                'AND slot in (%s) '
                'AND imgsize in (%s) '
                'order by filetime asc '
                % (tstart, tstart_pad, tstop, tstart, slot_str, imgsize_str))
    files = db.fetchall(db_query)
    return files


class Updater(object):
    def __init__(self,
                 db=None,
                 data_root=config['data_root'],
                 temp_root=config['temp_root'],
                 days_at_once=config['days_at_once'],
                 sql_def=config['sql_def'],
                 cda_table=config['cda_table'],
                 filetype=filetype,
                 start=None,
                 stop=None,
                 ):
        self.data_root = data_root
        self.temp_root = temp_root
        self.days_at_once = days_at_once
        self.sql_def = sql_def
        self.cda_table = cda_table
        self.filetype = filetype
        if not db:
            dbfile = os.path.join(data_root, 'archfiles.db3')
            db = Ska.DBI.DBI(dbi='sqlite', server=dbfile, autocommit=False)
        self.db = db
        self.start = start
        self.stop = stop

#def _rebuild_database(db=None, db_file=None,
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
#        db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
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
        # find the index in the cda archive list that matches
        # the first entry with the "start" date
        for idate, backcnt in izip(ingested_files[::-1]['ingest_date'],
                                   count(1)):
            if idate < startdate:
                break

        last_ok_date = None
        missing = []
        # for the entries after the start date, see if we have the
        # file or a later version
        for file, idx in izip(ingested_files[-backcnt:], count(0)):
            filename = file['filename']
            db_match = self.db.fetchall(
                "select * from archfiles where "
                + "filename = '%s' or filename = '%s.gz'"
                % (filename, filename))
            if len(db_match):
                continue
            file_re = re.search(
                r'acaf(\d+)N(\d{3})_(\d)_img0.fits(\.gz)?',
                filename)
            if not file_re:
                continue
            slot = int(file_re.group(3))
            filetime = int(file_re.group(1))
            version = int(file_re.group(2))
            # if there is an old file and we're just looking for new ones
            range_match = self.db.fetchall(
                """SELECT * from archfiles
                   WHERE filetime = %(filetime)d
                   and slot = %(slot)d""" % dict(filetime=filetime, slot=slot))
            if range_match and only_new:
                continue
            # if there is a newer file there already
            version_match = self.db.fetchall(
                """SELECT * from archfiles
                   WHERE filetime = %(filetime)d
                   and slot = %(slot)d
                   and revision >= %(version)d"""
                % dict(slot=slot, version=version, filetime=filetime))
            if version_match:
                continue
            # and if made it this far add the file to the list
            missing.append(file)
            # and update the date through which data is complete to
            # the time of the previous file ingest
            if last_ok_date is None:
                last_ok_date = \
                    ingested_files[-backcnt:][idx - 1]['ingest_date']

        if last_ok_date is None:
            last_ok_date = ingested_files[-1]['ingest_date']
        return missing, last_ok_date

    def _get_arc_ingested_files(self):
        table_file = os.path.join(self.data_root, self.cda_table)
        h5f = tables.openFile(table_file)
        tbl = h5f.getNode('/', 'data')
        arc_files = tbl[:]
        h5f.close()
        return arc_files

    def _get_archive_files(self, start, stop):
        """
        Update FITS file archive with arc5gl and ingest files into file archive

        :param filetype: a dictionary (or dictionary-like) object with keys for
                         'level', 'instrum', 'content', 'arc5gl_query'
                         and 'fileglob' for arc5gl.  For ACA0:
                         {'level': 'L0', 'instrum': 'PCAD',
                         'content': 'ACADATA', 'arc5gl_query': 'ACA0',
                         'fileglob': '*fits.gz'}
        :param start: start of interval to retrieve (Chandra.Time compatible)
        :param stop: end of interval to retrieve (Chandra.Time compatible)

        :returns: retrieved file names
        :rtype: list
        """

        filetype = self.filetype
        # Retrieve CXC archive files in a temp directory with arc5gl
        arc5 = Ska.arc5gl.Arc5gl()
        arc5.sendline('tstart=%s' % DateTime(start).date)
        arc5.sendline('tstop=%s' % DateTime(stop).date)
        arc5.sendline('get %s' % filetype['arc5gl_query'].lower())

        return sorted(glob(filetype['fileglob']))

    def _read_archfile(self, i, f, archfiles):
        """
        Read FITS filename ``f`` with index ``i`` (position within list of
        filenames) and get dictionary of values to store in file lookup
        database.  These values include all header items in
        ``archfiles_hdr_cols`` plus the header checksum, the image slot
        as determined by the filename, the imagesize as determined
        by IMGRAW, the year, day-of-year, and number of data rows in the file.

        :param i: index of file f within list of files archfiles
        :param f: filename
        :param archfiles: list of filenames for this batch
        :param db: database handle for file lookup database (Ska.DBI handle)

        :returns: info for a file.
        :rtype: dictionary
        """

        db = self.db
        # Check if filename is already in file lookup table
        # If so then delete temporary file and abort further processing.
        filename = os.path.basename(f)
        if db.fetchall('SELECT filename FROM archfiles WHERE filename=?',
                       (filename,)):
            logger.debug(
                'File %s already in archfiles - unlinking and skipping' % f)
            os.unlink(f)
            return None

        logger.info('Reading (%d / %d) %s' % (i, len(archfiles), filename))
        hdus = pyfits.open(f)
        hdu = hdus[1]

        # Accumlate relevant info about archfile that will be ingested
        # (this is borrowed from eng-archive but is set to archive to a
        # database table instead of an h5 table in this version)
        archfiles_row = dict((x, hdu.header.get(x.upper()))
                             for x in archfiles_hdr_cols)
        archfiles_row['checksum'] = hdu.header.get('checksum') or hdu._checksum
        imgsize = hdu.data[0].field('IMGRAW').shape[0]
        archfiles_row['imgsize'] = int(imgsize)
        archfiles_row['slot'] = int(re.search(
                r'acaf\d+N\d{3}_(\d)_img0.fits(\.gz)?',
                filename).group(1))
        archfiles_row['filename'] = filename
        archfiles_row['filetime'] = int(
            re.search(r'(\d+)', archfiles_row['filename']).group(1))
        filedate = DateTime(archfiles_row['filetime']).date
        year, doy = (int(x) for x in
                     re.search(r'(\d\d\d\d):(\d\d\d)', filedate).groups())
        archfiles_row['year'] = year
        archfiles_row['doy'] = doy
        archfiles_row['rows'] = len(hdu.data)
        hdus.close()

        # remove old versions of this file
        oldmatches = db.fetchall(
            """SELECT * from archfiles
               WHERE filetime = %(filetime)d
               and slot = %(slot)d
               and startmjf = %(startmjf)d and startmnf = %(startmnf)d
               and stopmjf = %(stopmjf)d and stopmnf = %(stopmnf)d """
            % archfiles_row)
        if len(oldmatches):
            self._arch_remove(oldmatches)

        interval_matches = get_interval_files(archfiles_row['tstart'],
                                              archfiles_row['tstop'],
                                              slots=[archfiles_row['slot']])
        if len(interval_matches):
            # if there are files there that still overlap the new file
            # and they are all older, remove the old files.
            if np.all(interval_matches['revision']
                      < archfiles_row['revision']):
                logger.info(
                    "removing overlapping files at older revision(s)")
                logger.info(interval_matches)
                self._arch_remove(interval_matches)
            # if the overlapping files are all from the same revision
            # just ingest them and hope for the best
            elif np.all(interval_matches['revision']
                        == archfiles_row['revision']):
                logger.info("ignoring overlap for same process revision")
            # if the files that overlap are all newer, let's not ingest
            # the "missing" file
            elif np.all(interval_matches['revision']
                        > archfiles_row['revision']):
                return None
            else:
                logger.error(archfiles_row)
                logger.error(interval_matches)
                # throw an error if there is still overlap
                raise ValueError("Cannot ingest %s, overlaps existing files"
                                 % filename)

        return archfiles_row

    def _arch_remove(self, defunct_matches):
        for file_record in defunct_matches:
            query = ("""delete from archfiles
                          WHERE filetime = %(filetime)d
                          and slot = %(slot)d
                          and startmjf = %(startmjf)d
                          and startmnf = %(startmnf)d
                          and stopmjf = %(stopmjf)d
                          and stopmnf = %(stopmnf)d """
                       % file_record)
            logger.info(query)
            self.db.execute(query)
            self.db.commit()
            archdir = os.path.abspath(os.path.join(
                    self.data_root,
                    str(file_record['year']),
                    "{0:03d}".format(file_record['doy'])
                    ))
            logger.info("deleting %s" %
                        os.path.join(archdir, file_record['filename']))
            real_file = os.path.join(archdir, file_record['filename'])
            if os.path.exists(real_file):
                #shutil.move(real_file, '/data/aca/archive/aca0/deleted')
                os.unlink(real_file)

    def _move_archive_files(self, archfiles):
        """
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
            tstart = re.search(r'(\d+)', str(basename)).group(1)
            datestart = DateTime(tstart).date
            year, doy = re.search(r'(\d\d\d\d):(\d\d\d)', datestart).groups()
            archdir = os.path.abspath(os.path.join(data_root,
                                                   year,
                                                   doy))
            # construct the destination filepath/name
            archfile = os.path.abspath(os.path.join(archdir, basename))
            if not os.path.exists(archdir):
                os.makedirs(archdir)
            if not os.path.exists(archfile):
                logger.info('mv %s %s' % (os.path.abspath(f), archfile))
                os.chmod(f, 0775)
                shutil.move(f, archfile)
            if os.path.exists(f):
                logger.info('Unlinking %s' % os.path.abspath(f))
                os.unlink(f)

    def _fetch_by_time(self, range_tstart, range_tstop):
        logger.info("Fetching %s from %s to %s" 
                    % ('ACA L0 Data',
                       DateTime(range_tstart).date,
                       DateTime(range_tstop).date))
        archfiles = self._get_archive_files(DateTime(range_tstart),
                                            DateTime(range_tstop))
        return archfiles

    def _fetch_individual_files(self, files):
        arc5 = Ska.arc5gl.Arc5gl(echo=True)
        logger.info('********** %s %s **********'
                    % (filetype['content'], time.ctime()))
        fetched_files = []
        ingest_dates = []
        # get the files, store in file archive, and record in database
        for file in files:
            # Retrieve CXC archive files in a temp directory with arc5gl
            missed_file = file['filename']
            arc5.sendline('dataset=flight')
            arc5.sendline('detector=pcad')
            arc5.sendline('subdetector=aca')
            arc5.sendline('level=0')
            arc5.sendline('filename=%s' % missed_file)
            arc5.sendline('version=last')
            arc5.sendline('operation=retrieve')
            arc5.sendline('go')
            have_files = sorted(glob(filetype['fileglob']))
            if not len(have_files):
                raise ValueError
            filename = have_files[0]
            # if it isn't gzipped, just gzip it
            if re.match('.*\.fits$', filename):
                import gzip
                f_in = open(file, 'rb')
                f_out = gzip.open("%s.gz" % filename, 'wb')
                f_out.writelines(f_in)
                f_out.close()
                f_in.close()
                filename = "%s.gz" % filename
            fetched_files.append(filename)
            ingest_dates.append(file['ingest_date'])
        return fetched_files, ingest_dates

    def _insert_files(self, files):
        for i, f in enumerate(files):
            arch_info = self._read_archfile(i, f, files)
            if arch_info:
                self._move_archive_files([f])
                self.db.insert(arch_info, 'archfiles')
        self.db.commit()

    def update(self):
        """
        Retrieve ACA0 telemetry files from the CXC archive, store in the
        Ska/ACA archive, and update database of files.
        """
        contentdir = self.data_root
        if not os.path.exists(contentdir):
            os.makedirs(contentdir)
        archdb = os.path.join(contentdir, 'archfiles.db3')
        # if the database of the archived files does not exist,
        # make it
        if not os.path.exists(archdb):
            logger.info("creating archfiles db from %s"
                        % self.sql_def)
            db_sql = os.path.join(os.environ['SKA_DATA'],
                                  'mica', self.sql_def)
            db_init_cmds = file(db_sql).read()
            self.db.execute(db_init_cmds, commit=True)
        if self.start:
            datestart = DateTime(self.start)
        else:
            # Get datestart as the most-recent file time from archfiles table
            # will need min-of-max-slot-datestart
            last_time = min([self.db.fetchone(
                        "select max(filetime) from archfiles where slot = %d"
                        % s)['max(filetime)'] for s in range(0, 8)])
            datestart = DateTime(last_time)
        datestop = DateTime(self.stop)
        padding_seconds = 10000
        # loop over the specified time range in chunks of
        # days_at_once in seconds with some padding
        for tstart in np.arange(datestart.day_start().secs,
                                datestop.day_end().secs,
                                self.days_at_once * 86400):
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

        timestamp_file = os.path.join(self.data_root, 'last_timestamp.txt')
        cda_checked_timestamp = open(timestamp_file).read().rstrip()
        missing_datetime = DateTime(cda_checked_timestamp)
        missing_files, last_ingest_date = \
            self._get_missing_archive_files(missing_datetime,
                                            only_new=True)

        # update the file to have up through the last confirmed good file
        # even before we try to fetch missing ones
        open(timestamp_file, 'w').write("%s" % last_ingest_date)

        if len(missing_files):
            logger.info("Found %d missing individual files"
                        % len(missing_files))
            # make a temporary directory
            tmpdir = Ska.File.TempDir(dir=self.temp_root)
            dirname = tmpdir.name
            logger.info("File save to temp dir %s" % dirname)
            with Ska.File.chdir(dirname):
                fetched_files, ingest_times = \
                    self._fetch_individual_files(missing_files)
                self._insert_files(fetched_files)

            last_ingest_date = missing_files[-1]['ingest_date']
            # update the file to have up through the last confirmed good file
            # even before we try to fetch missing ones
            open(timestamp_file, 'w').write("%s" % last_ingest_date)


def main():
    """
    Command line interface to fetch ACA L0 telemetry from the CXC Archive
    and store it in the Ska archive.
    """
    opt = get_options()
    kwargs = vars(opt)
    updater = Updater(**kwargs)
    updater.update()

if __name__ == '__main__':
    main()
