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

aca_dtype = [('TIME', '>f8'), ('QUALITY', '>i4'), ('MJF', '>i4'),
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
             ('IMGSIZE', '>i4'),
             ]

aca_dtype_names = [k[0] for k in aca_dtype]

config = dict(data_root='/data/aca/archive/aca0',
              temp_root='/data/aca/archive/temp',
              days_at_once=30.0,
              sql_def='archfiles_aca_l0_def.sql')


def get_options():
    parser = argparse.ArgumentParser(
        description="Fetch aca level 0 products and make a file archive")
    defaults = dict(config)
    parser.set_defaults(**defaults)
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--temp-root",
                        help="parent temp directory")
    parser.add_argument("--datestart",
                        help="start date for retrieve (defaults to max date of archived files)")
    parser.add_argument("--datestop",
                        help="stop date for retrieve (defaults to now)")
    parser.add_argument("--days-at-once",
                        type=float,
                        help="if over this number, bin to chunks of this number of days")
    opt = parser.parse_args()
    return opt


def get_slot_data(tstart, tstop, slot, imgsize=[4, 6, 8],
                  db=None, data_root=config['data_root'],
	              columns=aca_dtype_names):
    if not db:
        dbfile = os.path.join(config['data_root'], 'archfiles.db3')
        db = Ska.DBI.DBI(dbi='sqlite', server=dbfile)
    data_files = get_interval_files(tstart, tstop, slots=[slot],
                                    imgsize=imgsize, db=db)
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
    return all_rows


def get_interval_files(tstart, tstop,
                       slots=[0, 1, 2, 3, 4, 5, 6, 7], 
                       imgsize=[4, 6, 8], db=None):
    if not db:
        dbfile = os.path.join(config['data_root'], 'archfiles.db3')
        db = Ska.DBI.DBI(dbi='sqlite', server=dbfile)
    tstart = DateTime(tstart).secs
    tstop = DateTime(tstop).secs
    imgsize_str = ','.join([str(x) for x in imgsize])
    slot_str = ','.join([str(x) for x in slots])
    # a composite index isn't as fast as just doing a padded search on one
    # index first (tstart).  This gets extra files to make sure we don't miss the 
    # case where tstart is in the middle of an interval, but drastically
    # reduces the size of the bsearch on the tstop index
    tstart_pad = 10 * 86400
    db_query = ('SELECT * FROM archfiles '
                'WHERE tstart >= %f - %f '
                'AND tstart < %f '
                'AND tstop > %f '
                'AND slot in (%s) '
                'AND imgsize in (%s) '
                'order by filetime asc '
                % (tstart, tstart_pad, tstop, tstart, slot_str, imgsize_str))
    return db.fetchall(db_query)




def rebuild_database(db=None, db_file=None):
    if db is None and db_file is None:
        raise ValueError
    if db is None and db_file is not None:
        logger.info("creating archfiles db from %s"
                    % config['sql_def'])
        db_sql = os.path.join(os.environ['SKA_DATA'],
                              'mica', config['sql_def'])
        db_init_cmds = file(db_sql).read()
        db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
                         autocommit=False)
        db.execute(db_init_cmds, commit=True)
    year_dirs = sorted(glob(
            os.path.join(config['data_root'], '[12][0-9][0-9][0-9]')))
    for ydir in year_dirs:
        day_dirs = sorted(glob(
                os.path.join(ydir, '[0-3][0-9][0-9]')))
        for ddir in day_dirs:
            archfiles = sorted(glob(
                    os.path.join(ddir, '*_img0*')))
            db.execute("begin transaction")
            for i, f in enumerate(archfiles):
                arch_info = read_archfile(i, f, archfiles, db)
                if arch_info:
                    db.insert(arch_info, 'archfiles')
            db.commit()
                
                         

def get_archive_files(filetype, datestart, datestop):
    """Update FITS file archive with arc5gl and ingest files into msid (HDF5) archive"""

    # Retrieve CXC archive files in a temp directory with arc5gl
    arc5 = Ska.arc5gl.Arc5gl(echo=True)
    logger.info('********** %s %s **********' 
                % (filetype['content'], time.ctime()))

    arc5.sendline('tstart=%s' % DateTime(datestart).date)
    arc5.sendline('tstop=%s' % DateTime(datestop).date)
    arc5.sendline('get %s' % filetype['arc5gl_query'].lower())

    return sorted(glob(filetype['fileglob']))


def read_archfile(i, f, archfiles, db):
    """Read filename ``f`` with index ``i`` (position within list of filenames).  The
    file has type ``filetype`` and will be added to MSID file at row index ``row``.
    ``colnames`` is the list of column names for the content type (not used here).
    """
    # Check if filename is already in archfiles.
    # If so then abort further processing.
    filename = os.path.basename(f)
    if db.fetchall('SELECT filename FROM archfiles WHERE filename=?',
                   (filename,)):
        logger.info(
            'File %s already in archfiles - unlinking and skipping' % f)
        os.unlink(f)
        return None

    # Read FITS archive file and accumulate data into dats list and header into headers dict
    logger.info('Reading (%d / %d) %s' % (i, len(archfiles), filename))
    hdus = pyfits.open(f)
    hdu = hdus[1]

    # Accumlate relevant info about archfile that will be ingested into
    # MSID h5 files.  Commit info before h5 ingest so if there is a failure
    # the needed info will be available to do the repair.
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
    year, doy = (int(x)
                 for x in re.search(r'(\d\d\d\d):(\d\d\d)', filedate).groups())
    archfiles_row['year'] = year
    archfiles_row['doy'] = doy
    archfiles_row['rows'] = len(hdu.data)
    hdus.close()
    return  archfiles_row


def move_archive_files(filetype, archfiles, opt):

    if not os.path.exists(opt.data_root):
        os.makedirs(opt.data_root)

    for f in archfiles:
        if not os.path.exists(f):
            continue
        basename = os.path.basename(f)
        tstart = re.search(r'(\d+)', str(basename)).group(1)
        datestart = DateTime(tstart).date
        year, doy = re.search(r'(\d\d\d\d):(\d\d\d)', datestart).groups()

        archdir = os.path.abspath(os.path.join(opt.data_root,
                                               year,
                                               doy))
        archfile = os.path.abspath(os.path.join(archdir, basename))
        if not os.path.exists(archdir):
            os.makedirs(archdir)

        if not os.path.exists(archfile):
            logger.info('mv %s %s' % (os.path.abspath(f), archfile))
            shutil.move(f, archfile)

        if os.path.exists(f):
            logger.info('Unlinking %s' % os.path.abspath(f))
            os.unlink(f)




def update_archive(opt):
    contentdir = os.path.join(opt.data_root)
    if not os.path.exists(contentdir):
        os.makedirs(contentdir)
    archdb = os.path.join(contentdir, 'archfiles.db3')
    if not os.path.exists(archdb):
        logger.info("creating archfiles db from %s"
                    % config['sql_def'])
        db_sql = os.path.join(os.environ['SKA_DATA'],
                              'mica', config['sql_def'])
        db_init_cmds = file(db_sql).read()
        db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
                         autocommit=False)
        db.execute(db_init_cmds, commit=True)
    db = Ska.DBI.DBI(dbi='sqlite', server=archdb, autocommit=False)
    if opt.datestart:
        datestart = DateTime(opt.datestart)
    else:
        # Get datestart as the most-recent file time from archfiles table
        # will need min-of-max-slot-datestart
        last_time = min([db.fetchone(
                    "select max(filetime) from archfiles where slot = %d"
                    % s)['max(filetime)'] for s in range(0, 8)])
        datestart = DateTime(last_time)
    datestop = DateTime(opt.datestop)
    padding_seconds = 10000
    for tstart in np.arange(datestart.day_start().secs,
                            datestop.day_end().secs,
                            opt.days_at_once * 86400):

        # set times for a chunk
        range_tstart = tstart - padding_seconds
        range_tstop = tstart + opt.days_at_once * 86400
        if range_tstop > datestop.day_end().secs:
            range_tstop = datestop.day_end().secs
        range_tstop += padding_seconds
        print DateTime(range_tstart).date + " " + DateTime(range_tstop).date
        # make a temporary directory
        tmpdir = Ska.File.TempDir(dir=opt.temp_root)
        dirname = tmpdir.name
        print dirname
        # get the files
        with Ska.File.chdir(dirname):
            archfiles = get_archive_files(filetype,
                                          DateTime(range_tstart),
                                          DateTime(range_tstop))
            for i, f in enumerate(archfiles):
                arch_info = read_archfile(i, f, archfiles, db)
                if arch_info:
                    move_archive_files(filetype, [f], opt)
                    db.insert(arch_info, 'archfiles')
        db.commit()


def main():
    opt = get_options()
    update_archive(opt)


if __name__ == '__main__':
    main()
