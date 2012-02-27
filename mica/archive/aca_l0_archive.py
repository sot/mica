import os
import Ska.arc5gl
#import tempfile
from glob import glob
import re
import logging
import shutil
import Ska.DBI
#from Ska.Shell import bash
import numpy as np
from Chandra.Time import DateTime
import Ska.File
import time

arc5 = Ska.arc5gl.Arc5gl()
logger = logging.getLogger('asp1 fetch')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

# borrowed from eng_archive
archfiles_hdr_cols = ('tstart', 'tstop', 'startmjf', 'startmnf',
                      'stopmjf', 'stopmnf',
                      'tlmver', 'ascdsver', 'revision', 'date',
                      'imgsize')

#import itertools
import pyfits

dtype = [('level', '|S4'), ('instrum', '|S7'), ('content', '|S13'),
         ('arc5gl_query', '|S27'), ('fileglob', '|S9')]
filetypes = np.rec.fromrecords([('L0', 'PCAD', 'ACADATA', 'ACA0', '*fits.gz')],
                               dtype=dtype)


def get_archive_files(filetype, datestart, datestop):
#    """Update FITS file archive with arc5gl and ingest files into msid (HDF5) archive"""
#    
    # If running on the OCC GRETA network the cwd is a staging directory that
    # could already have files.  Also used in testing.
    # Don't return more than opt.max_arch_files files at once because of memory
    # issues on gretasot.  This only comes up when there has been some problem or stoppage.
    files = sorted(glob(filetype['fileglob']))
    #if opt.occ or files:
    #    return sorted(files)[:opt.max_arch_files]

    # Retrieve CXC archive files in a temp directory with arc5gl
    arc5 = Ska.arc5gl.Arc5gl(echo=True)

    # For *ephem0 the query needs to extend well into the future
    # to guarantee getting all available files.  This is the archives fault.
    if filetype['level'] == 'L0' and filetype['instrum'] == 'EPHEM':
        datestop = datestop + 50

    # For instrum==EPHEM break queries into time ranges no longer than
    # 100000 sec each.  EPHEM files are at least 7 days long and generated
    # no more often than every ~3 days so this should work.
    n_queries = (1 if filetype['instrum'] != 'EPHEM'
          else 1 + round((datestop.secs - datestart.secs) / 100000.))
    times = np.linspace(datestart.secs, datestop.secs, n_queries + 1)

    logger.info('********** %s %s **********' % (filetype['content'], time.ctime()))

    for t0, t1 in zip(times[:-1], times[1:]):
        arc5.sendline('tstart=%s' % DateTime(t0).date)
        arc5.sendline('tstop=%s' % DateTime(t1).date)
        arc5.sendline('get %s' % filetype['arc5gl_query'].lower())

    return sorted(glob(filetype['fileglob']))

def read_archfile(i, f, archfiles, db):
    """Read filename ``f`` with index ``i`` (position within list of filenames).  The
    file has type ``filetype`` and will be added to MSID file at row index ``row``.
    ``colnames`` is the list of column names for the content type (not used here).
    """
    # Check if filename is already in archfiles.  If so then abort further processing.
    filename = os.path.basename(f)
    colnames = archfiles_hdr_cols
    if db.fetchall('SELECT filename FROM archfiles WHERE filename=?', (filename,)):
        logger.info('File %s already in archfiles - unlinking and skipping' % f)
        os.unlink(f)
        return None

    # Read FITS archive file and accumulate data into dats list and header into headers dict
    logger.info('Reading (%d / %d) %s' % (i, len(archfiles), filename))
    hdus = pyfits.open(f)
    hdu = hdus[1]
    #dat = converters.convert(hdu.data, filetype['content'])

    # Accumlate relevant info about archfile that will be ingested into
    # MSID h5 files.  Commit info before h5 ingest so if there is a failure
    # the needed info will be available to do the repair.
    archfiles_row = dict((x, hdu.header.get(x.upper())) for x in archfiles_hdr_cols)
    archfiles_row['checksum'] = hdu._checksum
    #archfiles_row['rowstart'] = row
    #archfiles_row['rowstop'] = row + len(dat)
    imgsize_sq = hdu.data[0].field('IMGRAW').shape[0]
    archfiles_row['imgsize' ] = int(np.sqrt(imgsize_sq))
    archfiles_row['slot'] = int(re.search(r'acaf\d+N\d{3}_(\d)_img0.fits(\.gz)?', 
                                          filename).group(1))
    archfiles_row['filename'] = filename
    archfiles_row['filetime'] = int(re.search(r'(\d+)', archfiles_row['filename']).group(1))
    filedate = DateTime(archfiles_row['filetime']).date
    year, doy = (int(x) for x in re.search(r'(\d\d\d\d):(\d\d\d)', filedate).groups())
    archfiles_row['year'] = year
    archfiles_row['doy'] = doy
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
                                               filetype['content'].lower(),
                                               year,
                                               doy))
        archfile = os.path.abspath(os.path.join(archdir, basename))

        if not os.path.exists(archdir):
            os.makedirs(archdir)

        if not os.path.exists(archfile):
            logger.info('mv %s %s' % (os.path.abspath(f), archfile))
#            if not opt.dry_run:
#                if not opt.occ:
#            shutil.copy2(f, stagedir)
            shutil.move(f, archfile)

        if os.path.exists(f):
            logger.info('Unlinking %s' % os.path.abspath(f))
            os.unlink(f)




def get_arch_info(i, f, archfiles, db):
    """Read filename ``f`` with index ``i`` (position within list of filenames).  The
    file has type ``filetype`` and will be added to MSID file at row index ``row``.
    ``colnames`` is the list of column names for the content type (not used here).
    """

    filename = os.path.basename(f)

    # Read FITS archive file and accumulate data into dats list and header into headers dict
    logger.debug('Reading (%d / %d) %s' % (i, len(archfiles), filename))
    hdus = pyfits.open(f)
    hdu = hdus[1]

    # Accumlate relevant info about archfile that will be ingested into
    # MSID h5 files.  Commit info before h5 ingest so if there is a failure
    # the needed info will be available to do the repair.
    archfiles_row = dict((x, hdu.header.get(x.upper())) for x in archfiles_hdr_cols)
    archfiles_row['checksum'] = hdu._checksum
    archfiles_row['filename'] = filename
    archfiles_row['filetime'] = int(re.search(r'(\d+)', archfiles_row['filename']).group(1))
    filedate = DateTime(archfiles_row['filetime']).date
    year, doy = (int(x) for x in re.search(r'(\d\d\d\d):(\d\d\d)', filedate).groups())
    archfiles_row['year'] = year
    archfiles_row['doy'] = doy
    
    hdus.close()
    return archfiles_row




def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--data-root",
                      default="data")
    parser.add_option("--temp-root",
                      default="temp")
    opt, args = parser.parse_args()
    return opt, args


def main(opt):

    tmpdir = Ska.File.TempDir(dir=opt.temp_root)
    dirname = tmpdir.name

    filetype = filetypes[0]
    contentdir = os.path.join(opt.data_root, filetype['content'].lower())
    if not os.path.exists(contentdir):
        os.makedirs(contentdir)
    archdb = os.path.join(contentdir, 'archfiles.db3')

    # Get datestart as the most-recent file time from archfiles table
    db = Ska.DBI.DBI(dbi='sqlite', server=archdb, autocommit=False)
    # will need min-of-max-slot-datestart
    last_time = min([db.fetchone("select max(filetime) from archfiles where slot = %d" 
                                 % s)['max(filetime)'] for s in range(0, 8)])
    if last_time:
        datestart = DateTime(last_time)
    else:
        datestart = DateTime(-2)
    datestop = DateTime()

    print dirname
    with Ska.File.chdir(dirname):
        archfiles = get_archive_files(filetype, datestart, datestop)
        for i, f in enumerate(archfiles):
            arch_info = read_archfile(i, f, archfiles, db)
            if arch_info:
                move_archive_files(filetype, [f], opt)
                db.insert(arch_info, 'archfiles')
    db.commit()

if __name__ == '__main__':
    opt, args = get_options()
    main(opt)
