import cPickle as pickle
import numpy as np

dtype = [('level', '|S4'), ('instrum', '|S7'), ('content', '|S13'),
         ('arc5gl_query', '|S27'), ('fileglob', '|S9')]
filetypes = np.rec.fromrecords([('L0', 'PCAD', 'ACADATA', 'ACA0', '*fits.gz')],
                               dtype=dtype)


import ConfigParser
proj_config = ConfigParser.SafeConfigParser()
proj_config.read('archive.conf')

import os
import pyyaks.context
#from . import file_defs
from Ska.engarchive import units
import pyfits
import time

file_defs_msid_files = {'archfiles': '{{ft.content}}/archfiles.db3',
             'colnames': '{{ft.content}}/colnames.pickle',
             'colnames_all': '{{ft.content}}/colnames_all.pickle',
             'contentdir': '{{ft.content}}/',
             'data': '{{ft.content}}/{{ft.msid | upper}}.h5',
             'filetypes': 'filetypes.dat',
             'headers': '{{ft.content}}/headers.pickle',
             'msid': '{{ft.content}}/{{ft.msid | upper}}.h5',
             'msid_bad_times': 'msid_bad_times.dat',
             'stats': '{{ft.content}}/{{ft.interval}}/{{ft.msid | upper}}.h5',
             'statsdir': '{{ft.content}}/{{ft.interval}}/'}


SKA = os.getenv('SKA') or '/proj/sot/ska'
#ENG_ARCHIVE = os.getenv('ENG_ARCHIVE') or SKA + '/data/aca/'
ENG_ARCHIVE = '/data/aca/archive'
IGNORE_COLNAMES = ('TIME', 'MJF', 'MNF', 'TLM_FMT')
DIR_PATH = os.path.dirname(os.path.abspath(__file__))

# Maximum number of MSIDs that should ever match an input MSID spec
# (to prevent accidentally selecting a very large number of MSIDs)
MAX_GLOB_MATCHES = 10

## Global (eng_archive) definition of file names
msid_files = pyyaks.context.ContextDict('msid_files', basedir=ENG_ARCHIVE)
msid_files.update(file_defs_msid_files)


# Context dictionary to provide context for msid_files
ft = pyyaks.context.ContextDict('ft')

import collections
content = collections.OrderedDict()
for filetype in filetypes:
    ft['content'] = filetype['content'].lower()
    try:
        colnames = pickle.load(open(msid_files['colnames'].abs))
        content.update((x, ft['content'].val) for x in sorted(colnames)
                       if x not in IGNORE_COLNAMES)
    except IOError:
        pass


def msid_glob(msid):
    """Get the archive MSIDs matching ``msid``.

    The function returns a tuple of (msids, MSIDs) where ``msids`` is a list of
    MSIDs that is all lower case and (where possible) matches the input
    ``msid``.  The output ``MSIDs`` is all upper case and corresponds to the
    exact MSID names stored in the archive HDF5 files.

    :param msid: input MSID glob
    :returns: tuple (msids, MSIDs)
    """

    MSID = msid.upper()
    # First try MSID or DP_<MSID>.  If success then return the upper
    # case version and whatever the user supplied (could be any case).
    for match in (MSID, 'DP_' + MSID):
        if match in content:
            return [msid], [match]

    # Next try as a file glob.  If there is a match then return a
    # list of matches, all lower case and all upper case.  Since the
    # input was a glob the returned msids are just lower case versions
    # of the matched upper case MSIDs.
    for match in (MSID, 'DP_' + MSID):
        matches = fnmatch.filter(content.keys(), match)
        if matches:
            if len(matches) > MAX_GLOB_MATCHES:
                raise ValueError(
                    'MSID spec {} matches more than {} MSIDs.  '
                    'Refine the spec or increase fetch.MAX_GLOB_MATCHES'
                    .format(msid, MAX_GLOB_MATCHES))
            return [x.lower() for x in matches], matches

    raise ValueError('MSID {} is not in Archive'.format(MSID))

from Chandra.Time import DateTime

import Ska.engarchive.fetch

class SlotMSID(Ska.engarchive.fetch.MSID):
    
    def __init__(self, msid, slot, start, stop=None, filter_bad=False, stat=None):
        msids, MSIDs = msid_glob(msid)
        if len(MSIDs) > 1:
            raise ValueError('Multiple matches for {} in Eng Archive'
                             .format(msid))
        else:
            self.msid = msids[0]
            self.MSID = MSIDs[0]

        self.slot = slot
        self.unit = units.get_msid_unit(self.MSID)
        self.stat = stat
        if stat:
            self.dt = {'5min': 328, 'daily': 86400}[stat]

        self.tstart = DateTime(start).secs
        self.tstop = (DateTime(stop).secs if stop else
                      DateTime(time.time(), format='unix').secs)
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        try:
            self.content = content[self.MSID]
        except KeyError:
            raise ValueError('MSID %s is not in Eng Archive' % self.MSID)

        # Get the times, values, bad values mask from the HDF  files archive
        self._get_data()
#
#        # If requested filter out bad values and set self.bad = None
#        if filter_bad:
#            self.filter_bad()

    def _get_data(self):
        """Get data from the Eng archive"""
#        logger.info('Getting data for %s between %s to %s',
#                    self.msid, self.datestart, self.datestop)
        with Ska.engarchive.fetch._cache_ft():
            ft['content'] = self.content
            ft['msid'] = self.MSID

            gmd = self._get_slot_msid_data
            self.vals, self.times, self.bads, self.colnames = \
                    gmd(self.content, self.tstart, self.tstop, self.MSID, self.slot)

    @staticmethod
    def _get_slot_msid_data(content, tstart, tstop, msid, slot):
        acadata = get_slot_aca_data(content, tstart, tstop, slot)
        return acadata[msid], acadata['TIME'], \
            np.isnan(acadata[msid]), ['val', 'times', 'bads']


class memoized(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            self.cache[args] = value = self.func(*args)
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)
    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    
aca_dtype = [('TIME', '>f8'), ('QUALITY', '>i4'), ('MJF', '>i4'), ('MNF', '>i4'),
             ('END_INTEG_TIME', '>f8'), ('INTEG', '>f4'), ('GLBSTAT', '|u1'),
             ('COMMCNT', '|u1'), ('COMMPROG', '|u1'), ('IMGFID1', '|u1'),
             ('IMGNUM1', '|u1'), ('IMGFUNC1', '|u1'), ('IMGSTAT', '|u1'),
             ('IMGROW0', '>i2'), ('IMGCOL0', '>i2'), ('IMGSCALE', '>i2'),
             ('BGDAVG', '>i2'), ('IMGFID2', '|u1'), ('IMGNUM2', '|u1'),
             ('IMGFUNC2', '|u1'), ('BGDRMS', '>i2'), ('TEMPCCD', '>f4'),
             ('TEMPHOUS', '>f4'), ('TEMPPRIM', '>f4'), ('TEMPSEC', '>f4'),
             ('BGDSTAT', '|u1'), ('IMGFID3', '|u1'), ('IMGNUM3', '|u1'),
             ('IMGFUNC3', '|u1'), ('IMGFID4', '|u1'), ('IMGNUM4', '|u1'),
             ('IMGFUNC4', '|u1'), ('IMGRAW', '>f4', (64,)), ('HD3TLM62', '|u1'),
             ('HD3TLM63', '|u1'), ('HD3TLM64', '|u1'), ('HD3TLM65', '|u1'),
             ('HD3TLM66', '|u1'), ('HD3TLM67', '|u1'), ('HD3TLM72', '|u1'),
             ('HD3TLM73', '|u1'), ('HD3TLM74', '|u1'), ('HD3TLM75', '|u1'),
             ('HD3TLM76', '|u1'), ('HD3TLM77', '|u1')]



        

#@memoized
def get_slot_aca_data(content, tstart, tstop, slot):
    data_files = get_interval_files(content, tstart, tstop, slot)
    rows = np.sum(data_files['rows'])
    zero_row = np.zeros(1, dtype=aca_dtype)
    all_rows = zero_row.repeat(rows)
    all_rows['IMGRAW'].reshape(rows,8,8)[:] = np.nan
    rowcount = 0
    for f in data_files:
        fp = os.path.join(ENG_ARCHIVE, content,
                          str(f['year']), str(f['doy']), f['filename'])
        print fp
        hdu = pyfits.open(fp)
        chunk = hdu[1].data
        for fname in all_rows.dtype.names:
            if fname == 'IMGRAW':
                continue
            if fname in chunk.dtype.names:
                all_rows[fname][rowcount:rowcount+len(chunk)] = chunk.field(fname)
            else:
                all_rows[fname][rowcount:rowcount+len(chunk)] = np.nan
        imgsize = int(np.sqrt(chunk[0].field('IMGRAW').size))
        all_rows['IMGRAW'].reshape(rows,8,8)[
            rowcount:rowcount+len(chunk), 0:imgsize, 0:imgsize] = (
            chunk.field('IMGRAW').reshape(len(chunk), imgsize, imgsize))
        rowcount += len(chunk)
    return all_rows
    

#@memoized
def get_interval_files(content, tstart, tstop, slot, imgsize=[4,6,8]):

    import Ska.DBI
    ft['content'] = content
    db = Ska.DBI.DBI(dbi='sqlite', server=os.path.join(ENG_ARCHIVE, 
                                                       'acadata', 
                                                       'archfiles.db3'))
    imgsize_str = ','.join([str(x) for x in imgsize])
    db_query = ('SELECT * FROM archfiles '
                'WHERE tstop > %f '
                'AND tstart < %f '
                'AND slot == %d '
                'AND imgsize in (%s) '
                'order by filetime asc ' 
                % (tstart,tstop, slot, imgsize_str))
    time_list = db.fetchall(db_query)
    if not len(time_list):
        raise ValueError
    return time_list
        
