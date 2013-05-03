"""
Generalized module for fetching and archiving obsid-organized
telemetry such as asp_l1 and obspar products.
"""
import os
import tempfile
from glob import glob
import re
import logging
import shutil
import numpy as np
import pyfits

import Ska.arc5gl
import Ska.DBI
from Chandra.Time import DateTime
import Ska.File
import mica.version as mica_version

# borrowed from telem_archive
import csv
import gzip


def parse_obspar(file):
    """
    Return the rows of an IRAF formatted obspar as a dictionary.

    :param file: obspar file

    :returns: row of obspar
    :rtype: dictionary generator
    """
    convert = {'i': int,
               'r': float,
               's': str}
    try:
        lines = gzip.open(file).readlines()
    except IOError:
        lines = open(file).readlines()
    obs_read = csv.DictReader(lines,
                              fieldnames=('name', 'type', 'hidden', 'value',
                                          'def1', 'def2', 'descr'),
                              dialect='excel')

    for row in obs_read:
        # this empty-string '' hack is not present in the original
        if ((row['value'] == '')
                and ((row['type'] == 'r') or (row['type'] == 'i'))):
            row['value'] = None
        else:
            row['value'] = convert[row['type']](row['value'])
        row['name'] = row['name'].replace('-', '_')
        yield row

    return


def get_obspar(obsparfile):
    """Get the obspar for obsid starting at tstart.  Return as a dict."""
    obspar = dict()
    for row in parse_obspar(obsparfile):
        obspar.update({row['name']: row['value']})
    return obspar


class ProductVersionError(Exception):
    pass


class ObsArchive:
    """
    Object to store configuration, logging, and processing tasks
    to fetch obsid telemetry from the CXC archive and store in a Ska
    file archive, while logging the archive files to a file lookup
    database.

    :param config: configuration dictionary
    :returns: ObsArchive instance

    The configuration dictionary may have these key/values:

    * data_root: directory for products
                (example /data/aca/archive/asp1)
    * temp_root: directory for temporary storage of fetched
                telemetry
    * cols: headers that will be included in file lookup table
    * sql_def: sql file to build file lookup archfiles table
    * apstat_table: axafapstat database table from which to find
                   new processing (by id)
    * apstat_id: field in apstat_table to use as unique CXCDS
                processing id
    * label: label of product type for log messages
    * small: arc5gl keyword/filetype for small file from products
            (example asp1{fidprops}).  This will be retrieved with
            "get %s" % config['small'] and the retrieved files will
            be used to determine product version.
    * small_glob: glob to match files retrieved by
                 "get %s" % config[small]
                 (example '*fidpr*')
    * small_ver_regex: regular expression to search for version from
                     retrieved files (example 'pacdf\d+N(\d{3})_')
    * full: arc5gl keyword for products (example 'asp1')
    * rebuild: If True/set, allow update mode to rebuild the database
               from obsid 1.
    """

    def __init__(self, config):
        self._config = config
        self._logger = logging.getLogger('ObsArchive')

    @property
    def config(self):
        return self._config

    @property
    def logger(self):
        return self._logger

    def set_env(self):
        """
        Set environment included an arc5gl handle and
        and a handle to the axafapstat database
        """
        self._arc5 = Ska.arc5gl.Arc5gl()
        self._apstat = Ska.DBI.DBI(dbi='sybase', server='sqlsao',
                                   database='axafapstat')
        self._aca_db = Ska.DBI.DBI(dbi='sybase', server='sybase',
                                   user='aca_read')
        config = self.config
        db_file = os.path.join(os.path.abspath(config['data_root']),
                               'archfiles.db3')
        if not os.path.exists(db_file) or os.stat(db_file).st_size == 0:
            if not os.path.exists(config['data_root']):
                os.makedirs(config['data_root'])
            self.logger.info("creating archfiles db from %s"
                             % config['sql_def'])
            db_sql = os.path.join(os.environ['SKA_DATA'],
                                  'mica', config['sql_def'])
            db_init_cmds = file(db_sql).read()
            db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
                             autocommit=False)
            db.execute(db_init_cmds, commit=True)
        db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
                         autocommit=False)
        self._archfiles_db = db

    def set_read_env(self):
        """
        Set environment included an arc5gl handle and
        and a handle to the axafapstat database
        """
        config = self.config
        db_file = os.path.join(os.path.abspath(config['data_root']),
                               'archfiles.db3')
        db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
                         autocommit=False)
        self._archfiles_db = db


    def get_all_obspar_info(self, i, f, archfiles):
        """
        Read obspar and add 'obsid' and 'filename' keys to the dictionary
        i and archfiles are just passed to make the logging prettier.
        """
        logger = self.logger
        filename = os.path.basename(f)
        logger.debug('Reading (%d / %d) %s' % (i, len(archfiles), filename))
        obspar = get_obspar(f)
        obspar['obsid'] = obspar['obs_id']
        obspar['filename'] = filename
        return obspar

    def get_files(self, obsid=None, start=None, stop=None,
                  revision=None, content=None):
        data_root = self.config['data_root']
        if obsid is None:
            if start is None or stop is None:
                raise TypeError("Must supply either obsid or start and stop")
        file_records = self._get_file_records(obsid=obsid,
                                              start=start, stop=stop,
                                              revision=revision,
                                              content=content)
        files = [os.path.join(data_root,
                              ("%05d" % f['obsid'])[0:2],
                              "%05d_v%02d" % (f['obsid'], f['revision']),
                              str(f['filename']))
                 for f in file_records]
        return files

    def _get_file_records(self, obsid=None, start=None, stop=None,
                          revision=None, content=None):
        self.set_read_env()
        tstart_pad = 10 * 86400
        if content is None:
            content = self.config['content_types']
        if type(content) == str:
            content = [content]
        content_str = ','.join(["'%s'" % x for x in content])
        if obsid is None:
            if start is None or stop is None:
                raise TypeError("Must supply either obsid or start and stop")
            tstart = DateTime(start).secs
            tstop = DateTime(stop).secs
            db_query = ("SELECT * from archfiles "
                        "WHERE tstart >= %f - %f "
                        "AND tstart < %f "
                        "AND tstop > %f "
                        "AND content in (%s) "
                        % (tstart, tstart_pad, tstop, tstart, content_str))
        else:
            db_query = ('SELECT * from archfiles '
                        'WHERE obsid = %d '
                        'AND content in (%s) '
                        % (obsid, content_str))
        if revision is None:
            db_query += 'AND isdefault = 1 '
        else:
            if revision == 'last':
                db_query += """AND revision in
                               (SELECT max(revision) from archfiles
                                WHERE obsid = %d)""" % obsid
            elif revision == 'all':
                pass
            else:
                db_query += 'AND revision = %d ' % revision
        db_query += "order by tstart"
        files = self._archfiles_db.fetchall(db_query)
        return files

    def get_dir(self, obsid):
        """
        Return the latest released directory for an obsid
        Return None if there are no 'default' / released products.
        """
        dirmap = self.get_obs_dirs(obsid)
        if 'default' in dirmap:
            return dirmap['default']
        else:
            return None

    def get_obs_dirs(self, obsid):
        """
        Return a dictionary of the directories available for an obsid.
        This is just done with a glob in the data directories.
        """
        data_root = self._config['data_root']
        strobs = "%05d" % obsid
        chunk_dir = strobs[0:2]
        topdir = os.path.join(data_root, chunk_dir)
        dirmap = dict(revisions=[])
        verdirs = glob(os.path.join(topdir, "%s_v*" % strobs))
        if not verdirs:
            return None
        for v in verdirs:
            nmatch = re.search("%s_v(\d{2})" % strobs, v)
            if nmatch:
                dirmap[int(nmatch.group(1))] = v
                dirmap['revisions'].append(int(nmatch.group(1)))
        lastdirs = glob(os.path.join(topdir, "%s_last" % strobs))
        defdirs = glob(os.path.join(topdir, "%s" % strobs))
        if defdirs:
            dirmap['default'] = defdirs[0]
        if lastdirs:
            dirmap['last'] = lastdirs[0]
        else:
            if defdirs:
                dirmap['last'] = defdirs[0]
        return dirmap

    @staticmethod
    def get_file_ver(tempdir, fileglob, ver_regex):
        """
        Determine the version/revision of a set of archived files from
        their file names.

        :param tempdir: directory containing files
        :param fileglob: glob to match files in question
        :param ver_regex: regular expression to pull out version
                          from the set of files

        :returns: version number
        :rtype: integer
        """
        files = glob(os.path.join(tempdir, fileglob))
        if not files:
            return None
        versions = {}
        for f in files:
            fmatch = re.search(ver_regex, f)
            if fmatch:
                versions[int(fmatch.group(1))] = 1
            if len(versions) > 1:
                raise ValueError("Different version files in %s" % tempdir)
        # update version to number
        version = versions.keys()[0]
        return version

    def get_ver_num(self, obsid, version='default'):
        """
        Determine the version number associated with the current released
        products or with the products referenced by "version=last".

        :param obsid: obsid
        :param version: version string ('default'|'last')

        :returns: version
        :rtype: integer
        """
        arc5 = self._arc5
        apstat = self._apstat
        logger = self.logger
        config = self.config
        # this is obi agnostic, all obis should be same
        # version anyway...u
        tempdir = tempfile.mkdtemp()
        # if multi-obi, we'll need to be specific
        arc5.sendline("reset")
        arc5.sendline("cd %s" % tempdir)
        arc5.sendline("obsid=%d" % obsid)

        obis = apstat.fetchall(
            "select distinct obi from obidet_0_5 where obsid = %d" % obsid)
        if len(obis) > 1:
            minobi = np.min(obis['obi'])
            logger.info("limiting arc5gl to obi %d" % minobi)
            arc5.sendline("obi=%d" % minobi)
        if version != 'default':
            arc5.sendline("version=%s" % version)
        # just get a small file
        arc5.sendline("get %s" % config['small'])
        version = self.get_file_ver(tempdir,
                                    config['small_glob'],
                                    config['small_ver_regex'],
                                    )
        for afile in glob(os.path.join(tempdir, "*")):
            os.remove(afile)
        os.removedirs(tempdir)
        return version

    def get_fits_info(self, i, f, archfiles):
        """
        Read FITS file ``f`` with index ``i`` (position within list of
        filenames ``archfiles``) and get dictionary of values to store
        in file lookup database.  This values include all header key/value
        pairs with keys in ``config[cols]`` plus the header checksum,
        the filename, the year, and day-of-year.

        :param i: index of file f within list of files archfiles
        :param f: filename
        :param archfiles: list of filenames for this batch

        :returns: info for a file
        :rtype: dictionary
        """
        logger = self.logger
        config = self.config
        filename = os.path.basename(f)

        logger.debug('Reading (%d / %d) %s' % (i, len(archfiles), filename))
        hdus = pyfits.open(f)
        hdu = hdus[1]

        # Accumulate relevant info about archfile that will be ingested
        archfiles_row = dict((x, hdu.header.get(x.upper()))
                             for x in config['cols'])
        archfiles_row['checksum'] = hdu.header.get('checksum') or hdu._checksum
        archfiles_row['filename'] = filename
        archfiles_row['filetime'] = int(
            re.search(r'(\d+)', archfiles_row['filename']).group(1))
        filedate = DateTime(archfiles_row['filetime']).date
        year, doy = (
            int(x) for x
            in re.search(r'(\d\d\d\d):(\d\d\d)', filedate).groups())
        archfiles_row['year'] = year
        archfiles_row['doy'] = doy
        hdus.close()
        return archfiles_row

    def get_obspar_info(self, i, f, archfiles):
        """
        Wrap get_all_obspar_info() and just include columns in config['cols']
        """
        obspar = self.get_all_obspar_info(i, f, archfiles)
        return {col: obspar[col]
                for col in self.config['cols'] if col in obspar}

    def get_arch_info(self, i, f, archfiles):
        """
        Get information for a file for the file lookup table/database.
        For obspars, call the get_obspar_info() method.
        For FITS files, call get_fits_info() method.
        """
        config = self.config
        if config['full'] == 'obspar':
            return self.get_obspar_info(i, f, archfiles)
        else:
            return self.get_fits_info(i, f, archfiles)

    def get_arch(self, obsid, version='last'):
        """
        Retrieve telemetry for an observation from the CXC archive and
        store in the Ska file archive.

        :param obsid: obsid
        :param version: 'default', 'last', or revision/version number
        :returns: obsid directory in Ska file archive
        :rtype: directory string
        """

        arc5 = self._arc5
        apstat = self._apstat
        logger = self.logger
        config = self.config
        temp_root = config['temp_root']


        # get a numeric version
        # do this even if passed a numeric version, as this will
        # check to see if the version is available
        n_version = self.get_ver_num(obsid, version=version)
        if n_version is None:
            raise ProductVersionError("No %s data for ver %s"
                                      % (config['label'], version))
        # use the numeric version instead of 'last' or 'default'
        version = n_version

        # set up directory for data
        strobs = "%05d_v%02d" % (obsid, version)
        chunk_dir = strobs[0:2]
        obs_dir = os.path.join(config['data_root'], chunk_dir, strobs)
        if not os.path.exists(obs_dir):
            logger.info("making directory %s" % obs_dir)
            os.makedirs(obs_dir)
        else:
            logger.info("obsid dir %s already exists" % obs_dir)
            if not config['filecheck']:
                return obs_dir

        # get data
        if not os.path.exists(temp_root):
            os.makedirs(temp_root)
        tempdirobj = Ska.File.TempDir(dir=os.path.abspath(temp_root))
        tempdir = tempdirobj.name
        arc5.sendline("reset")
        arc5.sendline("cd %s" % tempdir)
        logger.info("retrieving data for %d in %s" % (obsid, tempdir))
        arc5.sendline("obsid=%d" % obsid)
        # if multi-obi, limit to just the first obi
        obis = apstat.fetchall(
            "select distinct obi from obidet_0_5 where obsid = %d" % obsid)
        if len(obis) > 1:
            minobi = np.min(obis['obi'])
            logger.info("limiting arc5gl to obi %d" % minobi)
            arc5.sendline("obi=%d" % minobi)
        arc5.sendline("version=%s" % version)
        arc5.sendline("get %s" % config['full'])
        # get the log too
        arc5.sendline("dataset=pipelog")
        arc5.sendline("go")
        archfiles = glob(os.path.join(tempdir, "*"))
        if not archfiles:
            raise ValueError("Retrieved no files")
        db = self._archfiles_db
        existing = db.fetchall(
            "select * from archfiles where obsid = %d and revision = '%s'"
            % (obsid, version))
        for i, f in enumerate(archfiles):
            # just copy the log files without inserting in db
            if re.match('.+log.gz', f):
                os.chmod(f, 0775)
                shutil.copy(f, obs_dir)
                os.remove(f)
                continue
            arch_info = self.get_arch_info(i, f, archfiles)
            arch_info['obsid'] = obsid
            if (len(existing)
                    and arch_info['filename'] in existing['filename']):
                logger.debug("skipping %s" % f)
                os.remove(f)
                continue

            db.insert(arch_info, 'archfiles')
            os.chmod(f, 0775)
            shutil.copy(f, obs_dir)
            os.remove(f)
        db.commit()
        return obs_dir

    # this needs help
    def update_link(self, obsid):
        """
        Create links in the obsid data directories to make it easy to find
        the current 'default'/released data, all versions that have been
        archived, and the 'last'/unreleased/provisional data if available.

        This is designed so that if obsid 5 has released data in version 1
        and provisional data in version 2, that the directories and links
        will look like:

        directory 00005_v01
        directory 00005_v02
        link      00005 -> 00005_v01
        link      00005_last -> 00005_v02
        """
        logger = self.logger
        config = self.config
        db = self._archfiles_db
        archive_dir = config['data_root']

        # for the obsid get the default and last version numbers
        default_ver = self.get_ver_num(obsid, 'default')
        last_ver = self.get_ver_num(obsid, 'last')

        # directories for this obsid are in 'chunk_dir'
        chunk_dir = ("%05d" % obsid)[0:2]
        chunk_dir_path = os.path.join(archive_dir, chunk_dir)
        obs_ln = os.path.join(archive_dir, chunk_dir, "%05d" % obsid)
        obs_ln_last = os.path.join(archive_dir, chunk_dir,
                                   "%05d_last" % obsid)

        if default_ver is not None:
            def_ver_dir = os.path.join(archive_dir, chunk_dir,
                                       '%05d_v%02d' % (obsid, default_ver))
            # link the default version
            logger.info("linking %s -> %s" % (
                os.path.relpath(def_ver_dir, chunk_dir_path),
                obs_ln))
            # don't use os.path.exists here because if we have a broken
            # link we want to remove it too
            if os.path.islink(obs_ln):
                os.unlink(obs_ln)
            os.symlink(
                os.path.relpath(def_ver_dir, chunk_dir_path),
                obs_ln)
            logger.info("updating archfiles default rev to %d for %d"
                        % (default_ver, obsid))
            db.execute("""UPDATE archfiles SET isdefault = 1
                          WHERE obsid = %d and revision = %d"""
                       % (obsid, default_ver))
            db.execute("""UPDATE archfiles SET isdefault = NULL
                          WHERE obsid = %d and revision != %d"""
                       % (obsid, default_ver))
            db.commit()

        # nothing to do if last_ver undefined
        if last_ver is None:
            return

        last_ver_dir = os.path.join(archive_dir, chunk_dir,
                                    '%05d_v%02d' % (obsid, last_ver))
        # if last and default are different, and we have last data,
        # make a link to it
        if last_ver != default_ver:
            if os.path.exists(last_ver_dir):
                obs_ln_last = os.path.join(archive_dir, chunk_dir,
                                           "%05d_last" % obsid)
                logger.info("linking %s -> %s" % (
                    os.path.relpath(last_ver_dir, chunk_dir_path),
                    obs_ln_last))
                if os.path.islink(obs_ln_last):
                    os.unlink(obs_ln_last)
                os.symlink(
                    os.path.relpath(last_ver_dir, chunk_dir_path),
                    obs_ln_last)
        else:
            # if default and last are now the same, delete last link
            if os.path.exists(obs_ln_last):
                logger.info("removing outdated link %s"
                            % obs_ln_last)
                os.remove(obs_ln_last)

    def get_todo_from_links(self, archive_dir):
        """
        Return a list of all of the *_last directories in the file archive
        (and specify revision=default to attempt to get new released products
        for them).
        """
        logger = self.logger
        logger.info("Checking for updates to obsids with provisional data")
        chunk_dirs = glob(os.path.join(archive_dir, "??"))
        todo_obs = []
        for cdir in chunk_dirs:
            last_links = glob(os.path.join(cdir, "?????_last"))
            for link in last_links:
                target = os.readlink(link)
                lmatch = re.search('(\d{5})_v(\d+)$', target)
                if lmatch:
                    obs = dict(obsid=int(lmatch.group(1)),
                               revision=int(lmatch.group(2)))
                    todo_obs.append(obs)
        return todo_obs

    def update(self):
        """
        Run the update process using the config already passed to the object.
        """

        self.set_env()
        config = self._config
        logger = self.logger
        archive_dir = config['data_root']
        apstat = self._apstat
        updated_obsids = []

        if (config['data_root'].startswith('/data/aca/archive')
                and not mica_version.release):
            raise ValueError(
                "non-release code attempting to write to official archive")

        last_id_file = os.path.join(archive_dir, 'last_id.txt')
        # if an obsid is requested, just do that
        if 'obsid' in config and config['obsid']:
            self.get_arch(config['obsid'], config['version'])
            try:
                self.update_link(config['obsid'])
            except ProductVersionError as ve:
                logger.debug(ve)
                logger.info("Could not determine link versions for %d"
                            % config['obsid'])
            return [config['obsid']]

        # if no obsid specified, try to retrieve all data
        # since tool last run
        # use saved id from the apstat table as a reference
        last_id = 0
        if os.path.exists(last_id_file):
            last_id_fh = open(last_id_file)
            last_id = int(last_id_fh.read().rstrip())
            logger.info("using %d for last_id" % last_id)
            last_id_fh.close()
        else:
            if not config.get('rebuild'):
                raise ValueError("last_id.txt not found.\n"
                                 + "To rebuild archive from obsid 1, "
                                 + "set rebuild=True in configuration")
        query_vars = {'apstat_table': config['apstat_table'],
                      'apstat_id': config['apstat_id'],
                      'last_id': last_id}
        apstat_query = ("""select * from %(apstat_table)s
                                  where %(apstat_id)s > %(last_id)d
                                  order by %(apstat_id)s"""
                        % query_vars)
        logger.debug(apstat_query)
        todo = apstat.fetchall(apstat_query)
        for obs in todo:
            logger.info("running get_arch for obsid %d run on %s"
                        % (obs['obsid'], obs['ap_date']))
            get_rev = obs['revision']
            # firstrun was the mode to initially populate the
            # file archive
            if config['firstrun']:
                # if I've already got a later version, skip this
                # obsid
                have = self.get_obs_dirs(obs['obsid'])
                if have:
                    max_rev = max(have['revisions'])
                    logger.info("%d vs %d" % (max_rev, obs['revision']))
                    if max_rev >= obs['revision']:
                        logger.info(
                            "skipping %d, firstrun and already have newer rev"
                            % obs['obsid'])
                        continue
                # otherwise override the revision from the apstat
                # table and just get the current 'default'
                else:
                    get_rev = 'default'

            try:
                self.get_arch(obs['obsid'], get_rev)
                self.update_link(obs['obsid'])
            except ProductVersionError as ve:
                logger.info("skipping %d, default ver not available"
                            % obs['obsid'])
                logger.debug(ve)
                continue
            # update the file that stores the last id from the apstat table
            last_id_fh = open(last_id_file, 'w')
            last_id_fh.write("%d" % obs[config['apstat_id']])
            last_id_fh.close()
            updated_obsids.append(obs['obsid'])
        if not len(todo):
            logger.info("No new data")

        # look for all previous "provisional" aspect solutions
        # and broken links and and update if possible
        prov_data = self.get_todo_from_links(archive_dir)
        for obs in prov_data:
            # check again for multi-obis and limit to first one
            obis = apstat.fetchall(
                "select distinct obi from obidet_0_5 where obsid = %d"
                % obs['obsid'])
            minobi = np.min(obis['obi'])
            logger.info("checking database status for obsid %d obi %d ver %s, "
                        % (obs['obsid'], minobi, obs['revision']), )
            query_vars = {'apstat_table': config['apstat_table'],
                          'apstat_id': config['apstat_id'],
                          'obsid': obs['obsid'],
                          'obi': minobi,
                          'revision': obs['revision']}
            apstat_query = ("""select * from %(apstat_table)s
                               where obsid = %(obsid)d and obi = %(obi)d
                               and revision = %(revision)d"""
                            % query_vars)
            logger.debug(apstat_query)
            current_status = apstat.fetchall(apstat_query)
            if len(current_status) == 0:
                raise ValueError(
                    "obsid %(obsid)d revision %(revision)d not in %(apstat_table)s"
                    % query_vars)
            if len(current_status) > 1:
                raise ValueError(
                    "obsid %(obsid)d revision %(revision)d multiple entries in %(apstat_table)s"
                    % query_vars)
            # a query to get the quality from the max science_2 data that
            # used this aspect_solution.  Ugh.
            if config['apstat_table'] == 'aspect_1':
                science_qual = apstat.fetchall(
                    """select quality from science_2 where science_2_id in (
                         select science_2_id from science_2_obi where science_1_id in (
                         select science_1_id from science_1 where aspect_1_id in (
                         select aspect_1_id from aspect_1
                         where obsid = {obsid} and obi = {obi}
                         and revision = {revision})))""".format(query_vars))
                # if any of the science_2 data associated with this obsid is
                # now not pending, set the quality to an arbitrary value 'X'
                # and try to update the obsid
                if np.any(science_qual['quality'] != 'P'):
                    current_status['quality'] = 'X'
            if current_status['quality'] != 'P':
                try:
                    self.get_arch(obs['obsid'], 'default')
                except ProductVersionError as ve:
                    logger.info("skipping %d, default ver not available"
                                % obs['obsid'])
                    logger.debug(ve)
                    continue
                self.update_link(obs['obsid'])
            updated_obsids.append(obs['obsid'])
        return updated_obsids
