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


# borrowed from telem_archive
import csv
import gzip

def parse_obspar(file):
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
        self._arc5 = Ska.arc5gl.Arc5gl()
        self._apstat = Ska.DBI.DBI(dbi='sybase', server='sqlsao',
                                   database='axafapstat')


    def get_obspar_info(self, i, f, archfiles):
        logger = self.logger
        filename = os.path.basename(f)
        logger.debug('Reading (%d / %d) %s' % (i, len(archfiles), filename))
        obspar = get_obspar(f)
        obspar['obsid'] = obspar['obs_id']
        obspar['filename'] = filename
        return obspar


    def get_dir(self, obsid):
        """Return the latest released directory for an obsid."""
        data_root = self._config['data_root']
        dirmap = get_obs_dirs(obsid, data_root)
        return dirmap['default']

    def get_obs_dirs(self, obsid):
        """Return a dictionary of all of the directories available for an obsid."""
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
        # just a small file to start
        arc5.sendline("get %s" % config['small'])
        version = self.get_file_ver(tempdir,
                               config['small_glob'],
                               config['small_ver_regex'],
                               )
        if not version:
            raise ProductVersionError("Version not defined")
        for afile in glob(os.path.join(tempdir, "*")):
            os.remove(afile)
        os.removedirs(tempdir)
        return version

    def get_fits_info(self, i, f, archfiles):
        """Read filename ``f`` with index ``i`` (position within list of
        filenames).  The file has type ``filetype`` and will be added to
        MSID file at row index ``row``.  ``colnames`` is the list of
        column names for the content type (not used here).
        """
        logger = self.logger
        config = self.config
        filename = os.path.basename(f)

        # Read FITS archive file and accumulate data into dats list and
        # header into headers dict
        logger.debug('Reading (%d / %d) %s' % (i, len(archfiles), filename))
        hdus = pyfits.open(f)
        hdu = hdus[1]

        # Accumlate relevant info about archfile that will be ingested into
        # MSID h5 files.  Commit info before h5 ingest so if there is a failure
        # the needed info will be available to do the repair.
        archfiles_row = dict((x, hdu.header.get(x.upper()))
                             for x in config['cols'])
        archfiles_row['checksum'] = hdu._checksum
        archfiles_row['filename'] = filename
        archfiles_row['filetime'] = int(
            re.search(r'(\d+)', archfiles_row['filename']).group(1))
        filedate = DateTime(archfiles_row['filetime']).date
        year, doy = (int(x)
                     for x in re.search(r'(\d\d\d\d):(\d\d\d)', filedate).groups())
        archfiles_row['year'] = year
        archfiles_row['doy'] = doy
        hdus.close()
        return archfiles_row

    
    def get_obspar_info(self, i, f, archfiles):
        obspar = self.get_obspar_info(i, f, archfiles)
        arch_info = dict()
        [arch_info.update({col: obspar[col]})
         for col in config['cols'] if col in obspar]
        return arch_info


    def get_arch_info(self, i, f, archfiles):
        config = self.config
        if config['full'] == 'obspar':
            return self.get_obspar_info(i, f, archfiles)
        else:
            return self.get_fits_info(i, f, archfiles)

    def get_arch(self, obsid, version='last'):
        arc5 = self._arc5
        apstat = self._apstat
        logger = self.logger
        config = self.config
        temp_root = config['temp_root']

        n_version = self.get_ver_num(obsid, version=version)
        if n_version is None:
            raise ProductVersionError("No %s data for ver %s" 
                             % (config['label'], version))
        # use numeric version instead of 'last' or 'default'
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

        # get data
        if not os.path.exists(temp_root):
            os.makedirs(temp_root)
        tempdirobj = Ska.File.TempDir(dir=os.path.abspath(temp_root))
        tempdir = tempdirobj.name
        arc5.sendline("reset")
        arc5.sendline("cd %s" % tempdir)
        logger.info("retrieving data for %d in %s" % (obsid, tempdir))
        arc5.sendline("obsid=%d" % obsid)
        # if multi-obi, we'll need to be specific
        obis = apstat.fetchall(
            "select distinct obi from obidet_0_5 where obsid = %d" % obsid)
        if len(obis) > 1:
            minobi = np.min(obis['obi'])
            logger.info("limiting arc5gl to obi %d" % minobi)
            arc5.sendline("obi=%d" % minobi)
        arc5.sendline("version=%s" % version)
        arc5.sendline("get %s" % config['full'])
        archfiles = glob(os.path.join(tempdir, "*"))
        if not archfiles:
            raise ValueError("Retrieved no files")
        if archfiles:
            db_file = os.path.join(os.path.abspath(config['data_root']),
                                   'archfiles.db3')
            if not os.path.exists(db_file):
                logger.info("creating archfiles db from %s"
                            % config['sql_def'])
                db_sql = os.path.join(os.environ['SKA_DATA'],
                                      'mica', config['sql_def'])
                db_init_cmds = file(db_sql).read()
                db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
                                 autocommit=False)
                db.execute(db_init_cmds, commit=True)
            db = Ska.DBI.DBI(dbi='sqlite', server=db_file,
                             autocommit=False)
            existing = db.fetchall(
                "select * from archfiles where obsid = %d and revision = '%s'"
                % (obsid, version))
            for i, f in enumerate(archfiles):
                arch_info = self.get_arch_info(i, f, archfiles)
                arch_info['obsid'] = obsid
                if (len(existing)
                   and arch_info['filename'] in existing['filename']):
                    print "skipping %s" % f
                    os.remove(f)
                    continue

                db.insert(arch_info, 'archfiles')
                shutil.copy(f, obs_dir)
                os.remove(f)
        #os.removedirs(tempdir)
            db.commit()
        return obs_dir, chunk_dir


    # this needs help
    def update_link(self, obsid):
        arc5 = self._arc5
        apstat = self._apstat
        logger = self.logger
        config = self.config
        archive_dir = config['data_root']

        # for the obsid get the default and last version numbers
        default_ver = self.get_ver_num(obsid, 'default')
        last_ver = self.get_ver_num(obsid, 'last')

        # set up directory for data
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
            # if default and last are the same and we have a last link,
            # delete it
            if os.path.exists(obs_ln_last):
                logger.info("removing outdated link %s"
                            % obs_ln_last)
                os.remove(obs_ln_last)

    def get_todo_from_links(self, archive_dir):
        logger = self.logger
        logger.info("Checking for updates to obsids with provisional data")
        chunk_dirs = glob(os.path.join(archive_dir, "??"))
        todo_obs = []
        for cdir in chunk_dirs:
            last_links = glob(os.path.join(cdir, "?????_last"))
            for link in last_links:
                lmatch = re.search('(\d{5})(_last)?$', link)
                if lmatch:
                    obs = dict(obsid=int(lmatch.group(1)),
                               revision='default')
                    todo_obs.append(obs)
        return todo_obs

    def update(self):
        self.set_env()
        config = self._config
        logger = self.logger
        archive_dir = config['data_root']
        apstat = self._apstat

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
            return

        # look for all previous "provisional" aspect solutions
        # and broken links and and update if possible
        prov_data = self.get_todo_from_links(archive_dir)
        for obs in prov_data:
            logger.info("running get_arch for obsid %d ver %s, "
                        % (obs['obsid'], obs['revision']))
            try:
                self.get_arch(obs['obsid'], obs['revision'])
            except ProductVersionError as ve:
                logger.info("skipping %d, default ver not available"
                            % obs['obsid'])
                logger.debug(ve)
                continue
            self.update_link(obs['obsid'])
         # if no obsid specified, try to retrieve all asp_1 runs
        # since tool last run
        # use an aspect_1_id in a file
        last_id = 0
        if os.path.exists(last_id_file):
            last_id_fh = open(last_id_file)
            last_id = int(last_id_fh.read().rstrip())
            logger.info("using %d for last_id" % last_id)
            last_id_fh.close()
        query_vars = dict(config)
        query_vars.update({'last_id': last_id})
        apstat_query = ("""select * from %(apstat_table)s
                                  where %(apstat_id)s > %(last_id)d
                                  order by %(apstat_id)s"""
                        % query_vars)
        logger.debug(apstat_query)
        todo = apstat.fetchall(apstat_query)
        for obs in todo:
            logger.info("running get_arch for obsid %d run on %s"
                        % (obs['obsid'], obs['ap_date']))
            if config['firstrun']:
                # if I've already got a later version, skip
                have = self.get_obs_dirs(obs['obsid'])
                if have:
                    max_rev = max(have['revisions'])
                    print "%d vs %d" % (max_rev, obs['revision'])
                    if max_rev >= obs['revision']:
                        logger.info("skipping %d, firstrun and already have newer rev"
                                    % obs['obsid'])
                        continue

                # ignore aspect_1 table versions and just try for default
                # also ignore errors in this case
                try:
                    self.get_arch(obs['obsid'], 'default')
                except ProductVersionError as ve:
                    logger.info("skipping %d, default ver not available"
                                % obs['obsid'])
                    logger.debug(ve)
                    continue
            else:
                self.get_arch(obs['obsid'], obs['revision'])
            self.update_link(obs['obsid'])
            last_id_fh = open(last_id_file, 'w')
            last_id_fh.write("%d" % obs[config['apstat_id']])
            last_id_fh.close()
        if not len(todo):
            logger.info("No new data")
