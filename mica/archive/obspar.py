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
# these are set as globals in main()
archive_dir = None
arc5 = None
#aca_db = None
apstat = None
SKA = os.environ['SKA']

logger = logging.getLogger('fetch')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


#from configobj import ConfigObj
#config = ConfigObj("obspar.conf")
config = dict(data_root='./data/aca/archive/obspar',
              temp_root='./data/aca/archive/tempobs',
              apstat_table='obidet_0_5',
              apstat_id='obidet_0_5_id',
              label='obspar',
              small='obspar',
              small_glob='axaff*par*',
              small_ver_regex='axaff\d{5}_\d{3}N(\d{3}).*',
              full='obspar')


archfiles_hdr_cols = [
'filename',
'obsid',
'title',
'observer',
'ao',
'object',
'ss_object',
'obs_id',
'obi_num',
'seq_num',
'instrume',
'grating',
'detector',
'detnam',
'si_mode',
'optical_monitor',
'raster_scan',
'dither',
'dither_y_amp',
'dither_y_freq',
'dither_y_phase',
'dither_z_amp',
'dither_z_freq',
'dither_z_phase',
'ra_pnt',
'dec_pnt',
'roll_pnt',
'ra_targ',
'dec_targ',
'y_det_offset',
'z_det_offset',
'radial_offset',
'defocus',
'sim_z_offset',
'pre_id',
'uninterrupted',
'seg_max_num',
'ra_nom',
'dec_nom',
'roll_nom',
'date_obs',
'date_end',
'tstart',
'tstop',
'sched_start',
'sched_stop',
'sched_exp_time',
'obs_mode',
'maneuver',
'maneuver_v1',
'maneuver_v2',
'maneuver_v3',
'maneuver_angle',
'maneuver_ref',
'mjdref',
'timezero',
'timeunit',
'timesys',
'timversn',
'datamode',
'readmode',
'ccdi0_on',
'ccdi1_on',
'ccdi2_on',
'ccdi3_on',
'ccds0_on',
'ccds1_on',
'ccds2_on',
'ccds3_on',
'ccds4_on',
'ccds5_on',
'dropped_chip_count',
'onchip_sum',
'sumrow',
'sumcol',
'subarray',
'startrow',
'rowcnt',
'subarray_frame_time',
'duty_cycle',
'dtycycle',
'exptimea',
'exptimeb',
'most_efficient',
'eventfilter',
'phamin',
'pharange',
'bias_request',
'mode',
'mission',
'telescop',
'sim_x',
'sim_y',
'sim_z',
'foc_len',
'py_shutter_position',
'range_switch_level',
'antico_enable',
'upper_level_disc',
'timing_mode',
'trigger_level',
'u_blank_hi',
'u_blank_low',
'v_blank_low',
'v_blank_hi',
'uld_enable',
'zero_block',
'blank_enable',
'width_enable',
'spect_mode',
'my_shutter_position',
'width_threshold',
'mt_a',
'mt_aop',
'mt_e',
'mt_epoch',
'mt_i',
'mt_ma',
'mt_raan',
'timeref',
'tassign',
'origin',
'ascdsver',
'obi_data',
'revision',
'obspar_ver',
'obspar_type',
'obspar_stat']

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


def get_obspar_info(i, f, archfiles):
    filename = os.path.basename(f)
    logger.debug('Reading (%d / %d) %s' % (i, len(archfiles), filename))
    obspar = get_obspar(f)
    obspar['obsid'] = obspar['obs_id']
    obspar['filename'] = filename
    return obspar
    

def get_dir(obsid, data_root=config['data_root']):
    """Return the latest released directory for an obsid."""
    dirmap = get_obs_dirs(obsid, data_root)
    return dirmap['default']


def get_obs_dirs(obsid, data_root=config['data_root']):
    """Return a dictionary of all of the directories available for an obsid."""

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


def get_ver_num(obsid, version='default'):
    # this is obi agnostic, all obis should be same
    # version anyway...
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
    # just fetch the aspect solution to start
    arc5.sendline("get %s" % config['small'])
    version = get_file_ver(tempdir,
                           config['small_glob'],
                           config['small_ver_regex'],
                           )
    if not version:
        raise ValueError("Version not defined")
    for afile in glob(os.path.join(tempdir, "*")):
        os.remove(afile)
    os.removedirs(tempdir)
    return version


def get_arch(obsid, version='last', temp_root=None):
    n_version = get_ver_num(obsid, version=version)
    if n_version is None:
        raise ValueError("No %s data for ver %s" 
                         % (config['label'], version))
    # use numeric version instead of 'last' or 'default'
    version = n_version

    # set up directory for data
    strobs = "%05d_v%02d" % (obsid, version)
    chunk_dir = strobs[0:2]
    obs_dir = os.path.join(archive_dir, chunk_dir, strobs)
    if not os.path.exists(obs_dir):
        logger.info("making directory %s" % obs_dir)
        os.makedirs(obs_dir)
    else:
        logger.info("obsid dir %s already exists" % obs_dir)

    # get data
    tempdirobj = Ska.File.TempDir(dir=temp_root)
    tempdir = os.path.abspath(tempdirobj.name)
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
        db_file = os.path.join(os.path.abspath(archive_dir), 'archfiles.db3')
        if not os.path.exists(db_file):
            db_sql = os.path.join(os.environ['SKA_DATA'],
                                  'mica', 'obspar_def.sql')
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
            obspar = get_obspar_info(i, f, archfiles)
            arch_info = dict()
            [arch_info.update({col: obspar[col]})
             for col in archfiles_hdr_cols if col in obspar]
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
def update_link(obsid):
    # for the obsid get the default and last version numbers
    default_ver = get_ver_num(obsid, 'default')
    last_ver = get_ver_num(obsid, 'last')

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


def get_todo_from_links(archive_dir):
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





def set_env(opt):

    global archive_dir
    archive_dir = opt.data_root
    global arc5
    arc5 = Ska.arc5gl.Arc5gl()
    global apstat
    apstat = Ska.DBI.DBI(dbi='sybase', server='sqlsao', database='axafapstat')


def update_archive(opt):

    set_env(opt)

    last_id_file = os.path.join(archive_dir, 'last_id.txt')
    # if an obsid is requested, just do that
    if opt.obsid:
        get_arch(opt.obsid, opt.version, opt.temp_root)
        update_link(opt.obsid)
        return

    # look for all previous "provisional" aspect solutions
    # and broken links and and update if possible
    prov_data = get_todo_from_links(archive_dir)
    for obs in prov_data:
        logger.info("running get_arch for obsid %d ver %s, "
                    % (obs['obsid'], obs['revision']))
        try:
            get_arch(obs['obsid'], obs['revision'], opt.temp_root)
        except ValueError as ve:
            logger.info("skipping %d, default ver not available"
                        % obs['obsid'])
            logger.debug(ve)
            continue
        update_link(obs['obsid'])
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
        if opt.firstrun:
            # if I've already got a later version, skip
            have = get_obs_dirs(obs['obsid'], archive_dir)
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
                get_arch(obs['obsid'], 'default', opt.temp_root)
            except ValueError as ve:
                logger.info("skipping %d, default ver not available"
                            % obs['obsid'])
                logger.debug(ve)
                continue
        else:
            get_arch(obs['obsid'], obs['revision'], opt.temp_root)
        update_link(obs['obsid'])
        last_id_fh = open(last_id_file, 'w')
        last_id_fh.write("%d" % obs[config['apstat_id']])
        last_id_fh.close()
    if not len(todo):
        logger.info("No new data")

def get_options():
#    from optparse import OptionParser
#    parser = OptionParser()
    import argparse
    parser = argparse.ArgumentParser(
        description="Fetch aspect level 1 products and make a file archive")
    defaults = dict(config)
    parser.set_defaults(**defaults)
    parser.add_argument("--obsid",
                        type=int,
                        help="specific obsid to process")
    parser.add_argument("--version",
                        default='last',
                        help="specific processing version to retrieve")
    parser.add_argument("--firstrun",
                        action='store_true',
                        help="for archive init., ignore rev in aspect_1 table")
    parser.add_argument("--data-root",
                        help="parent directory for all data")
    parser.add_argument("--temp-root",
                        help="parent temp directory")
    parser.add_argument("--proctype",
                        default="asp1")
    opt = parser.parse_args()
    return opt 


def main():
    opt = get_options()
    update_archive(opt)

if __name__ == '__main__':
    main()
