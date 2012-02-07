#!/usr/bin/env python

import os
import re
from glob import glob
import pyfits
import json
import string
import random
import csv
import gzip

import Ska.arc5gl
from Ska.Shell import bash
import Ska.Table

job_dir = os.path.abspath("./jobs")


def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--obsid",
                      type='int',
                      help="obsid to process")
    parser.add_option("--archive",
                      type='int',
                      help="obsid to fetch from archive and process")
    #parser.add_option("--version",
    #                  help="version of products )
    #parser.add_option("--release")
    #parser.add_option("--clean",
    #                  action='store_true')
    parser.add_option("--revision",
                      default="1",
                      type="int",
                      help="integer revision label for output")
    parser.add_option("--skip-slot",
                      type='int',
                      action='append',
                      help="slots to skip in processing")
    parser.add_option('--skip-slot-method',
                      default='telem',
                      help="method to cut slots")
    parser.add_option("--range",
                      help="processing range specifier")
    parser.add_option("--dir",
                      default=".",
                      help="directory for telem fetching")
    parser.add_option("--jobdir",
                      default=job_dir,
                      help="directory for job queue json files")
    opt, args = parser.parse_args()
    return opt, args


pcad0_files = ['pcad0/pcad*_%d_*fits*' % x for x in [3, 5, 7, 8, 14, 15]]

pipe_config = dict(pipe_ped='asp_l1_std',
                   infiles=pcad0_files
                   + ['asp05/pcad*aipr0a.fits*',
                      'asp05/pcad*cai0a.*',
                      'aca0/*',
                      'sim05/sim*coor0a*',
                      'obspar/axaf*obs0a.par*',
                      'acis2eng/acis*fits*',
                      'obc4eng/obc*fits*'],
                   outfiles=[],
                   archfiles=[('aca0', 'aca0'),
                              ('pcad0', 'pcad0'),
                              ('sim05', 'sim05'),
                              ('asp05', 'asp05'),
                              ('obspar', 'obspar'),
                              ('acis2eng', 'acis_eng_0{acis2eng}'),
                              ('obc4eng', 'obc_eng_0{obc4eng}')])

# the cai files have non-compliant IRAF par files.
# override the data types
cai_override = {'obs_id': 'i',
                'obi_num': 'i',
                'ascdsver': 's',
                'num_cai': 'i'}


def parse_obspar(file, override=None):
# borrowed from telem_archive ... override is new
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
        if override and (row['name'] in override):
            row['value'] = convert[override[row['name']]](row['value'])
        else:
            row['value'] = convert[row['type']](row['value'])
        row['name'] = row['name'].replace('-', '_')
        yield row

    return


def get_par(parfile, override=None):
    """
    Read an IRAF-style par file and return a dictionary.
    """
    par = {}
    for row in parse_obspar(parfile, override):
        par.update({row['name']: row['value']})
    return par


def get_obspar(obsparfile):
    """Get the obspar for obsid starting at tstart.  Return as a dict."""
    obspar = {'num_ccd_on': 0}
    for row in parse_obspar(obsparfile):
        obspar.update({row['name']: row['value']})
        if re.match(r'^ccd[is]\d_on$', row['name']) and row['value'] == 'Y':
            obspar['num_ccd_on'] += 1
    return obspar


def dir_setup(dir, istart):
    """
    Makes

      - a directory to contain processing for each aspect interval (AI)
      - at least one input directory in each AI directory which will be
        populated with telemetry or links to telemetry for processing.
      - at least one output directory in each AI directory

    The input and output directories have trailing integer revision numbers
    that are incremented in successive processing attempts.
    """

    workdir = os.path.join(dir,
                           'ASP_L1_STD_%09d' % istart)
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    rev = 1
    indir = os.path.join(workdir, "in%d" % rev)
    indirs = glob(os.path.join(workdir, "in*"))
    if len(indirs):
        rev = len(indirs) + 1
        indir = os.path.join(workdir, "in%d" % rev)
        if os.path.exists(indir):
            raise ValueError("Bad in directory sequence (%s exists)"
                             % indir)
    os.makedirs(indir)
    outdir = os.path.join(workdir, "out%d" % rev)
    if os.path.exists(outdir):
        raise ValueError("Bad in directory sequence (%s exists)"
                         % outdir)
    os.makedirs(outdir)
    return workdir, indir, outdir


def link_files(dir, indir, outdir, istart, istop, obiroot, skip_slot=None):
    """
    Creates symbolic links from the specified indir to the available telemetry.
    Fits files are only linked in if their header time keywords are relevant.
    ACA0 image files may be skipped if the slot is in skip_slot list.
    obspars must be for the correct obi.
    """
    dirmap = dict(infiles=indir,
                  outfiles=outdir)
    for filetype in ['infiles', 'outfiles']:
        ldir = dirmap[filetype]
        for fileglob in pipe_config[filetype]:
            match = glob(os.path.join(opt.dir, fileglob))
            for mfile in match:
                fitsmatch = re.match('.*fits', mfile)
                if fitsmatch:
                    header = pyfits.getheader(mfile)
                    if ((istart >= header['tstop'])
                        or (istop <= header['tstart'])):
                        #print "skipping file out of timerange %s" % mfile
                        continue
                    aca0 = re.search('aca.*_(\d)_img0', mfile)
                    if skip_slot and aca0:
                        aca_file_slot = int(aca0.group(1))
                        if aca_file_slot in skip_slot:
                            #print "skipping slot file on %s" % mfile
                            continue
                obsparmatch = re.match('.*obs0a\.par', mfile)
                if obsparmatch:
                    obimatch = re.match('.*axaf%s_obs0a\.par' % obiroot, mfile)
                    if not obimatch:
                        #print "skipping obspar for different obi"
                        continue
                #print("ln -s %s %s" % (os.path.relpath(mfile,ldir), ldir))
                bash("ln -s %s %s" % (os.path.relpath(mfile, ldir), ldir))


def make_list_files(dir, indir, outdir, root):
    """
    Create .lis files for the pipeline.
    """
    # remove any present already
    inlists = glob(os.path.join(indir, "*.lis"))
    [os.unlink(x) for x in inlists]
    # remove any present already
    outlists = glob(os.path.join(outdir, "*.lis"))
    [os.unlink(x) for x in outlists]

    for listend, listglob in (('sim.lis', 'sim*coor0a*fits'),
                              ('pcad.lis', 'pcad*eng0*fits'),
                              ('acis.list', 'acis*eng0*fits'),
                              ('obc.list', 'obc*eng0*fits')):
        lfile = open(os.path.join(indir, "%(root)s_%(listend)s"
                                  % dict(root=root,
                                         listend=listend)), 'w')
        sglob = sorted(glob(os.path.join(indir, listglob)))        
        lfile.write("\n".join([os.path.basename(x) for x in sglob]))
        lfile.close()
    # aca0
    lfile = open(os.path.join(indir, "%(root)s_tel.lis"
                              % dict(root=root)), 'w')
    for slot in [3, 4, 5, 6, 7, 0, 1, 2]:
        sglob = sorted(glob(os.path.join(indir, 'aca*_%d_*0.fits' % slot)))
        telem_lines = '\n'.join([os.path.basename(x) for x in sglob])
        lfile.write(telem_lines)
        lfile.write("\n")
    lfile.close()
    # pcad adat in outdir if present
    lfile = open(os.path.join(outdir, "%(root)s_dat.lis"
                              % dict(root=root)), 'w')
    sglob = sorted(glob(os.path.join(indir, 'pcad*adat*fits')))
    lfile.write('\n'.join([os.path.basename(x) for x in sglob]))
    lfile.close()


def get_range_ai(ai_cmds, proc_range):
    """
    Limit the array of "todo" aspect to those requested in the --range
    specifier
    """
    # if no --range, return all
    if not proc_range:
        return ai_cmds
    # if a single integer, return on that aspect interval
    intmatch = re.match('^(\d+)$', proc_range)
    if intmatch:
        interv = int(intmatch.group(1))
        return [ai_cmds[int(intmatch.group(1))]]
    # if of the form 0:1, return that range of intervals
    # (python form, not inclusive)
    imatch = re.match('^(\d+):(\d+)$', proc_range)
    if imatch:
        return ai_cmds[int(imatch.group(1)):int(imatch.group(2))]
    # if of the form 1: , return range 1 -> end
    omatch = re.match('^(\d+):$', proc_range)
    if omatch:
        return ai_cmds[int(omatch.group(1)):]
    # if of the form 0:+3000, find a tstop corresponding
    # to tstart of aspect interval 0 plus 3000 seconds
    tmatch = re.match('^(\d+):\+(\d+)$', proc_range)
    if tmatch:
        # get n seconds of specified interval
        interv = int(tmatch.group(1))
        seconds = int(tmatch.group(2))
        tstart = ai_cmds[interv]['istart']
        tstop = tstart + seconds
        cut_ai = ai_cmds[interv].copy()
        cut_ai['istop'] = tstop
        print "attempted to process first %d sec of ai %d" % (
            seconds, interv)
        return [cut_ai]


def queue_ai(ais, skip_slot=None):
    """
    Write JSON dictionaries for each job out to a directory.
    """
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)
    for ai in ais:
        # make a random string for a job
        id = ''.join(random.choice(string.ascii_letters) for x in range(6))
        job_id = "%s_%d.job" % (id, int(ai['istart']))
        f = open(os.path.join(job_dir, job_id), 'w')
        f.write(json.dumps(ai, sort_keys=True, indent=4))
        f.close()

#def make_pipe_cmds(ais, skip_slot=None):
#    for cmd in ais:
#        print 'flt_run_pipe -i %(indir)s -o %(outdir)s ' % cmd \
#            + '-r %(root)s -t %(pipe_ped)s ' % cmd \
#            + '-a "INTERVAL_START=%(istart)f" ' % cmd \
#            + '-a "INTERVAL_STOP=%(istop)f" ' % cmd \
#            + '-a obiroot=%(obiroot)s 2>&1 >> %(log)s' % cmd


def main(opt):
    job_dir = opt.jobdir
    # get files
    if (opt.obsid or opt.archive):
        arc5 = Ska.arc5gl.Arc5gl()
        for (prod, query) in pipe_config['archfiles']:
            proddir = os.path.join(opt.dir, prod)
            if not os.path.exists(proddir):
                os.makedirs(proddir)
            else:
                #print "%s exists skipping..." % prod
                continue
            obsid = opt.obsid or opt.archive
            #arc5.echo = 1
            arc5.sendline("cd %s" % os.path.abspath(proddir))
            arc5.sendline("obsid=%d" % int(obsid))
            arc5.sendline("get %s" % query)

    # check files
    for filetype in ['infiles', 'outfiles']:
        for fileglob in pipe_config[filetype]:
            match = glob(os.path.join(opt.dir, fileglob))
            if not len(match):
                raise ValueError("No files found for glob %s"
                                 % fileglob)
            for mfile in match:
                if re.match(".*\.gz", mfile):
                    bash("gunzip -f %s" % os.path.abspath(mfile))

    # constant aspect interval files
    obi = {}
    caiprops_files = glob(os.path.join(opt.dir, "asp05",
                                       "pcad*cai0a*"))
    for cai_file in caiprops_files:
        cai = get_par(cai_file, cai_override)
        if not cai['obi_num'] in obi:
            obi[cai['obi_num']] = {}
        interval = 0
        while ("istart_%d" % interval in cai
               and "istop_%d" % interval in cai):
            obi[cai['obi_num']][interval] = \
                {'istart': cai['istart_%d' % interval],
                 'istop':  cai['istop_%d' % interval]}
            interval += 1

    ai_cmds = []
    for obi_num in obi:
        # read possible obspars
        obspar_files = glob(os.path.join(opt.dir, "obspar/*.par"))
        for ofile in obspar_files:
            obspar = get_obspar(ofile)
            if obspar['obi_num'] == obi_num:
                obsmatch = re.search('axaf(.+)_obs0a.par', ofile)
                obiroot = obsmatch.group(1)
        if not obiroot:
            raise ValueError("no obspar for obi %d" % obi_num)

        for ai_num in obi[obi_num]:
            aspect_interval = obi[obi_num][ai_num]
            istart = aspect_interval['istart']
            istop = aspect_interval['istop']
            root = "f%09d" % istart
            # directory setup
            workdir, indir, outdir = dir_setup(opt.dir,
                                               int(istart))

            # if skipping the slot by chucking the telem
            telem_skip_slot = []
            process_skip_slot = []
            if opt.skip_slot_method == 'telem':
                telem_skip_slot = opt.skip_slot
            else:
                process_skip_slot = opt.skip_slot

            # link relevant files
            link_files(opt.dir, indir, outdir, istart, istop,
                       obiroot, telem_skip_slot)

            # make list files
            make_list_files(opt.dir, indir, outdir, root)

            # spec
            cmd = dict(dir=os.path.abspath(opt.dir),
                       obi=obi_num,
                       indir=indir,
                       outdir=outdir,
                       root=root,
                       pipe_ped="%s.ped" % pipe_config['pipe_ped'],
                       istart=istart,
                       istop=istop,
                       obiroot=obiroot,
                       log="%s_f%09d.log" % (pipe_config['pipe_ped'], istart))
            if len(process_skip_slot):
                cmd['skip_slot'] = process_skip_slot

            ai_cmds.append(cmd)

    range_ais = get_range_ai(ai_cmds, opt.range)
    queue_ai(range_ais)

if __name__ == '__main__':
    opt, args = get_options()
    main(opt)
