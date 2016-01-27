#!/usr/bin/env

import os
import sys
import re
import subprocess
import json
from itertools import izip, count
import numpy as np
import asciitable
import Ska.Shell
import Ska.DBI


SC1 = ' IDX SLOT        ID  TYPE   SZ  MINMAG    MAG   MAXMAG   YANG   ZANG DIM RES HALFW'
SC2 = ' IDX SLOT        ID  TYPE   SZ  MINMAG    MAG   MAXMAG   YANG   ZANG DIM RES HALFW NOTES'
SC3 = ' IDX SLOT        ID  TYPE   SZ  MINMAG    MAG   MAXMAG   YANG   ZANG DIM RES HALFW PASS NOTES'
SC4 = ' IDX SLOT        ID  TYPE   SZ   P_ACQ    MAG   MAXMAG   YANG   ZANG DIM RES HALFW PASS NOTES'

HDRS = [
    dict(
        pattern=SC1,
        hdrs=['idx', 'slot', 'id', 'type', 'sz', 'minmag', 'mag', 'maxmag',
              'yang', 'zang', 'dim', 'res', 'halfw'],
        col_starts=(1, 4, 7, 19, 25, 30, 38, 46, 54, 61, 68, 72, 76),
        col_ends=(2, 7, 18, 24, 29, 37, 45, 53, 60, 67, 71, 75, 82),
        final_dtype = [('idx', '<i8'), ('slot', '<i8'), ('idnote', '|S3'), ('id', '<i8'),
                       ('type', '|S3'), ('sz', '|S3'),
                       ('minmag', '<f8'), ('mag', '<f8'), ('maxmag', '<f8'),
                       ('yang', '<i8'), ('zang', '<i8'), ('dim', '<i8'),
                       ('res', '<i8'), ('halfw', '<i8')]),
    dict(
        pattern=SC2,
        hdrs=['idx', 'slot', 'id', 'type', 'sz', 'minmag', 'mag', 'maxmag',
              'yang', 'zang', 'dim', 'res', 'halfw', 'pass'],
        # notes explicitly renamed 'pass' for these catalogs
        col_starts=(1, 4, 7, 19, 25, 30, 38, 46, 54, 61, 68, 72, 76, 81),
        col_ends=(2, 7, 18, 24, 29, 37, 45, 53, 60, 67, 71, 75, 80, 92),
        final_dtype = [('idx', '<i8'), ('slot', '<i8'), ('idnote', '|S3'), ('id', '<i8'),
                       ('type', '|S3'), ('sz', '|S3'),
                       ('minmag', '<f8'), ('mag', '<f8'), ('maxmag', '<f8'),
                       ('yang', '<i8'), ('zang', '<i8'), ('dim', '<i8'),
                       ('res', '<i8'), ('halfw', '<i8'), ('pass', '|S10')]),
    dict(
        pattern=SC3,
        hdrs=['idx', 'slot', 'id', 'type', 'sz', 'minmag', 'mag', 'maxmag',
              'yang', 'zang', 'dim', 'res', 'halfw', 'pass', 'notes'],
        col_starts=(1, 4, 7, 19, 25, 30, 38, 46, 54, 61, 68, 72, 76, 81, 87),
        col_ends=(2, 7, 18, 24, 29, 37, 45, 53, 60, 67, 71, 75, 80, 86, 92),
        final_dtype = [('idx', '<i8'), ('slot', '<i8'), ('idnote', '|S3'), ('id', '<i8'),
                       ('type', '|S3'), ('sz', '|S3'),
                       ('minmag', '<f8'), ('mag', '<f8'), ('maxmag', '<f8'),
                       ('yang', '<i8'), ('zang', '<i8'), ('dim', '<i8'),
                       ('res', '<i8'), ('halfw', '<i8'), ('pass', '|S5'), ('notes', '|S5')]),
    dict(
        pattern=SC4,
        hdrs=['idx', 'slot', 'id', 'type', 'sz', 'p_acq', 'mag', 'maxmag',
              'yang', 'zang', 'dim', 'res', 'halfw', 'pass', 'notes'],
        col_starts=(1, 4, 7, 19, 25, 30, 38, 46, 54, 61, 68, 72, 76, 81, 87),
        col_ends=(2, 7, 18, 24, 29, 37, 45, 53, 60, 67, 71, 75, 80, 86, 92),
        final_dtype = [('idx', '<i8'), ('slot', '<i8'), ('idnote', '|S3'), ('id', '<i8'),
                       ('type', '|S3'), ('sz', '|S3'),
                       ('p_acq', '<f8'), ('mag', '<f8'), ('maxmag', '<f8'),
                       ('yang', '<i8'), ('zang', '<i8'), ('dim', '<i8'),
                       ('res', '<i8'), ('halfw', '<i8'), ('pass', '|S5'), ('notes', '|S5')]),
    ]


# Expected type for each of the columns in the catalog table
OKTYPE = dict(idx=int, slot=int, idnote=str, id=int, type=str, sz=str,
              minmag=float, mag=float, maxmag=float,
              yang=int, zang=int,
              dim=int, res=int, halfw=int, notes=str)
OKTYPE['pass'] = str


def get_targ(obs_text):
    targ_search = re.search("(OBSID((.*)\n){2,4}\n)", obs_text)
    if not targ_search:
        raise ValueError("No OBSID found for this catalog")
    for oline in targ_search.group(1).split("\n"):
        short_re = re.match("^OBSID:\s(\S{1,5})\s*$", oline)
        if short_re:
            return {}
        long_re = re.match(
            "OBSID:\s*(\S{1,5})\s+(.*)\s+(\S+)\s+SIM\sZ\soffset:\s*(-*\d+)\s.*\sGrating:\s*(\S+)\s*", 
            oline)
        if long_re:
            form0 = re.match(
                "OBSID:\s*\S{1,5}\s+(.*)\s+(\S+)\s+SIM\sZ\soffset:\s*(-*\d+)\s+Grating:\s*(\S+)\s*",
                oline)
            if form0:
                return dict(target_id=form0.group(1).strip(),
                            sci_instr=form0.group(2),
                            sim_z_offset_steps=int(form0.group(3)),
                            grating=form0.group(4))
            form1 = re.match(
                "OBSID:\s*\S{1,5}\s+(.*)\s+(\S+)\s+SIM\sZ\soffset:\s*(-*\d+)\s+\((-*.+)mm\)\s+Grating:\s*(\S+)\s*",
                oline)
            if form1:
                return dict(target_id=form1.group(1).strip(),
                            sci_instr=form1.group(2),
                            sim_z_offset_steps=int(form1.group(3)),
                            sim_z_offset_mm=float(form1.group(4)),
                            grating=form1.group(5))


def get_dither(obs_text):
    targ_search = re.search("(OBSID((.*)\n){2,4}\n)", obs_text)
    if not targ_search:
        raise ValueError("No OBSID found for this catalog")
    for oline in targ_search.group(1).split("\n"):
        dither_re = re.match(
            "Dither:\s(\S+)\s+Y_amp=\s*(\S+)\s+Z_amp=\s*(\S+)\s+Y_period=\s*(\S+)\s+Z_period=\s*(\S+)\s*",
            oline)
        if dither_re:
            return dict(dither_state=dither_re.group(1),
                        dither_y_amp=float(dither_re.group(2)),
                        dither_z_amp=float(dither_re.group(3)),
                        dither_y_period=float(dither_re.group(4)),
                        dither_z_period=float(dither_re.group(5)))
        dis_dither_re = re.match(
            "Dither:\sOFF\s*",
            oline)
        if dis_dither_re:
            return dict(dither_state='OFF',
                        dither_y_amp=float(0),
                        dither_z_amp=float(0),
                        dither_y_period=None,
                        dither_z_period=None)

    return {}


def get_coords(obs_text):
    coord_search = re.search(
        "RA, Dec, Roll \(deg\):\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s*",
        obs_text)
    if not coord_search:
        return {}
    return dict(point_ra=float(coord_search.group(1)),
                point_dec=float(coord_search.group(2)),
                point_roll=float(coord_search.group(3)))


def get_starcat_header(obs_text):
    starcat_header = re.search(
        "MP_STARCAT\sat\s(\S+)\s\(VCDU\scount\s=\s(\d+)\)",
        obs_text)
    if not starcat_header:
        return {}
    return dict(mp_starcat_time=starcat_header.group(1),
                mp_starcat_vcdu_cnt=int(starcat_header.group(2)))


def get_manvrs(obs_text):
    man_search = re.search("(MP_TARGQUAT((.*)\n){2,}\n)MP_STARCAT", obs_text)
    if not man_search:
        return {}
    manvr_block = man_search.group(1)
    manvrs = []
    for mline in manvr_block.split("\n"):
        targquat_re = re.match(
            "MP_TARGQUAT at (\S+) \(VCDU count = (\d+)\).*",
            mline)
        if targquat_re:
            curr_manvr = dict(mp_targquat_time=targquat_re.group(1),
                              mp_targquat_vcdu_cnt=int(targquat_re.group(2)))
        quat_re = re.match(
            "\s+Q1,Q2,Q3,Q4:\s+(-?\d\.\d+)\s+(-?\d\.\d+)\s+(-?\d\.\d+)\s+(-?\d\.\d+).*",
            mline)
        if quat_re:
            curr_manvr.update(dict(target_Q1=float(quat_re.group(1)),
                                   target_Q2=float(quat_re.group(2)),
                                   target_Q3=float(quat_re.group(3)),
                                   target_Q4=float(quat_re.group(4))))
        angle_re = re.match(
            "\s+MANVR: Angle=\s+(\d+\.\d+)\sdeg\s+Duration=\s+(\d+)\ssec\s+Slew\serr=\s+(\d+\.\d)\sarcsec",
            mline)
        if angle_re:
            curr_manvr.update(dict(angle_deg=float(angle_re.group(1)),
                                   duration_sec=int(angle_re.group(2)),
                                   slew_err_arcsec=float(angle_re.group(3)),
                                   instance=len(manvrs)))
            if 'target_Q1' not in curr_manvr:
                raise ValueError("No Q1,Q2,Q3,Q4 line found when parsing starcheck.txt")
            manvrs.append(curr_manvr)
            del curr_manvr
    return manvrs


def get_catalog(obs_text):
    catmatch = re.search("MP_STARCAT", obs_text)
    if not catmatch:
        return {}
    # hdr line is bracked by two lines of dashes
    hdrmatch = re.search("-{20,}\n(.*)\n-{20,}", obs_text)
    hdr = hdrmatch.group(1)
    for posshdr in HDRS:
        if posshdr['pattern'] == hdr:
            hdrformat = posshdr
    if hdrformat is None:
        raise ValueError
    # get the lines that start with an [
    catlines = [t for t in obs_text.split("\n")
                if re.compile("^\[.*").match(t)]
    rawcat = asciitable.read(catlines, Reader=asciitable.FixedWidthNoHeader,
                             col_starts=hdrformat['col_starts'],
                             col_ends=hdrformat['col_ends'],
                             names=hdrformat['hdrs'])
    cat = []
    for row, idx in izip(rawcat, count()):
        catrow = dict()
        for field in row.dtype.names:
            if field not in OKTYPE.keys():
                continue
            if field == 'id':
                idmatch = re.match("(\D+)?(\d+)?", str(row[field]))
                catrow['idnote'] = idmatch.group(1)
                catrow['id'] = idmatch.group(2)
                if catrow['id'] is not None:
                    catrow['id'] = int(catrow['id'])
                continue
            if row[field] == '---':
                catrow[field] = None
                continue

            catrow[field] = OKTYPE[field](row[field])
            #if not isinstance(row[field], oktype[field]):
            #    raise TypeError("%s not %s" % (field, oktype[field]))
            #catrow[field] = row[field]
        cat.append(catrow)
    return cat


def get_warnings(obs_text):
    warnlines = [t for t in obs_text.split("\n")
                 if re.compile("^\>\>\s+(WARNING|INFO).*").match(t)]
    warn = []
    for wline in warnlines:
        form0 = re.match(
            '^\>\>\s+(WARNING|INFO)\s*:\s+(.+)\.\s+(\[\s?(\d+)\]\-\[\s?(\d+)\]).?\s*(.*)',
            wline)
        if form0:
            # append two warnings for this format
            warn.append(dict(warning_type=form0.group(2),
                             idx=form0.group(4),
                             warning="%s -> %s" % (form0.group(3), form0.group(6))))
            warn.append(dict(warning_type=form0.group(2),
                             idx=form0.group(5),
                             warning="%s -> %s" % (form0.group(3), form0.group(6))))
            continue
        form1 = re.match('.*(WARNING|INFO).*\[\s?(\d+)\]([\w\s]+)\.(.*)$',
                         wline)
        if form1:
            warn.append(dict(warning_type=form1.group(3),
                             idx=form1.group(2),
                             warning=form1.group(4)))
            continue
        form2 = re.match(
            '^\>\>\s+(WARNING|INFO)\s*:\s+(.+)\.\s+(?:\[ |\[)(\d+)\].?\s*(.*)',
            wline)
        if form2:
            warn.append(dict(warning_type=form2.group(2),
                             idx=form2.group(3),
                             warning=form2.group(4)))
            continue
        form3 = re.match(
            '^\>\>\s+(WARNING|INFO)\s*:\s+(.+)',
            wline)
        if form3:
            warn.append(dict(warning_type=None,
                             idx=None,
                             warning=form3.group(2)))
    return warn


def get_pred_temp(obs_text):
    pred_line = re.search(
        "Predicted Max CCD temperature: (-?\d+\.\d)\sC",
        obs_text)
    if not pred_line:
        return None
    return float(pred_line.group(1))


def fix_obs(obs):
    # Fix broken starcheck.txt manually
    # obsid 11866 has non-ascii characters in the starcheck.txt file
    # (and python complains if I copy them here for reference)
    if obs['obsid'] == 11866:
        obs['target_id'] = 'cl0422-5009'


def get_cat(obs_text):
    obsmatch = re.match("^OBSID:\s(\d+).*", obs_text)
    if not obsmatch:
        return {}
    # join the top level stuff into one dictionary
    obs = dict(obsid=int(obsmatch.group(1)))
    obs.update(get_targ(obs_text))
    obs.update(get_coords(obs_text))
    obs.update(get_dither(obs_text))
    obs.update(get_starcat_header(obs_text))
    # the whole catalog should be broken up into
    # the obs/manvrs/catalog/warnings keys
    fix_obs(obs)
    return dict(
        obsid=obs['obsid'],
        obs=obs,
        manvrs=get_manvrs(obs_text),
        catalog=get_catalog(obs_text),
        warnings=get_warnings(obs_text),
        pred_ccd_temp=get_pred_temp(obs_text))


def read_starcheck(starcheck_file):
    sc_text = open(starcheck_file, 'r').read()
    chunks = re.split("={20,}\s?\n?\n", sc_text)
    catalogs = []
    for chunk, idx in izip(chunks, count()):
        obs = get_cat(chunk)
        if obs:
            catalogs.append(obs)
    return catalogs
