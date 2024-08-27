# script to update obspar and asp_l1 tables with "isdefault" column values
# assumes that the archfiles.db3 files have already been copied to archfiles_redo.db3
# files in the same directories
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import re

import ska_dbi

for t in ["obspar", "asp1"]:
    db = ska_dbi.DBI(dbi="sqlite", server="/data/aca/archive/%s/archfiles_redo.db3" % t)
    db.execute("ALTER table archfiles add column isdefault int")
    if t == "obspar":
        db.execute("ALTER table archfiles add column content text default 'OBSPAR'")
    obsids = db.fetchall("select distinct obsid from archfiles")
    for obs in obsids:
        obsid = obs["obsid"]
        chunk_dir = ("%05d" % obsid)[0:2]
        def_dir = os.path.join("/data/aca/archive", t, chunk_dir, "%05d" % obsid)
        realdefault = os.path.realpath(def_dir)
        lmatch = re.search(r"(\d{5})_v(\d+)$", realdefault)
        if lmatch:
            default_ver = int(lmatch.group(2))
        else:
            print("could not find default for %d" % obsid)
            continue
        print("updating %d default %d" % (obsid, default_ver))
        db.execute(
            """UPDATE archfiles SET isdefault = 1
                  WHERE obsid = %d and revision = %d"""
            % (obsid, default_ver)
        )
        db.execute(
            """UPDATE archfiles SET isdefault = NULL
                  WHERE obsid = %d and revision != %d"""
            % (obsid, default_ver)
        )
