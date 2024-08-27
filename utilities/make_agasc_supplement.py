# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Create the initial version of agasc_supplement.h5.

This file is a supplement to the stable AGASC to inform star selection
and star catalog checking.

This script simply creates the initial file that has only bad
stars from two sources:
 - starcheck bad star list
   https://github.com/sot/starcheck/blob/master/starcheck/data/agasc.bad
 - GAIA high proper motion file $SKA/analysis/gaia/agasc_gaia_xmatch_PM_gt_50mas.fits.gz
   See: https://nbviewer.jupyter.org/url/cxc.cfa.harvard.edu/mta/ASPECT/ipynb/star_selection/gaia
        GAIA guide star crossmatch.ipynb

"""

import os
from pathlib import Path

import numpy as np
from astropy.table import Table

HOME = Path(os.environ["HOME"])
SKA = Path(os.environ["SKA"])

agasc_ids = []
sources = []

# Starcheck bad star list is not installed anywhere so just grab from local git repo
lines = open(
    HOME / "git" / "starcheck" / "starcheck" / "data" / "agasc.bad", "r"
).readlines()
for line in lines:
    line = line.strip()
    if line.startswith("#"):
        continue
    agasc_ids.append(line.split()[0])
    sources.append(1)  # source=1 implies this is from the starcheck agasc.bad file

# GAIA
dat = Table.read(SKA / "analysis" / "gaia" / "agasc_gaia_xmatch_PM_gt_50mas.fits.gz")
agasc_ids.extend(dat["AGASC_ID"].tolist())
sources.extend([2] * len(dat))

agasc_ids = np.array(agasc_ids, dtype=np.int32)
sources = np.array(sources, dtype=np.int16)

out = Table([agasc_ids, sources], names=["agasc_id", "source"])
out.write("agasc_supplement.h5", format="hdf5", path="bads")
