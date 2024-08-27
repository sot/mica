# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import Quaternion
import agasc
from Chandra.Time import DateTime
from chandra_aca.plot import plot_stars

from mica.starcheck import get_starcheck_catalog
from mica.common import MICA_ARCHIVE

CONFIG = dict(
    data_root=os.path.join(MICA_ARCHIVE, 'starcheck'),
    dbi='sqlite',
    server_name='starcheck.db3',
)
CONFIG['server'] = os.path.join(CONFIG['data_root'], CONFIG['server_name'])


def plot(obsid, mp_dir=None):
    sc = get_starcheck_catalog(obsid, mp_dir)
    quat = Quaternion.Quat(
        (sc['obs']['point_ra'], sc['obs']['point_dec'], sc['obs']['point_roll'])
    )
    field = agasc.get_agasc_cone(
        sc['obs']['point_ra'],
        sc['obs']['point_dec'],
        radius=1.5,
        date=DateTime(sc['obs']['mp_starcat_time']).date,
    )
    fig = plot_stars(
        catalog=sc['cat'],
        attitude=quat,
        stars=field,
        title="RA %.2f Dec %.2f" % (sc['obs']['point_ra'], sc['obs']['point_dec']),
    )
    return fig, sc['cat'], sc['obs']
