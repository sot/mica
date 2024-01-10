# Licensed under a 3-clause BSD style license - see LICENSE.rst
from parse_cm.paths import load_dir_from_load_name

DEFAULT_CONFIG = {}


def load_name_to_mp_dir(load_name):
    """Convert ``load_name`` like DEC2506C to /2006/DEC2506/oflsc/.

    :param load_name: str load name
    :returns: str mica-format mission planning dir
    """
    # Get the last 3 parts of the full load directory path YEAR/LOAD/ofls{REV}
    dir_parts = load_dir_from_load_name(load_name).parts
    out = "/" + "/".join(dir_parts[-3:]) + "/"
    return out
