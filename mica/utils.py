# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

from kadi.commands.core import ska_load_dir

DEFAULT_CONFIG = dict(
    timelines_db=dict(
        dbi="sqlite",
        server=os.path.join(os.environ["SKA"], "data", "cmd_states", "cmd_states.db3"),
    )
)


def load_name_to_mp_dir(load_name):
    """Convert ``load_name`` like DEC2506C to /2006/DEC2506/oflsc/.

    :param load_name: str load name
    :returns: str mica-format mission planning dir
    """
    # Get the last 3 parts of the full load directory path YEAR/LOAD/ofls{REV}
    dir_parts = ska_load_dir(load_name).parts
    out = "/" + "/".join(dir_parts[-3:]) + "/"
    return out
