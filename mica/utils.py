# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

from Chandra.Time import DateTime
from kadi.commands.core import ska_load_dir
import Ska.DBI

DEFAULT_CONFIG = dict(
    timelines_db=dict(
        dbi="sqlite",
        server=os.path.join(os.environ["SKA"], "data", "cmd_states", "cmd_states.db3"),
    )
)


def get_timeline_at_date(date, timelines_db=None):
    """
    Return timeline that contains a given date.  The 'timeline_loads' query is used
    to give the mp_dir as well as datestart and datestop for the timeline.

    :param date: Chandra.Time compatible date
    :param timelines_db: optional already-open handle to cmd_states/timelines database.
    :returns: dictionary of appropriate record from timeline_loads table
    """

    date = DateTime(date).date
    if timelines_db is None:
        timelines_db = Ska.DBI.DBI(**DEFAULT_CONFIG["timelines_db"])
    return timelines_db.fetchone(
        "select * from timeline_loads where datestop >= '%s' "
        " and datestart <= '%s' and scs <= 130 order by datestart desc" % (date, date)
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
