import os
import tables

SKA = os.environ.get("SKA", "/proj/sot/ska")
TABLE_FILE = os.path.join(SKA, "data", "acq_stats", "acq_stats.h5")


def get_stats(filter=True):
    """
    Retrieve numpy array of acq stats

    :param filter: True filters out 'known_bad' rows from the table
    :returns acq_stats: numpy.ndarray
    """

    with tables.open_file(TABLE_FILE, "r") as h5:
        stats = h5.root.data[:]
        if filter:
            stats = stats[~stats["known_bad"]]
    return stats


def get_star_stats(id, filter=True):
    """
    Retrieve acq stats for agasc id

    :param id: agasc id
    :param filter: True filters out rows marked 'known_bad' in table
    :returns acq_stats: numpy.ndarray
    """
    with tables.open_file(TABLE_FILE, "r") as h5:
        stats = h5.root.data.read_where(f"agasc_id == {id}")
        if filter:
            stats = stats[~stats["known_bad"]]
    return stats
