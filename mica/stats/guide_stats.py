import os
import tables

SKA = os.environ.get('SKA', '/proj/sot/ska')
TABLE_FILE = os.path.join(SKA, 'data', 'guide_stats', 'guide_stats.h5')


def get_stats(filter=True):
    """
    Retrieve numpy array of guide stats

    :param filter: True filters out 'known_bad' rows from the table
    :returns gui_stats: numpy.ndarray
    """

    h5 = tables.open_file(TABLE_FILE, 'r')
    stats = h5.root.data[:]
    h5.close()
    if filter:
        stats = stats[~stats['known_bad']]
    return stats


