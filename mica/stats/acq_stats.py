import os
import tables

SKA = os.environ.get('SKA', '/proj/sot/ska')
TABLE_FILE = os.path.join(SKA, 'data', 'acq_stats', 'acq_stats.h5')


def get_stats(filter=True):
    """
    Retrieve numpy array of acq stats

    :param filter: True filters out 'known_bad' rows from the table
    :returns acq_stats: numpy.ndarray
    """

    h5 = tables.open_file(TABLE_FILE, 'r')
    stats = h5.root.data[:]
    h5.close()
    if filter:
        stats = stats[~stats['known_bad']]
    return stats
