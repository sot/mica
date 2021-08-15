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
    with tables.open_file(TABLE_FILE, 'r') as h5:
        stats = h5.root.data[:]
        if filter:
            stats = stats[~stats['known_bad']]
    return stats


def get_star_stats(id, filter=True):
    """
    Retrieve guide stats for agasc id

    :param id: agasc id
    :param filter: True filters out rows marked 'known_bad' in table
    :returns stats: numpy.ndarray
    """
    with tables.open_file(TABLE_FILE, 'r') as h5:
        stats = h5.root.data.read_where(f'agasc_id == {id}')
        if filter:
            stats = stats[~stats['known_bad']]
    return stats
