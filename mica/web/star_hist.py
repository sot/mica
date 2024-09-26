# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.table import Table
from Chandra.Time import DateTime

from mica.stats import acq_stats, guide_stats


def get_acq_data(agasc_id):
    """
    Fetch acquisition history from mica acq stats for an agasc id

    :param agasc_id: AGASC id
    :returns: list of dicts of acquisitions
    """

    acq = acq_stats.get_star_stats(agasc_id)
    # make list of dicts for use in light templates in kadi web app
    if not len(acq):
        return []
    acq = Table(acq)
    acq.sort("guide_start")
    acq_table = []
    for s in acq:
        srec = {}
        # Use these columns as they are named from the mica acq stats table
        for col in ["type", "obsid", "obi", "slot", "mag", "mag_obs", "star_tracked"]:
            srec[col] = s[col]
        # rename these columns in the dictionary
        srec["date"] = s["guide_start"]
        srec["acq_dy"] = s["cdy"]
        srec["acq_dz"] = s["cdz"]
        srec["id"] = s["acqid"]
        acq_table.append(srec)
    return acq_table


def get_gui_data(agasc_id):
    """
    Fetch guide/track history for an agasc id

    :param agasc_id: AGASC id
    :returns: list of dicts of uses as guide stars
    """
    gui = guide_stats.get_star_stats(agasc_id)
    if not len(gui):
        return []
    gui = Table(gui)
    gui.sort("kalman_datestart")
    # make list of dicts for use in light templates in kadi web app
    gui_table = []
    for s in gui:
        srec = {}
        # Use these columns as they are named from the mica acq stats table
        for col in ["type", "obsid", "obi", "slot"]:
            srec[col] = s[col]
        # rename these columns in the dictionary
        srec["date"] = s["kalman_datestart"]
        srec["mag"] = s["mag_aca"]
        srec["mag_obs"] = s["aoacmag_mean"]
        srec["perc_not_track"] = (1 - s["f_track"]) * 100.0
        gui_table.append(srec)
    return gui_table


def get_star_stats(agasc_id, start=None, stop=None):
    """
    Fetch acq and gui history of a star

    :param agasc_id: AGASC id
    :param start: start of optional time filter (>=) (Chandra.Time compatible)
    :param stop: stop time of optional time filter (<) (Chandra.time compatible)
    :returns: 2 lists, first of acq attempts, second of guide attempts
    """
    acq_table = get_acq_data(agasc_id)
    gui_table = get_gui_data(agasc_id)
    if start is not None:
        acq_table = [s for s in acq_table if s["date"] >= DateTime(start).date]
        gui_table = [s for s in gui_table if s["date"] >= DateTime(start).date]
    if stop is not None:
        acq_table = [s for s in acq_table if s["date"] < DateTime(stop).date]
        gui_table = [s for s in gui_table if s["date"] < DateTime(stop).date]
    return acq_table, gui_table
