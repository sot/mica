# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Experimental/alpha code to work with ACA L0 Header 3 data.

This module provides tools for reading and processing ACA (Aspect Camera Assembly)
Level 0 Header 3 telemetry data.
"""

import functools
import warnings

import numpy as np
from cxotime import CxoTime, CxoTimeLike
from scipy.interpolate import interp1d

from mica.archive import aca_l0

TWO_TO_15 = np.uint16(2**15)


def two_byte_sum(byte_msids, scale=1, as_readout_offset=False) -> callable:
    """
    Create a function to combine two bytes into a 16-bit signed integer.

    Parameters
    ----------
    byte_msids : list of str
        List of two MSID names representing the byte values to combine.
    scale : float, optional
        Scale factor to apply to the result. Default is 1.
    as_readout_offset : bool, optional
        If True, apply special readout offset formula. Default is False.

    Returns
    -------
    callable
        Function that takes slot_data and returns combined 16-bit values.
    """

    def func(slot_data) -> np.ndarray:
        # For each pair bytes0[i], bytes1[i], return the 16-bit signed integer
        # corresponding to those two bytes. The input bytes are unsigned.
        bytes0 = slot_data[byte_msids[0]].astype(np.uint8)
        bytes1 = slot_data[byte_msids[1]].astype(np.uint8)

        # Make a 2xN array, then transpose to Nx2, then flatten to 2N, then copy to
        # get values continous in memory.
        bytes8_2xN = np.vstack([bytes0, bytes1], dtype=np.uint8)
        bytes8 = bytes8_2xN.transpose().flatten().copy()

        if as_readout_offset:
            # 16-bit readout offset values need this magic formula from M. Baski (see
            # PEA sampled background patch notes). Note that scale is ignored here.
            uints16 = bytes8.view(">u2")
            out = (uints16 + TWO_TO_15).view("<i2")
        else:
            # Now view the 2N bytes as N 16-bit signed integers.
            ints16 = bytes8.view(">i2")
            out = ints16 * scale

        return out

    return func


# 8x8 header values (largest possible ACA0 header set)
ACA_DTYPE = [
    ("TIME", ">f8"),
    ("HD3TLM62", "|u1"),
    ("HD3TLM63", "|u1"),
    ("HD3TLM64", "|u1"),
    ("HD3TLM65", "|u1"),
    ("HD3TLM66", "|u1"),
    ("HD3TLM67", "|u1"),
    ("HD3TLM72", "|u1"),
    ("HD3TLM73", "|u1"),
    ("HD3TLM74", "|u1"),
    ("HD3TLM75", "|u1"),
    ("HD3TLM76", "|u1"),
    ("HD3TLM77", "|u1"),
]

ACA_DTYPE_NAMES = [k[0] for k in ACA_DTYPE]


a_to_d = [
    ("3D70", -40),
    ("3C9E", -35),
    ("3B98", -30),
    ("3A55", -25),
    ("38CD", -20),
    ("36FD", -15),
    ("34E1", -10),
    ("3279", -5),
    ("2FCB", 0),
    ("2CDF", 5),
    ("29C5", 10),
    ("2688", 15),
    ("2340", 20),
    ("2000", 25),
    ("1CD3", 30),
    ("19CC", 35),
    ("16F4", 40),
    ("1454", 45),
    ("11EF", 50),
    ("0FC5", 55),
    ("0DD8", 60),
    ("0C22", 65),
    ("0A9F", 70),
    ("094D", 75),
    ("0825", 80),
]

# reverse this end up with increasing hex values
a_to_d = a_to_d[::-1]
a_to_d = np.rec.fromrecords(a_to_d, names=["hex", "tempC"])
x = np.array([int(a, 16) for a in a_to_d["hex"]])
ad_func = interp1d(x, a_to_d["tempC"], kind="cubic", bounds_error=False)


def ad_temp(msids):
    """
    Create a function to convert A/D readings to temperature.

    Parameters
    ----------
    msids : list of str
        List of MSID names for the A/D readings.

    Returns
    -------
    callable
        Function that takes slot_data and returns temperature values.
    """

    def func(slot_data):
        sum = two_byte_sum(msids)(slot_data)
        out = ad_func(sum)
        return out

    return func


# Define a dictionary that defines the header 3 'MSID's. This includes a value key that
# describes how to determine the value of the MSID.
#
# For an MSID of HD3TLM<I><W> in ACA image data for slot <S>, that maps to table 11.1
# like:
# Image No.   = <S>
# Image type  = <I>
# Hdr 3 Word  = <W>

HDR3_DEF = {
    "062": {
        "desc": "AC status word",
        "msid": "ac_status_word",
        "longdesc": """
AC status word.  A status word read from the AC.  The bits in the word are
defined as follows:

xxxx xxrr rrtt eccp

X = spare bits
R = CCD readout mode
   1 => S/H input is grounded during pixel readout
   2 => CCD reset is pulsed during column flush
   4 => CCD reset/sw is pulsed during row shifts
   8 => MUX is switched to ground when A/D not in use
T = Test signal select (test signals not available in flight)
E = Cmd error.  AC cmd input buffer was overwritten
C = Clock period for parallel shifts
   0 => 6 microsec
   1 => 12 microsec
   2 => 24 microsec
   3 => 48 microsec
P = PromWrite? flag: true when AC EEPROM is in the write mode.
""",
    },
    "064": {
        "desc": "Misc status bits",
        "msid": "misc_status_bits",
        "longdesc": """
Miscellaneous status bits showing the status of the following 16 flag variables
starting with the LSB and ending with the MSB:

bit 0 (LSB): AcSendTimeOut?
bit 1: AcIdleTimeOut?
bit 2: TecActive?
bit 3: TecHeat?
bit 4: DecAcTable?
bit 5: AcTableCkSumOK?
bit 6: StackError?
bit 7: WarmBoot?
bit 8: IdleCode LSB
bit 9: CalMode?
bit 10: CalModePending?
bit 11: IuData?
bit 12: IuDataPending?
bit 13: DsnFixed?
bit 14: InitialCalFillOK?
bit 15 (MSB): IoUpdTimeout?
""",
    },
    "066": {
        "desc": "A/D CCD molyb therm 1",
        "msid": "ccd_molyb_therm_1",
        "value": ad_temp(["HD3TLM66", "HD3TLM67"]),
        "longdesc": """
A/D converter reading for the CCD moly base thermistor number 1
""",
    },
    "072": {
        "desc": "A/D CCD molyb therm 2",
        "msid": "ccd_molyb_therm_2",
        "value": ad_temp(["HD3TLM72", "HD3TLM73"]),
        "longdesc": """
A/D converter reading for the CCD moly base thermistor number 2
""",
    },
    "074": {
        "desc": "A/D CCD detector therm",
        "msid": "ccd_det_therm",
        "value": ad_temp(["HD3TLM74", "HD3TLM75"]),
        "longdesc": """
A/D converter reading for the CCD detector thermistor
""",
    },
    "076": {
        "desc": "A/D +5 volt PS (V)",
        "msid": "ad_5v_ps",
        "value": two_byte_sum(["HD3TLM76", "HD3TLM77"], scale=0.20518e-3),
        "longdesc": """
A/D converter reading for the +5 volt power supply; 1 LSB=0.30518 mv
""",
    },
    "162": {
        "desc": "A/D +15 volt PS (V)",
        "msid": "ad_15v_ps",
        "value": two_byte_sum(["HD3TLM62", "HD3TLM63"], scale=0.61035e-3),
        "longdesc": """
A/D converter reading for the +15 volt power supply; 1 LSB=0.61035 mv
""",
    },
    "164": {
        "desc": "A/D -15 volt PS (V)",
        "msid": "ad_m15v_ps",
        "value": two_byte_sum(["HD3TLM64", "HD3TLM65"], scale=0.61035e-3),
        "longdesc": """
A/D converter reading for the -15 volt power supply; 1 LSB=0.61035 mv
""",
    },
    "166": {
        "desc": "A/D +27 volt PS (V)",
        "msid": "ad_27v_ps",
        "value": two_byte_sum(["HD3TLM66", "HD3TLM67"], scale=1.04597e-3),
        "longdesc": """
A/D converter reading for the +27 volt power supply; 1 LSB=1.04597 mv
""",
    },
    "172": {
        "desc": "A/D analog ground (V)",
        "msid": "ad_analog_gnd",
        "value": two_byte_sum(["HD3TLM72", "HD3TLM73"], scale=0.30518e-3),
        "longdesc": """
A/D converter reading for analog ground; 1 LSB=0.30518 mv
""",
    },
    "174": {
        "desc": "A/D for A/D convertor therm",
        "msid": "ad_converter_therm",
        "value": ad_temp(["HD3TLM74", "HD3TLM75"]),
        "longdesc": """
A/D converter reading for the A/D converter thermistor.
""",
    },
    "176": {
        "desc": "A/D secondary mirror therm. HRMA side",
        "msid": "ad_smhs_therm",
        "value": ad_temp(["HD3TLM76", "HD3TLM77"]),
        "longdesc": """
A/D converter reading for the secondary mirror thermistor, HRMA side
""",
    },
    "262": {
        "desc": "A/D secondary mirror therm. Opp HRMA side",
        "msid": "ad_smohs_therm",
        "value": ad_temp(["HD3TLM62", "HD3TLM63"]),
        "longdesc": """
A/D converter reading for the secondary mirror thermistor, Opposite from the
HRMA side
""",
    },
    "264": {
        "desc": "A/D primary mirror therm. HRMA side",
        "msid": "ad_pmhs_therm",
        "value": ad_temp(["HD3TLM64", "HD3TLM65"]),
        "longdesc": """
A/D converter reading for the primary mirror thermistor, HRMA side
""",
    },
    "266": {
        "desc": "A/D primary mirror therm. Opp HRMA side",
        "msid": "ad_pmohs_therm",
        "value": ad_temp(["HD3TLM66", "HD3TLM67"]),
        "longdesc": """
A/D converter reading for the primary mirror thermistor, opposite from the
HRMA side
""",
    },
    "272": {
        "desc": "A/D AC housing therm.  HRMA side",
        "msid": "ad_achhs_therm",
        "value": ad_temp(["HD3TLM72", "HD3TLM73"]),
        "longdesc": """
A/D converter reading for the AC housing thermistor, HRMA side
""",
    },
    "274": {
        "desc": "A/D AC housing therm.  Opp HRMA side",
        "msid": "ad_achohs_therm",
        "value": ad_temp(["HD3TLM74", "HD3TLM75"]),
        "longdesc": """
A/D converter reading for the AC housing thermistor, opposite HRMA side
""",
    },
    "276": {
        "desc": "A/D lens cell therm.",
        "msid": "ad_lc_therm",
        "value": ad_temp(["HD3TLM76", "HD3TLM77"]),
        "longdesc": """
A/D converter reading for the lens cell thermistor
""",
    },
    "362": {
        "desc": "Processor stack pointer and telem update counter",
        "msid": "proc_stack_telem_ctr",
        "longdesc": """
A word containing the processor data stack pointer in the high byte, and
an update counter in the low byte that increments once for every 1.025
second telemetry update.
""",
    },
    "364": {
        "desc": "Science header pulse period",
        "msid": "sci_hdr_pulse_period",
        "longdesc": """
The science header pulse period, as measured by the PEA; 1 LSB = 2 microseconds
""",
        "nbytes": 4,
    },
    "372": {
        "desc": "16-bit zero offset for pixels from CCD quad A",
        "msid": "zero_off16_quad_a",
        "value": two_byte_sum(["HD3TLM72", "HD3TLM73"], as_readout_offset=True),
        "longdesc": """
A 16-bit zero offset for pixels read from CCD quadrant A; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
""",
    },
    "374": {
        "desc": "16-bit zero offset for pixels from CCD quad B",
        "msid": "zero_off16_quad_b",
        "value": two_byte_sum(["HD3TLM74", "HD3TLM75"], as_readout_offset=True),
        "longdesc": """
A 16-bit zero offset for pixels read from CCD quadrant B; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
""",
    },
    "376": {
        "desc": "16-bit zero offset for pixels from CCD quad C",
        "msid": "zero_off16_quad_c",
        "value": two_byte_sum(["HD3TLM76", "HD3TLM77"], as_readout_offset=True),
        "longdesc": """
A 16-bit zero offset for pixels read from CCD quadrant C; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
""",
    },
    "462": {
        "desc": "16-bit zero offset for pixels from CCD quad D",
        "msid": "zero_off16_quad_d",
        "value": two_byte_sum(["HD3TLM62", "HD3TLM63"], as_readout_offset=True),
        "longdesc": """
A 16-bit zero offset for pixels read from CCD quadrant D; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
""",
    },
    "464": {
        "desc": "32-bit zero offset for pixels from CCD quad A",
        "msid": "zero_off32_quad_a",
        "longdesc": """
A 32-bit zero offset for pixels read from CCD quadrant A; 1 LSB = 2^-16
A/D converter counts
""",
        "nbytes": 4,
    },
    "472": {
        "desc": "32-bit zero offset for pixels from CCD quad B",
        "msid": "zero_off32_quad_b",
        "longdesc": """
A 32-bit zero offset for pixels read from CCD quadrant B; 1 LSB = 2^-16
A/D converter counts
""",
        "nbytes": 4,
    },
    "476": {
        "desc": "32-bit zero offset for pixels from CCD quad C",
        "msid": "zero_off32_quad_c",
        "longdesc": """
A 32-bit zero offset for pixels read from CCD quadrant C; 1 LSB = 2^-16
A/D converter counts
""",
        "nbytes": 4,
    },
    "564": {
        "desc": "32-bit zero offset for pixels from CCD quad D",
        "msid": "zero_off32_quad_d",
        "longdesc": """
A 32-bit zero offset for pixels read from CCD quadrant D; 1 LSB = 2^-16
A/D converter counts
""",
        "nbytes": 4,
    },
    "572": {
        "desc": "last CCD flush duration",
        "msid": "ccd_flush_dur",
        "longdesc": """
The time required for the most recent flush of the CCD; 1 LSB=2 microseconds
""",
    },
    "574": {
        "desc": "CCD row shift clock period",
        "msid": "ccd_row_shift_period",
        "longdesc": """
The CCD row shift clock period currently in effect; 1 LSB = 1 microsecond
""",
    },
    "576": {
        "desc": "average backround reading",
        "msid": "avg_bkg",
        "longdesc": """
An overall average background reading derived from the most recent CCD readout.
This is an average from all tracked images and from all search readout
segments. One LSB = 1 A/D converter count (nominally 5 electrons).
""",
        "value": two_byte_sum(["HD3TLM76", "HD3TLM77"], scale=5),
    },
    "662": {
        "desc": "Header 1 for DSN records",
        "msid": "dsn_hdr1",
        "longdesc": """
Header 1 for Deep Space Network record.
""",
        "nbytes": 6,
    },
    "672": {
        "desc": "record counter and Header 2 for DSN",
        "msid": "dsn_hdr2",
        "longdesc": """
The record counter and Header 2 for Deep Space Network records.
The record counter occupies the three high order bytes, and Header 2
occupies the low order byte.
""",
        "nbytes": 4,
    },
    "676": {
        "desc": "CCD temperature",
        "msid": "ccd_temp",
        "longdesc": """
CCD temperature.  1 LSB=0.01/(2^16) degrees C.  The high order 16 bits give the
CCD temperature in units of 1 LSB = 0.01 degrees C.
""",
        "nbytes": 4,
        "value": two_byte_sum(["HD3TLM76", "HD3TLM77"], scale=0.01),
    },
    "764": {
        "desc": "CCD setpoint",
        "msid": "ccd_setpoint",
        "longdesc": """
The CCD temperature control setpoint; 1 LSB=0.01 degrees C
""",
        "value": two_byte_sum(["HD3TLM64", "HD3TLM65"], scale=0.01),
    },
    "766": {
        "desc": "temperature for position/angle cal",
        "msid": "aca_temp",
        "longdesc": """
The temperature used in the angle calibration equations that convert star
positions from CCD row and column coordinates to Y and Z angles for OBC
telemetry; 1 LSB = 1/256 degrees C.
""",
        "nbytes": 4,
        "value": two_byte_sum(["HD3TLM72", "HD3TLM73"], scale=1 / 256.0),
    },
    "774": {
        "desc": "RAM address of last write-read test failure",
        "msid": "last_ram_fail_addr",
        "longdesc": """
The address in RAM of the failure most recently detected by the RAM
write-and-read test
""",
    },
    "776": {
        "desc": "TEC DAC number",
        "msid": "dac",
        "value": two_byte_sum(["HD3TLM76", "HD3TLM77"]),
        "longdesc": """
The number most recently written to the TEC power control DAC.
""",
    },
}

MSID_DEFS = {
    value["msid"]: value | {"slot": key[0]} for (key, value) in HDR3_DEF.items()
}


@functools.lru_cache(maxsize=8)
def get_hdr3_slot_data(tstart: float, tstop: float, slot: int):
    """
    Get data for HDR3 telemetry for slot and time range.

    This allows efficient data access for an MSIDset query that fetches multiple MSIDs
    (slots) over the same time range.

    Parameters
    ----------
    tstart : float
        Start time in seconds since epoch.
    tstop : float
        Stop time in seconds since epoch.
    slot : int
        ACA slot number.

    Returns
    -------
    numpy.ndarray
        Array of telemetry data for the specified slot and time range.
    """
    data = aca_l0.get_slot_data(
        tstart,
        tstop + 33,
        slot,
        imgsize=[8],
        columns=ACA_DTYPE_NAMES,
    )
    return data


class MSID(object):
    """
    ACA header 3 data object.

    ACA header 3 data object to work with header 3 data from available 8x8 ACA L0
    telemetry::

      >>> from mica.archive import aca_hdr3
      >>> ccd_temp = aca_hdr3.MSID('ccd_temp', '2012:001', '2012:020')
      >>> type(ccd_temp.vals)
      'numpy.ndarray'

    When given an ``msid`` and ``start`` and ``stop`` range, the object will query the
    ACA L0 archive to populate the object, which includes the MSID values (``vals``) at
    the given times (``times``).

    The parameter ``msid_data`` is used to create an MSID object from the data of
    another MSID object.

    Parameters
    ----------
    msid : str
        MSID name.
    start : CxoTimeLike
        Chandra.Time compatible start time.
    stop : CxoTimeLike
        Chandra.Time compatible stop time.
    filter_bad : deprecated
        This parameter is ignored and will be removed in a future version.
    clear_cache : bool, optional
        Clear the cache for HDR3 slot data after reading the data for this MSID. This is
        useful to set to False when reading multiple MSIDs that are in the same slot and
        time range, to avoid redundant reads of the same data. Default is True.

    Attributes
    ----------
    msid : str
        MSID name.
    vals : numpy.ndarray
        MSID values.
    times : numpy.ndarray
        Time stamps corresponding to the MSID values.
    desc : str
        Short description of the MSID.
    longdesc : str
        Long description of the MSID.
    tstart : float
        Start time in seconds since epoch.
    tstop : float
        Stop time in seconds since epoch.
    datestart : str
        Start date string.
    datestop : str
        Stop date string.
    hdr3_msid : dict
        Header 3 MSID definition dictionary.
    """

    def __init__(self, msid, start, stop, filter_bad=..., clear_cache=True):
        if filter_bad is not ...:
            warnings.warn("'filter_bad' parameter is ignored", UserWarning)

        self.msid: str = msid
        start = CxoTime(start)
        stop = CxoTime(stop)
        self.tstart: float = start.secs
        self.tstop: float = stop.secs
        self.datestart: str = start.date
        self.datestop: str = stop.date

        slot = MSID_DEFS[self.msid]["slot"]

        # Get the 8x8 data with some padding on each end that gets cut later
        slot_data = get_hdr3_slot_data(self.tstart, self.tstop, slot)
        if clear_cache:
            get_hdr3_slot_data.cache_clear()

        # Find samples where the time stamp changes by a value other than 4.1 secs
        # (which is the value for 8x8 readouts).  In that case there must have been a
        # break in L0 decom, typically due to a change to 4x4 or 6x6 data.
        #  t[0] = 1.0
        #  t[1] = 5.1   <= This record could be bad, as indicated by the gap afterward
        #  t[2, 3] = 17.4, 21.5
        # For the diffs add final time stamp of 0.0 so the length matches that of
        # slot_data. The final slot_data record is always chopped but this is OK.
        dt = np.diff(np.concatenate([slot_data["TIME"], [0.0]]))
        ok = np.abs(dt - 4.1) < 1e-3
        slot_data = slot_data[ok]

        # Chop off the padding
        i_stop = np.searchsorted(slot_data["TIME"], self.tstop, side="right")
        slot_data = slot_data[:i_stop]

        # Since we only requested 8x8 image data there should never by any masked
        # values, so convert to normal ndarray (after checking just to be sure).
        slot_data_nomask = np.empty(len(slot_data), dtype=slot_data.dtype)
        for name in slot_data.dtype.names:
            if np.any(slot_data[name].mask):
                raise ValueError(f"unexpected masked values in {name} for {slot}")
            slot_data_nomask[name] = slot_data[name].data

        msid_def = MSID_DEFS[self.msid]
        if "value" not in msid_def:
            raise NotImplementedError(
                f"function to compute {self.msid} from HDR3 telemetry is not defined"
            )

        self.vals = msid_def["value"](slot_data_nomask)
        self.desc = msid_def["desc"]
        self.longdesc = msid_def["longdesc"]
        self.times = slot_data_nomask["TIME"]
        self.hdr3_msid = msid_def

    def copy(self):
        from copy import deepcopy

        return deepcopy(self)

    def plot(self, ax=None, **kwargs):
        import matplotlib.pyplot as plt
        from ska_matplotlib import plot_cxctime

        if ax is None:
            _, ax = plt.subplots()

        plot_cxctime(self.times, self.vals, ax=ax, **kwargs)
        ax.set_title(f"{self.msid.upper()}")


class Msid(MSID):
    """
    ACA header 3 data object (alias for MSID).

    ACA header 3 data object to work with header 3 data from available 8x8 ACA L0
    telemetry. This is an alias for the MSID class.

    Parameters
    ----------
    msid : str
        MSID name.
    start : CxoTimeLike
        Chandra.Time compatible start time.
    stop : CxoTimeLike
        Chandra.Time compatible stop time.

    Notes
    -----
    When given an ``msid`` and ``start`` and ``stop`` range, the object will
    query the ACA L0 archive to populate the object, which includes the MSID
    values (``vals``) at the given times (``times``). Only valid data values
    are returned.

    Examples
    --------
    >>> from mica.archive import aca_hdr3
    >>> ccd_temp = aca_hdr3.Msid('ccd_temp', '2012:001', '2012:020')
    >>> type(ccd_temp.vals)
    <class 'numpy.ndarray'>
    """

    def __init__(self, msid, start, stop):
        super(Msid, self).__init__(msid, start, stop)


class MSIDset(dict):
    """
    ACA header 3 data object for multiple MSIDs.

    ACA header 3 data object to work with header 3 data from
    available 8x8 ACA L0 telemetry. An MSIDset works with multiple
    MSIDs simultaneously.

    Parameters
    ----------
    msids : list of str
        List of MSID names.
    start : CxoTimeLike
        CxoTime compatible start time.
    stop : CxoTimeLike
        CxoTime compatible stop time.

    Attributes
    ----------
    tstart : float
        Start time in seconds since epoch.
    tstop : float
        Stop time in seconds since epoch.
    datestart : str
        Start date string.
    datestop : str
        Stop date string.

    Examples
    --------
    >>> from mica.archive import aca_hdr3
    >>> perigee_data = aca_hdr3.MSIDset(['ccd_temp', 'aca_temp', 'dac'],
    ...                                 '2012:001', '2012:030')
    """

    def __init__(self, msids: list[str], start: CxoTimeLike, stop: CxoTimeLike):
        super(MSIDset, self).__init__()
        start = CxoTime(start)
        stop = CxoTime(stop)
        self.tstart: float = start.secs
        self.tstop: float = stop.secs
        self.datestart: str = start.date
        self.datestop: str = stop.date

        for msid in msids:
            self[msid] = MSID(msid, start, stop, clear_cache=False)

        # Clear cache for memory after getting all the MSIDs
        get_hdr3_slot_data.cache_clear()
