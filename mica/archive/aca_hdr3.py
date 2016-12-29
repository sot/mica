"""
Experimental/alpha code to work with ACA L0 Header 3 data
"""
import re
import numpy as np
import numpy.ma as ma
import collections
from scipy.interpolate import interp1d
from Chandra.Time import DateTime
from Ska.Numpy import search_both_sorted

from mica.archive import aca_l0
from mica.common import MissingDataError


#In case it isn't obvious, for an MSID of HD3TLM<I><W> in ACA image data
#for slot <S>, that maps to table 11.1 like:
# Image No.   = <S>
# Image type  = <I>
# Hdr 3 Word  = <W>


def two_byte_sum(byte_msids, scale=1):
    def func(slot_data):
        return ((slot_data[byte_msids[0]].astype('int') >> 7) * (-1 * 65535)
                + (slot_data[byte_msids[0]].astype('int') << 8)
                + (slot_data[byte_msids[1]].astype('int'))) * scale
    return func


# 8x8 header values (largest possible ACA0 header set)
ACA_DTYPE = [('TIME', '>f8'), ('QUALITY', '>i4'), ('IMGSIZE', '>i4'),
             ('HD3TLM62', '|u1'),
             ('HD3TLM63', '|u1'), ('HD3TLM64', '|u1'), ('HD3TLM65', '|u1'),
             ('HD3TLM66', '|u1'), ('HD3TLM67', '|u1'), ('HD3TLM72', '|u1'),
             ('HD3TLM73', '|u1'), ('HD3TLM74', '|u1'), ('HD3TLM75', '|u1'),
             ('HD3TLM76', '|u1'), ('HD3TLM77', '|u1'), ('FILENAME', '<U128')
             ]

ACA_DTYPE_NAMES = [k[0] for k in ACA_DTYPE]


a_to_d = [
('3D70', -40),
('3C9E', -35),
('3B98', -30),
('3A55', -25),
('38CD', -20),
('36FD', -15),
('34E1', -10),
('3279',  -5),
('2FCB',   0),
('2CDF',   5),
('29C5',  10),
('2688',  15),
('2340',  20),
('2000',  25),
('1CD3',  30),
('19CC',  35),
('16F4',  40),
('1454',  45),
('11EF',  50),
('0FC5',  55),
('0DD8',  60),
('0C22',  65),
('0A9F',  70),
('094D',  75),
('0825',  80)]

# reverse this end up with increasing hex values
a_to_d = a_to_d[::-1]
a_to_d = np.rec.fromrecords(a_to_d, names=['hex', 'tempC'])
x = np.array([int(a, 16) for a in a_to_d['hex']])
ad_func = interp1d(x, a_to_d['tempC'], kind='cubic', bounds_error=False)


def ad_temp(msids):

    def func(slot_data):
        sum = two_byte_sum(msids)(slot_data)
        # As of scipy 0.17 cannot interpolate a masked array.  In this
        # case we can temporarily fill with some value that will always
        # be in the range, then re-mask afterward.
        masked = isinstance(sum, np.ma.MaskedArray)
        if masked:
            mask = sum.mask
            sum = sum.filled(16000)
        out = ad_func(sum)
        if masked:
            out = np.ma.MaskedArray(out, mask=mask)
        return out

    return func

# dictionary that defines the header 3 'MSID's.
# Also includes a value key that describes how to determine the value
# of the MSID

HDR3_DEF = {
    '062': {'desc': 'AC status word',
            'msid': 'ac_status_word',
            'longdesc': """
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
"""},

    '064': {'desc': 'Misc status bits',
            'msid': 'misc_status_bits',
            'longdesc': """
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
"""},

    '066': {'desc': 'A/D CCD molyb therm 1',
            'msid': 'ccd_molyb_therm_1',
            'value': ad_temp(['HD3TLM66', 'HD3TLM67']),
            'longdesc': """
A/D converter reading for the CCD moly base thermistor number 1
"""},

    '072': {'desc': 'A/D CCD molyb therm 2',
            'msid': 'ccd_molyb_therm_2',
            'value': ad_temp(['HD3TLM72', 'HD3TLM73']),
            'longdesc': """
A/D converter reading for the CCD moly base thermistor number 2
"""},

    '074': {'desc': 'A/D CCD detector therm',
            'msid': 'ccd_det_therm',
            'value': ad_temp(['HD3TLM74', 'HD3TLM75']),
            'longdesc': """
A/D converter reading for the CCD detector thermistor
"""},

    '076': {'desc': 'A/D +5 volt PS',
            'msid': 'ad_5v_ps',
            'value': two_byte_sum(['HD3TLM76', 'HD3TLM77'], scale=0.20518),
            'longdesc': """
A/D converter reading for the +5 volt power supply; 1 LSB=0.30518 mv
"""},

    '162': {'desc': 'A/D +15 volt PS',
            'msid': 'ad_15v_ps',
            'value': two_byte_sum(['HD3TLM62', 'HD3TLM63'], scale=0.61035),
            'longdesc': """
A/D converter reading for the +15 volt power supply; 1 LSB=0.61035 mv
"""},

    '164': {'desc': 'A/D -15 volt PS',
            'msid': 'ad_m15v_ps',
            'value': two_byte_sum(['HD3TLM64', 'HD3TLM65'], scale=0.61035),
            'longdesc': """
A/D converter reading for the -15 volt power supply; 1 LSB=0.61035 mv
"""},

    '166': {'desc': 'A/D +27 volt PS',
            'msid': 'ad_27v_ps',
            'value': two_byte_sum(['HD3TLM66', 'HD3TLM67'], scale=1.04597),
            'longdesc': """
A/D converter reading for the +27 volt power supply; 1 LSB=1.04597 mv
"""},

    '172': {'desc': 'A/D analog ground',
            'msid': 'ad_analog_gnd',
            'value': ad_temp(['HD3TLM72', 'HD3TLM73']),
            'longdesc': """
A/D converter reading for analog ground; 1 LSB=0.30518 mv
"""},

    '174': {'desc': 'A/D for A/D convertor therm',
            'msid': 'ad_converter_therm',
            'value': ad_temp(['HD3TLM74', 'HD3TLM75']),
            'longdesc': """
A/D converter reading for the A/D converter thermistor.
"""},

    '176': {'desc': 'A/D secondary mirror therm. HRMA side',
            'msid': 'ad_smhs_therm',
            'value': ad_temp(['HD3TLM76', 'HD3TLM77']),
            'longdesc': """
A/D converter reading for the secondary mirror thermistor, HRMA side
"""},

    '262': {'desc': 'A/D secondary mirror therm. Opp HRMA side',
            'msid': 'ad_smohs_therm',
            'value': ad_temp(['HD3TLM62', 'HD3TLM63']),
            'longdesc': """
A/D converter reading for the secondary mirror thermistor, Opposite from the
HRMA side
"""},

    '264': {'desc': 'A/D primary mirror therm. HRMA side',
            'msid': 'ad_pmhs_therm',
            'value': ad_temp(['HD3TLM64', 'HD3TLM65']),
            'longdesc': """
A/D converter reading for the primary mirror thermistor, HRMA side
"""},

    '266': {'desc': 'A/D primary mirror therm. Opp HRMA side',
            'msid': 'ad_pmohs_therm',
            'value': ad_temp(['HD3TLM66', 'HD3TLM67']),
            'longdesc': """
A/D converter reading for the primary mirror thermistor, opposite from the
HRMA side
"""},

    '272': {'desc': 'A/D AC housing therm.  HRMA side',
            'msid': 'ad_achhs_therm',
            'value': ad_temp(['HD3TLM72', 'HD3TLM73']),
            'longdesc': """
A/D converter reading for the AC housing thermistor, HRMA side
"""},

    '274': {'desc': 'A/D AC housing therm.  Opp HRMA side',
            'msid': 'ad_achohs_therm',
            'value': ad_temp(['HD3TLM74', 'HD3TLM75']),
            'longdesc': """
A/D converter reading for the AC housing thermistor, opposite HRMA side
"""},

    '276': {'desc': 'A/D lens cell therm.',
            'msid': 'ad_lc_therm',
            'value': ad_temp(['HD3TLM76', 'HD3TLM77']),
            'longdesc': """
A/D converter reading for the lens cell thermistor
"""},

    '362': {'desc': 'Processor stack pointer and telem update counter',
            'msid': 'proc_stack_telem_ctr',
            'longdesc': """
A word containing the processor data stack pointer in the high byte, and
an update counter in the low byte that increments once for every 1.025
second telemetry update.
"""},

    '364': {'desc': 'Science header pulse period',
            'msid': 'sci_hdr_pulse_period',
            'longdesc': """
The science header pulse period, as measured by the PEA; 1 LSB = 2 microseconds
""",
            'nbytes': 4},

    '372': {'desc': '16-bit zero offset for pixels from CCD quad A',
            'msid': 'zero_off16_quad_a',
            'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant A; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},

    '374': {'desc': '16-bit zero offset for pixels from CCD quad B',
            'msid': 'zero_off16_quad_b',
            'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant B; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},

    '376': {'desc': '16-bit zero offset for pixels from CCD quad C',
            'msid': 'zero_off16_quad_c',
            'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant C; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},

    '462': {'desc': '16-bit zero offset for pixels from CCD quad D',
            'msid': 'zero_off16_quad_d',
            'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant D; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},

    '464': {'desc': '32-bit zero offset for pixels from CCD quad A',
            'msid': 'zero_off32_quad_a',
            'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant A; 1 LSB = 2^-16
A/D converter counts
""",
            'nbytes': 4},

    '472': {'desc': '32-bit zero offset for pixels from CCD quad B',
            'msid': 'zero_off32_quad_b',
            'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant B; 1 LSB = 2^-16
A/D converter counts
""",
            'nbytes': 4},

    '476': {'desc': '32-bit zero offset for pixels from CCD quad C',
            'msid': 'zero_off32_quad_c',
            'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant C; 1 LSB = 2^-16
A/D converter counts
""",
            'nbytes': 4},

    '564': {'desc': '32-bit zero offset for pixels from CCD quad D',
            'msid': 'zero_off32_quad_d',
            'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant D; 1 LSB = 2^-16
A/D converter counts
""",
            'nbytes': 4},

    '572': {'desc': 'last CCD flush duration',
            'msid': 'ccd_flush_dur',
            'longdesc': """
The time required for the most recent flush of the CCD; 1 LSB=2 microseconds
"""},

    '574': {'desc': 'CCD row shift clock period',
            'msid': 'ccd_row_shift_period',
            'longdesc': """
The CCD row shift clock period currently in effect; 1 LSB = 1 microsecond
"""},

    '576': {'desc': 'average backround reading',
            'msid': 'avg_bkg',
            'longdesc': """
An overall average background reading derived from the most recent CCD readout.
This is an average from all tracked images and from all search readout
segments. One LSB = 1 A/D converter count (nominally 5 electrons).
""",
            'value': two_byte_sum(['HD3TLM76', 'HD3TLM77'], scale=5)},

    '662': {'desc': 'Header 1 for DSN records',
            'msid': 'dsn_hdr1',
            'longdesc': """
Header 1 for Deep Space Network record.
""",
            'nbytes': 6},

    '672': {'desc': 'record counter and Header 2 for DSN',
            'msid': 'dsn_hdr2',
            'longdesc': """
The record counter and Header 2 for Deep Space Network records.
The record counter occupies the three high order bytes, and Header 2
occupies the low order byte.
""",
            'nbytes': 4},

    '676': {'desc': 'CCD temperature',
            'msid': 'ccd_temp',
            'longdesc': """
CCD temperature.  1 LSB=0.01/(2^16) degrees C.  The high order 16 bits give the
CCD temperature in units of 1 LSB = 0.01 degrees C.
""",
            'nbytes': 4,
            'value': two_byte_sum(['HD3TLM76', 'HD3TLM77'], scale=0.01)},

    '764': {'desc': 'CCD setpoint',
            'msid': 'ccd_setpoint',
            'longdesc': """
The CCD temperature control setpoint; 1 LSB=0.01 degrees C
""",
            'value': two_byte_sum(['HD3TLM64', 'HD3TLM65'])},

    '766': {'desc': 'temperature for position/angle cal',
            'msid': 'aca_temp',
            'longdesc': """
The temperature used in the angle calibration equations that convert star
positions from CCD row and column coordinates to Y and Z angles for OBC
telemetry; 1 LSB = 1/256 degrees C.
""",
            'nbytes': 4,
            'value': two_byte_sum(['HD3TLM72', 'HD3TLM73'],
                                  scale=1 / 256.)},
    '774': {'desc': 'RAM address of last write-read test failure',
            'msid': 'last_ram_fail_addr',
            'longdesc': """
The address in RAM of the failure most recently detected by the RAM
write-and-read test
"""},

    '776': {'desc': 'TEC DAC number',
            'msid': 'dac',
            'value': two_byte_sum(['HD3TLM76', 'HD3TLM77']),
            'longdesc': """
The number most recently written to the TEC power control DAC.
"""}}

MSID_ALIASES = {HDR3_DEF[key]['msid']: key for key in HDR3_DEF}


class MSID(object):
    """
    ACA header 3 data object to work with header 3 data from
    available 8x8 ACA L0 telemetry::

      >>> from mica.archive import aca_hdr3
      >>> ccd_temp = aca_hdr3.MSID('ccd_temp', '2012:001', '2012:020')
      >>> type(ccd_temp.vals)
      'numpy.ma.core.MaskedArray'

    When given an ``msid`` and ``start`` and ``stop`` range, the object will
    query the ACA L0 archive to populate the object, which includes the MSID
    values (``vals``) at the given times (``times``).

    The parameter ``msid_data`` is used to create an MSID object from
    the data of another MSID object.

    When ``filter_bad`` is supplied then only valid data values are stored
    and the ``vals`` and ``times`` attributes are `np.ndarray` instead of
    `ma.MaskedArray`.

    :param msid: MSID name
    :param start: Chandra.Time compatible start time
    :param stop: Chandra.Time compatible stop time
    :param msid_data: data dictionary or object from another MSID object
    :param filter_bad: remove missing values
    """

    def __init__(self, msid, start, stop, msid_data=None, filter_bad=False):
        if msid_data is None:
            msid_data = MSIDset([msid], start, stop)[msid]
        self.msid = msid
        self.tstart = DateTime(start).secs
        self.tstop = DateTime(stop).secs
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        # msid_data may be dictionary or object with these
        # attributes
        for attr in ('hdr3_msid', 'vals', 'times', 'desc', 'longdesc'):
            if hasattr(msid_data, attr):
                setattr(self, attr, getattr(msid_data, attr))
            else:
                setattr(self, attr, msid_data.get(attr))

        # If requested filter out bad values and set self.bad = None
        if filter_bad:
            self.filter_bad()

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)

    def filter_bad(self, copy=False):
        """Filter out any missing values.

        After applying this method the ``vals`` attributes will be a
        plain np.ndarray object instead of a masked array.

        :param copy: return a copy of MSID object with bad values filtered
        """
        obj = self.copy() if copy else self

        if isinstance(obj.vals, ma.MaskedArray):
            obj.times = obj.times[~obj.vals.mask]
            obj.vals = obj.vals.compressed()

        if copy:
            return obj


class Msid(MSID):
    """
    ACA header 3 data object to work with header 3 data from available 8x8 ACA L0
    telemetry.

    >>> from mica.archive import aca_hdr3
    >>> ccd_temp = aca_hdr3.Msid('ccd_temp', '2012:001', '2012:020')
    >>> type(ccd_temp.vals)
    'numpy.ndarray'

    When given an ``msid`` and ``start`` and ``stop`` range, the object will
    query the ACA L0 archive to populate the object, which includes the MSID
    values (``vals``) at the given times (``times``).  Only valid data values
    are returned.

    :param msid: MSID
    :param start: Chandra.Time compatible start time
    :param stop: Chandra.Time compatible stop time
    """

    def __init__(self, msid, start, stop):
        super(Msid, self).__init__(msid, start, stop, filter_bad=True)


def confirm_msid(req_msid):
    """
    Check to see if the 'MSID' is an alias or is in the HDR3_DEF
    dictionary.  If in the aliases, return the unaliased value.
    :param req_msid: requested msid
    :return: hdr3_def MSID name
    """
    if req_msid in MSID_ALIASES:
        return MSID_ALIASES[req_msid]
    else:
        if req_msid not in HDR3_DEF:
            raise MissingDataError("msid %s not found" % req_msid)
        else:
            return req_msid


def slot_for_msid(msid):
    """
    For a given 'MSID' return the slot number that contains those data.
    """
    mmatch = re.match('(\d)\d\d', msid)
    slot = int(mmatch.group(1))
    return slot


class MSIDset(collections.OrderedDict):
    """
    ACA header 3 data object to work with header 3 data from
    available 8x8 ACA L0 telemetry.  An MSIDset works with multiple
    MSIDs simultaneously.

    >>> from mica.archive import aca_hdr3
    >>> perigee_data = aca_hdr3.MSIDset(['ccd_temp', 'aca_temp', 'dac'],
    ...                                 '2012:001', '2012:030')

    :param msids: list of MSIDs
    :param start: Chandra.Time compatible start time
    :param stop: Chandra.Time compatible stop time
    """
    def __init__(self, msids, start, stop):
        super(MSIDset, self).__init__()
        self.tstart = DateTime(start).secs
        self.tstop = DateTime(stop).secs
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        slot_datas = {}
        slots = set(slot_for_msid(confirm_msid(msid)) for msid in msids)
        for slot in slots:
            # get the 8x8 data
            tstop = self.tstop + 33.0  # Major frame of padding
            slot_data = aca_l0.get_slot_data(
                self.tstart, tstop, slot,
                imgsize=[8], columns=ACA_DTYPE_NAMES)

            # Find samples where the time stamp changes by a value other than 4.1 secs
            # (which is the value for 8x8 readouts).  In that case there must have been a
            # break in L0 decom, typically due to a change to 4x4 or 6x6 data.
            #  t[0] = 1.0
            #  t[1] = 5.1   <= This record could be bad, as indicated by the gap afterward
            #  t[2, 3] = 17.4, 21.5
            # To form the time diffs first add `tstop` to the end so that if 8x8 data
            # does not extend through `tstop` then the last record gets chopped.
            dt = np.diff(np.concatenate([slot_data['TIME'], [tstop]]))
            bad = np.abs(dt - 4.1) > 1e-3
            slot_data[bad] = ma.masked

            # Chop off the padding
            i_stop = np.searchsorted(slot_data['TIME'], self.tstop, side='right')
            slot_data = slot_data[:i_stop]

            # explicitly unmask useful columns
            slot_data['TIME'].mask = ma.nomask
            slot_data['IMGSIZE'].mask = ma.nomask
            slot_data['FILENAME'].mask = ma.nomask
            slot_datas[slot] = slot_data
        # make a shared time ndarray that is the union of the time sets in the
        # slots.  The ACA L0 telemetry has the same timestamps across slots,
        # so the only differences here are caused by different times in
        # non-TRAK across the slots (usually SRCH differences at the beginning
        # of the observation)
        shared_time = np.unique(np.concatenate([
            slot_datas[slot]['TIME'].data for slot in slots]))
        for msid in msids:
            hdr3_msid = confirm_msid(msid)
            slot = slot_for_msid(hdr3_msid)
            full_data = ma.zeros(len(shared_time),
                                 dtype=slot_datas[slot].dtype)
            full_data.mask = ma.masked
            fd_idx = search_both_sorted(shared_time,
                                        slot_datas[slot]['TIME'])
            full_data[fd_idx] = slot_datas[slot]
            # make a data dictionary to feed to the MSID constructor
            slot_data = {'vals': HDR3_DEF[hdr3_msid]['value'](full_data),
                         'desc': HDR3_DEF[hdr3_msid]['desc'],
                         'longdesc': HDR3_DEF[hdr3_msid]['longdesc'],
                         'times': shared_time,
                         'hdr3_msid': hdr3_msid}
            self[msid] = MSID(msid, start, stop, slot_data)
