
"""
Experimental code to work with ACA L0 Header 3 data
"""
import re
import numpy as np
import collections
from scipy.interpolate import interp1d
from Chandra.Time import DateTime
import mica.archive.aca_l0

class MissingDataError(Exception):
    pass

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
aca_dtype = [('TIME', '>f8'), ('QUALITY', '>i4'), 
             ('HD3TLM62', '|u1'),
             ('HD3TLM63', '|u1'), ('HD3TLM64', '|u1'), ('HD3TLM65', '|u1'),
             ('HD3TLM66', '|u1'), ('HD3TLM67', '|u1'), ('HD3TLM72', '|u1'),
             ('HD3TLM73', '|u1'), ('HD3TLM74', '|u1'), ('HD3TLM75', '|u1'),
             ('HD3TLM76', '|u1'), ('HD3TLM77', '|u1'),
             ]

aca_dtype_names = [k[0] for k in aca_dtype]


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
#adpoly = np.polyfit(x, a_to_d['tempC'], 6)
ad_func = interp1d(x, a_to_d['tempC'], kind='cubic', bounds_error=False)


def ad_temp(msids):
    def func(slot_data):
        sum = two_byte_sum(msids)(slot_data)
        return ad_func(sum)
    return func
    

hdr3_def = \
{'062': {'desc': 'AC status word',
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
               'value': ad_temp(['HD3TLM66', 'HD3TLM67']),
               'longdesc': """
A/D converter reading for the CCD moly base thermistor number 1
"""},
       '072': {'desc': 'A/D CCD molyb therm 2',
               'value': ad_temp(['HD3TLM72', 'HD3TLM73']),
               'longdesc': """
A/D converter reading for the CCD moly base thermistor number 2
"""},
       '074': {'desc': 'A/D CCD detector therm',
               'value': ad_temp(['HD3TLM74', 'HD3TLM75']),
               'longdesc': """
A/D converter reading for the CCD detector thermistor
"""},
       '076': {'desc': 'A/D +5 volt PS',
               'value': two_byte_sum(['HD3TLM76', 'HD3TLM77'], scale=0.20518),
               'longdesc': """
A/D converter reading for the +5 volt power supply; 1 LSB=0.30518 mv
"""},
       '162': {'desc': 'A/D +15 volt PS',
               'value': two_byte_sum(['HD3TLM62', 'HD3TLM63'], scale=0.61035),
               'longdesc': """
A/D converter reading for the +15 volt power supply; 1 LSB=0.61035 mv
"""},
       '164': {'desc': 'A/D -15 volt PS',
               'value': two_byte_sum(['HD3TLM64', 'HD3TLM65'], scale=0.61035),
               'longdesc': """
A/D converter reading for the -15 volt power supply; 1 LSB=0.61035 mv
"""},
       '166': {'desc': 'A/D +27 volt PS',
               'value': two_byte_sum(['HD3TLM66', 'HD3TLM67'], scale=1.04597),
               'longdesc': """
A/D converter reading for the +27 volt power supply; 1 LSB=1.04597 mv
"""},
       '172': {'desc': 'A/D analog ground',
               'value': ad_temp(['HD3TLM72', 'HD3TLM73']),
               'longdesc': """
A/D converter reading for analog ground; 1 LSB=0.30518 mv
"""},
       '174': {'desc': 'A/D for A/D convertor therm',
               'value': ad_temp(['HD3TLM74', 'HD3TLM75']),
               'longdesc': """
A/D converter reading for the A/D converter thermistor.
"""},
       '176': {'desc': 'A/D secondary mirror therm. HRMA side',
               'value': ad_temp(['HD3TLM76', 'HD3TLM77']),
               'longdesc': """
A/D converter reading for the secondary mirror thermistor, HRMA side
"""},
       '262': {'desc': 'A/D secondary mirror therm. Opp HRMA side',
               'value': ad_temp(['HD3TLM62', 'HD3TLM63']),
               'longdesc': """
A/D converter reading for the secondary mirror thermistor, Opposite from the
HRMA side
"""},
       '264': {'desc': 'A/D primary mirror therm. HRMA side',
               'value': ad_temp(['HD3TLM64', 'HD3TLM65']),
               'longdesc': """
A/D converter reading for the primary mirror thermistor, HRMA side
"""},
       '266': {'desc': 'A/D primary mirror therm. Opp HRMA side',
               'value': ad_temp(['HD3TLM66', 'HD3TLM67']),
               'longdesc': """
A/D converter reading for the primary mirror thermistor, opposite from the
HRMA side
"""},
       '272': {'desc': 'A/D AC housing therm.  HRMA side',
               'value': ad_temp(['HD3TLM72', 'HD3TLM73']),
               'longdesc': """
A/D converter reading for the AC housing thermistor, HRMA side
"""},
       '274': {'desc': 'A/D AC housing therm.  Opp HRMA side',
               'value': ad_temp(['HD3TLM74', 'HD3TLM75']),
               'longdesc': """
A/D converter reading for the AC housing thermistor, opposite HRMA side
"""},
       '276': {'desc': 'A/D lens cell therm.',
               'value': ad_temp(['HD3TLM76', 'HD3TLM77']),
               'longdesc': """
A/D converter reading for the lens cell thermistor
"""},
       '362': {'desc': 'Processor stack pointer and telem update counter',
               'longdesc': """
A word containing the processor data stack pointer in the high byte, and
an update counter in the low byte that increments once for every 1.025
second telemetry update.
"""},
       '364': {'desc': 'Science header pulse period',
               'longdesc': """
The science header pulse period, as measured by the PEA; 1 LSB = 2 microseconds
""",
               'nbytes': 4},
       '372': {'desc': '16-bit zero offset for pixels from CCD quad A',
               'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant A; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},
       '374': {'desc': '16-bit zero offset for pixels from CCD quad B',
               'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant B; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},
       '376': {'desc': '16-bit zero offset for pixels from CCD quad C',
               'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant C; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},
       '462': {'desc': '16-bit zero offset for pixels from CCD quad D',
               'longdesc': """
A 16-bit zero offset for pixels read from CCD quadrant D; 1 LSB = 1 A/D
converter count (nominally 5 electrons)
"""},
       '464': {'desc': '32-bit zero offset for pixels from CCD quad A',
               'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant A; 1 LSB = 2^-16
A/D converter counts
""",
               'nbytes': 4},
       '472': {'desc': '32-bit zero offset for pixels from CCD quad B',
               'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant B; 1 LSB = 2^-16
A/D converter counts
""",
               'nbytes': 4},
       '476': {'desc': '32-bit zero offset for pixels from CCD quad C',
               'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant C; 1 LSB = 2^-16
A/D converter counts
""",
               'nbytes': 4},
       '564': {'desc': '32-bit zero offset for pixels from CCD quad D',
               'longdesc': """
A 32-bit zero offset for pixels read from CCD quadrant D; 1 LSB = 2^-16
A/D converter counts
""",
               'nbytes': 4},
       '572': {'desc': 'last CCD flush duration',
               'longdesc': """
The time required for the most recent flush of the CCD; 1 LSB=2 microseconds
"""},
       '574': {'desc': 'CCD row shift clock period',
               'longdesc': """
The CCD row shift clock period currently in effect; 1 LSB = 1 microsecond
"""},
       '576': {'desc': 'average backround reading',
               'longdesc': """
An overall average background reading derived from the most recent CCD readout.
This is an average from all tracked images and from all search readout
segments. One LSB = 1 A/D converter count (nominally 5 electrons).
""",
               'value': two_byte_sum(['HD3TLM76', 'HD3TLM77'], scale=5)},
       '662': {'desc': 'Header 1 for DSN records',
               'longdesc': """
Header 1 for Deep Space Network record.
""",
               'nbytes': 6},
       '672': {'desc': 'record counter and Header 2 for DSN',
               'longdesc': """
The record counter and Header 2 for Deep Space Network records.
The record counter occupies the three high order bytes, and Header 2
occupies the low order byte.
""",
               'nbytes': 4},
       '676': {'desc': 'CCD temperature',
               'longdesc': """
CCD temperature.  1 LSB=0.01/(2^16) degrees C.  The high order 16 bits give the
CCD temperature in units of 1 LSB = 0.01 degrees C.
""",
               'nbytes': 4,
               'value': two_byte_sum(['HD3TLM76', 'HD3TLM77'], scale=0.01)},
       '764': {'desc': 'CCD setpoint',
               'longdesc': """
The CCD temperature control setpoint; 1 LSB=0.01 degrees C
""",
               'value': two_byte_sum(['HD3TLM64', 'HD3TLM65'])},
       '766': {'desc': 'temperature for position/angle cal',
               'longdesc': """
The temperature used in the angle calibration equations that convert star
positions from CCD row and column coordinates to Y and Z angles for OBC
telemetry; 1 LSB = 1/256 degrees C.
""",
               'nbytes': 4,
               'value': two_byte_sum(['HD3TLM72', 'HD3TLM73'],
                                     scale=1 / 256.)},
       '774': {'desc': 'RAM address of last write-read test failure',
               'longdesc': """
The address in RAM of the failure most recently detected by the RAM
write-and-read test
"""},
       '776': {'desc': 'TEC DAC number',
               'value': two_byte_sum(['HD3TLM76', 'HD3TLM77']),
               'longdesc': """
The number most recently written the the [sic] TEC power control DAC.
"""}}

msid_aliases = {'dac': {'hdr3': '776'},
                'aca_temp': {'hdr3': '766'},
                'ccd_temp': {'hdr3': '676'}}


default_nbytes = 2



class MSID(object):
    def __init__(self, msid, start, stop, msid_data=None):
        self.tstart = DateTime(start).secs
        self.tstop = DateTime(stop).secs
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        if msid_data is None:
            msid_data = MSIDset([msid], start, stop)[msid]
        self.hdr3_msid = msid_data['hdr3_msid']
        self.msid = msid
        self.vals = msid_data['vals']
        self.times = msid_data['times']


def confirm_msid(req_msid):
    if req_msid in msid_aliases:
        return msid_aliases[req_msid]['hdr3']
    else:
        if req_msid not in hdr3_def:
            raise MissingDataError("msid %s not found" % req_msid)
        else:
            return req_msid

def slot_for_msid(msid):
    mmatch = re.match('(\d)\d\d', msid)
    slot = int(mmatch.group(1))
    return slot


class MSIDset(collections.OrderedDict):
    def __init__(self, msids, start, stop):
        super(MSIDset, self).__init__()
        self.tstart = DateTime(start).secs
        self.tstop = DateTime(stop).secs
        self.datestart = DateTime(self.tstart).date
        self.datestop = DateTime(self.tstop).date
        self.slot_data = {}
        slots = [slot_for_msid(confirm_msid(msid)) for msid in msids]
        for slot in slots:
            slot_data = mica.archive.aca_l0.get_slot_data(
                self.tstart, self.tstop, slot,
                imgsize=[8], columns=aca_dtype_names)
            self.slot_data[slot] = slot_data
        for msid in msids:
            hdr3_msid = confirm_msid(msid)
            slot = slot_for_msid(hdr3_msid)
            slot_data = {'vals': hdr3_def[hdr3_msid]['value'](self.slot_data[slot]),
                         'times': self.slot_data[slot]['TIME'],
                         'hdr3_msid': hdr3_msid}
            self[msid] = MSID(msid, start, stop, slot_data)
