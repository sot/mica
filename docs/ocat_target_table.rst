OCAT Target Table
---------------------

The :mod:`mica.archive.ocat_target_table` module includes code
to fetch a mica version of the content from the OCAT details 
aka Chaser aka target table.

This is the type of content that can be seen directly at:

https://cda.harvard.edu/srservices/ocatDetails.do?obsid=2121

Using this module to get the data for all science observations:

   >>> from mica.archive.ocat_target_table import get_ocat_target_table
   >>> dat = get_ocat_target_table()
   >>> dat[dat['OBSID'] == 5438][0]['TARGET_NAME']
   'R Aqr'

Some of these fields may be described in the chaser help at:

https://cda.harvard.edu/chaser/chaserFieldHelp.html


Data table fields
^^^^^^^^^^^^^^^^^

============================= ================================================================
 Column                       Description
============================= ================================================================
SEQ_NUM                       Sequence Number
STATUS                        Status (unobserved, archived, untriggered, etc)
OBSID                         Obsid
PR_NUM                        Proposal Number
TARGET_NAME                   Target name
GRID_NAME                     database id of grid name if grid observation
INSTR                         Instrument
GRAT                          Grating (HETG, LETG, NONE)
TYPE                          Observation type (TOO, GO, GTO)
OBS_CYCLE                     Observation cycle
PROP_CYCLE                    Proposal cycle
CHARGE_CYCLE                  Charge cycle
START_DATE                    Observation start date
PUBLIC_AVAIL                  Date publicly available
READOUT_DETECTOR              Detector (which HRC detector or string of actual ACIS ccds)
DATAMODE                      Instrument data mode
JOINT                         Joint observatories (string)
HST                           HST time (orbits?)
NOAO                          NOAO time
NRAO                          NRAO time
RXTE                          RXTE time
SPITZER                       SPITZER time
SUZAKU                        SUZAKU time
XMM                           XMM time
SWIFT                         SWIFT time
NUSTAR                        NUSTAR time
CATEGORY                      Science category
SEG_MAX_NUM                   
PROP_TITLE                    Proposal Title
PI_NAME                       Principal investigator name
OBSERVER                      Observer name
APP_EXP                       Approved exposure time (ks)
EXP_TIME                      Actual exposure time (ks)
RA                            Right Ascension
Dec                           Declination
SOE_ROLL                      Roll from SOE
TIME_CRIT                     Time critical flag
Y_OFF                         Y offset 
Z_OFF                         Z offset
X_SIM                         X SIM
Z_SIM                         Z SIM
RASTER                        Raster flag
OBJ_TYPE                      Object type
OBJ                           Solar system object name
PHOTO                         Photometry flag
VMAG                          V Mag of object
EST_CNT_RATE                  Estimated count rate
Forder_CNT_RATE               First order count rage
COUNT_RATE                    
EVENT_COUNT
DITHER                        Dither flag
Y_AMP                         Dither Y amplitude if custom dither
Y_FREQ                        Dither Y frequency if custom dither
Y_PHASE                       Dither Y phase if custom dither
Z_AMP                         Dither Z amplitude if custom dither
Z_FREQ                        Dither Z frequency if custom dither
Z_PHASE                       Dither Z phase if custom dither
ROLL                          Roll constraint flag
WINDOW                        Window constraint flag
UNINT                         Uninterrupt constraint flag
POINTING_UPDATE               Pointing update constraint flag
MONITOR                       Monitor series flag
PRE_ID                        Obsid of previous observation in monitor series
MON_MIN                       Min days from pre_id for monitor observation
MON_MAX                       Max days from pre_id for monitor observation
GROUP_ID                      Database group id
CONSTR
EPOCH
PERIOD
PSTART
PS_MARG
PEND
PE_MARG
MULTITEL
MULTITEL_OBS
MULTITEL_INT
CONSTR_RMK                    Constraint in remarks flag
TOO_TYPE
TOO_START
TOO_STOP
ALT_GROUP
ALT_TRIG
SIMODE                        Science Instrument (SI) mode
HRC
SPECT_MODE
BLANK_EN
U_HI
V_HI
U_LO
V_LO
TIMING
Z_BLK                  
ACIS  
MODE                          ACIS mode (CC or TE)
BEP_PACK                      ACIS BEP PACK (F, G, VF, F+B)
DROPPED_CHIP_CNT              Dropped chip count
I0                            ACIS I0 ccd status (Y, N, optional with number, or D if dropped)
I1                            ACIS I1 ccd status (Y, N, optional with number, or D if dropped)
I2                            ACIS I2 ccd status (Y, N, optional with number, or D if dropped)
I3                            ACIS I3 ccd status (Y, N, optional with number, or D if dropped)
S0                            ACIS S0 ccd status (Y, N, optional with number, or D if dropped)
S1                            ACIS S1 ccd status (Y, N, optional with number, or D if dropped)
S2                            ACIS S2 ccd status (Y, N, optional with number, or D if dropped)
S3                            ACIS S3 ccd status (Y, N, optional with number, or D if dropped)
S4                            ACIS S4 ccd status (Y, N, optional with number, or D if dropped)
S5                            ACIS S5 ccd status (Y, N, optional with number, or D if dropped)
SPECTRA_MAX_COUNT             Spectra Max Count
MULTIPLE_SPECTRAL_LINES       Multiple spectral lines expected (Y, N)
SUBARY                        ACIS subarray (CUSTOM, NONE)
STRT_ROW                      Start row of ACIS subarray
ROW_CNT                       Number of rows of ACIS subarray
D_CYC
SEC_CNT
PR_TIME
SEC_TIME
F_TIME
OC_SUM
OC_ROW
OC_COL
EVFIL
EVFIL_LO
EVFIL_RA
EFFICIENT                     ACIS use most efficient (Y, N)
SPWIN                         Spatial window (Y, N)
============================= ================================================================

The HDF5 in-kernel searches may be faster working with the table directly for some
operations.
