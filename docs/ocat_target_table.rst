Chandra Data Archive and Ocat
=============================

The :mod:`mica.archive.cda.services` module provides a Python interface to the
Chandra Data Archive (CDA) web services and an interface to a local disk copy of
the Observation Catalog (Ocat) if that is available.

The CDA web services interface provides access to the following. Note however
that the CDA is not accessible from the GRETA network.

- Ocat details: full details from the Ocat (124 fields) for the mission, plus
  access to ACIS windows, roll requirements and time requirements for
  applicable observations.
- Ocat summary: summary data from the Ocat (26 fields).
- Proposal abstracts: abstract information for each observation.
- Archive file list: list of raw data files in the Chandra archive for each
  observation.

On the HEAD and GRETA networks, a copy of the Ocat details table (124 fields) is
maintained and updated daily. This is a compressed HDF5 file which is
approximately 5 Mb and can be synced to a local (laptop) computer.
It is located at::

  ${SKA}/data/mica/archive/ocat_target_table.h5

This local file version can be used for faster, network-free queries of the
Ocat.

This is the type of content that can be seen directly at:

https://cda.harvard.edu/srservices/ocatDetails.do?obsid=2121

Local Ocat access
-----------------

The :func:`~mica.archive.cda.services.get_ocat_local` function provides fast
access to the Ocat if the local Ocat HDF5 file is available on disk. In
particular reading the entire Ocat or querying on a target name substring is at
least 10 times faster than using the `Web access to CDA`_.

The local version of the Ocat is cached in memory the first time it is read in
full, so this means subsequent full reads or searches for target names will be
much faster. The cache expires in one hour in case you have a long-running job.

Any of the fields in the returned table can be used for filtering the query
via exact matches of the parameters::

   >>> from mica.archive.cda import get_ocat_local
   >>> dat = get_ocat_local(obsid=5438)
   >>> dat['target_name']
   'R Aqr'

   >>> dat = get_ocat_local()  # Get the whole Ocat for the mission
   >>> dat = get_ocat_local(instr='HRC-S', grat='LETG')  # All HRC-S LETG obs
   >>> dat
   <Table length=734>
   seq_num  status  obsid  pr_num  target_name ... evfil evfil_lo evfil_ra efficient spwin
                                               ...         keV      keV
     str6   str11   int64   str8      str29    ...  str1 float64  float64     str1    str1
   ------- -------- ----- -------- ----------- ... ----- -------- -------- --------- -----
           archived 62635                      ...     0       --       --         0     0
           archived 62636                      ...     0       --       --         0     0
       ...      ...   ...      ...         ... ...   ...      ...      ...       ...   ...
    901342 archived 20944 18900560     Mkn 421 ...     0       --       --         0     0
    901410 archived 21389 20900088 PKS0405-123 ...     0       --       --         0     0
    901410 archived 21955 20900088 PKS0405-123 ...     0       --       --         0     0

The ``target_name`` is handled specially. By default, this parameter matches
any substring in the Ocat target name field, ignoring spaces::

   # Targets with 'jet' in name.
   >>> dat = get_ocat_local(target_name='jet')

However, if the target name is a valid catalog source name that can be resolved
by the CDS name resolver, then the target name can be used for a radial cone
search around the target position by specifying ``resolve_name=True``::

   # Observations within 4 arcmin of 3C273
   >>> dat = get_ocat_local(target_name='3c273', resolve_name=True, radius=4)

   # Observations within 1 degree of RA, Dec = 10, 10
   >>> dat = get_ocat_local(ra=10, dec=10, radius=60)

One convenience in this function is that when the ``obsid`` is specified, then
by default the result is returned as a ``dict`` instead of a ``Table``. This
saves the boilerplate step of selecting element-0 of a Table.

Advanced users may consider use of the ``where`` keyword argument, which allows
fast in-kernel filtering of the table contents.
See https://www.pytables.org/usersguide/condition_syntax.html for details.
This feature is used internally for the exact parameter matching as well as
for the positional cone search around a coordinate.

Web access to CDA
-----------------

The Chandra Data Archive hosts a web service that provides access to various
elements of the archive. Examples follow.

Ocat
^^^^
The :func:`~mica.archive.cda.services.get_ocat_web` function for web access to
the Ocat is similar to the local file access documented above, and all of those
examples will work just by changing ``get_ocat_local`` to ``get_ocat_web``.

Additional functionality is described in the function docs
:func:`~mica.archive.cda.services.get_ocat_web`. In particular the Web API parameters
can be used via the function call.

Proposal abstracts
^^^^^^^^^^^^^^^^^^
The :func:`~mica.archive.cda.services.get_proposal_abstract` function allows
getting information about the proposal and the abstract::

    >>> from mica.archive.cda import get_proposal_abstract
    >>> get_proposal_abstract(obsid=8000)
    {'abstract': 'We propose the Chandra-COSMOS survey which will provide an '
                 'unprecedented combination of contiguous area, depth and '
                 'resolution. 36 densely tiled observations will cover the central '
                 '0.7 sq.deg. COSMOS field to a uniform 200ksec depth. COSMOS '
                 'explores the coupled evolution of galaxies, dark matter halos '
                 'and AGNs (massive black holes) largely free of cosmic variance. '
                 'COSMOS is a comprehensive survey including: HST, Spitzer, '
                 'Subaru, VLT, Magellan, VLA, MAMBO, GALEX, & potentially EVLA & '
                 'ALMA. Chandra resolution & sensitivity enables the study of '
                 'large scale phenomena: (1) influence of the surrounding '
                 'environment; (2) interaction between galaxies; (3) influence of '
                 'groups and clusters',
     'principal_investigator': 'Martin Elvis',
     'proposal_number': '08900073',
     'proposal_title': 'THE CHANDRA-COSMOS SURVEY'}

Archive file list
^^^^^^^^^^^^^^^^^
The :func:`~mica.archive.cda.services.get_archive_file_list` function allows
getting a list of archive files for given ``obsid``, ``detector``, ``level``,
and ``dataset`` (and possibly other parameters).

Example::

    >>> get_archive_file_list(obsid=2365, detector='pcad',
    ...                           subdetector='aca', level=1, obi=2)
    <Table length=27>
            Filename            Filesize      Timestamp
                str30               int64          str19
    ------------------------------ -------- -------------------
    pcadf126690624N007_asol1.fits  7300800 2021-04-09 08:04:29
    pcadf02365_002N001_asol1.fits  4728960 2021-04-09 08:04:30
                            ...      ...                 ...
    pcadf126695890N007_adat61.fits  1293120 2021-04-09 08:04:28
    pcadf126695890N007_adat71.fits  1293120 2021-04-09 08:04:28

    >>> get_archive_file_list(obsid=400, detector='acis', level=2, filetype='evt2')
    <Table length=1>
            Filename         Filesize      Timestamp
            str24            int64          str19
    ------------------------ -------- -------------------
    acisf00400N007_evt2.fits  4619520 2011-07-08 13:52:57


Ocat table fields
-----------------

Some of these fields may be described in the chaser help at:

https://cda.harvard.edu/chaser/chaserFieldHelp.html


============================= ================================================================
 Column                       Description
============================= ================================================================
seq_num                       Sequence Number
status                        Status (unobserved, archived, untriggered, etc)
obsid                         Obsid
pr_num                        Proposal Number
target_name                   Target name
grid_name                     database id of grid name if grid observation
instr                         Instrument
grat                          Grating (HETG, LETG, NONE)
type                          Observation type (TOO, GO, GTO)
obs_cycle                     Observation cycle
prop_cycle                    Proposal cycle
charge_cycle                  Charge cycle
start_date                    Observation start date
public_avail                  Date publicly available
readout_detector              Detector (which HRC detector or string of actual ACIS ccds)
datamode                      Instrument data mode
joint                         Joint observatories (string)
hst                           HST time (orbits?)
noao                          NOAO time
nrao                          NRAO time
rxte                          RXTE time
spitzer                       SPITZER time
suzaku                        SUZAKU time
xmm                           XMM time
swift                         SWIFT time
nustar                        NUSTAR time
category                      Science category
seg_max_num
prop_title                    Proposal Title
pi_name                       Principal investigator name
observer                      Observer name
app_exp                       Approved exposure time (ks)
exp_time                      Actual exposure time (ks)
ra                            Target Right Ascension
dec                           Target Declination
soe_roll                      Roll from SOE
time_crit                     Time critical flag
y_off                         Y offset
z_off                         Z offset
x_sim                         X SIM
z_sim                         Z SIM
raster                        Raster flag
obj_type                      Object type
obj                           Solar system object name
photo                         Photometry flag
vmag                          V Mag of object
est_cnt_rate                  Estimated count rate
forder_cnt_rate               First order count rate
count_rate
event_count
dither                        Dither flag
y_amp                         Dither Y amplitude if custom dither
y_freq                        Dither Y frequency if custom dither
y_phase                       Dither Y phase if custom dither
z_amp                         Dither Z amplitude if custom dither
z_freq                        Dither Z frequency if custom dither
z_phase                       Dither Z phase if custom dither
roll                          Roll constraint flag
window                        Window constraint flag
unint                         Uninterrupt constraint flag
pointing_update               Pointing update constraint flag
monitor                       Monitor series flag
pre_id                        Obsid of previous observation in monitor series
mon_min                       Min days from pre_id for monitor observation
mon_max                       Max days from pre_id for monitor observation
group_id                      Database group id
constr
epoch
period
pstart
ps_marg
pend
pe_marg
multitel
multitel_obs
multitel_int
constr_rmk                    Constraint in remarks flag
too_type
too_start
too_stop
alt_group
alt_trig
simode                        Science Instrument (SI) mode
hrc
spect_mode
blank_en
u_hi
v_hi
u_lo
v_lo
timing
z_blk
acis
mode                          ACIS mode (CC or TE)
bep_pack                      ACIS BEP PACK (F, G, VF, F+B)
dropped_chip_cnt              Dropped chip count
i0                            ACIS I0 ccd status (Y, N, optional with number, or D if dropped)
i1                            ACIS I1 ccd status (Y, N, optional with number, or D if dropped)
i2                            ACIS I2 ccd status (Y, N, optional with number, or D if dropped)
i3                            ACIS I3 ccd status (Y, N, optional with number, or D if dropped)
s0                            ACIS S0 ccd status (Y, N, optional with number, or D if dropped)
s1                            ACIS S1 ccd status (Y, N, optional with number, or D if dropped)
s2                            ACIS S2 ccd status (Y, N, optional with number, or D if dropped)
s3                            ACIS S3 ccd status (Y, N, optional with number, or D if dropped)
s4                            ACIS S4 ccd status (Y, N, optional with number, or D if dropped)
s5                            ACIS S5 ccd status (Y, N, optional with number, or D if dropped)
spectra_max_count             Spectra Max Count
multiple_spectral_lines       Multiple spectral lines expected (Y, N)
subary                        ACIS subarray (CUSTOM, NONE)
strt_row                      Start row of ACIS subarray
row_cnt                       Number of rows of ACIS subarray
d_cyc
sec_cnt
pr_time
sec_time
f_time
oc_sum
oc_row
oc_col
evfil
evfil_lo
evfil_ra
efficient                     ACIS use most efficient (Y, N)
spwin                         Spatial window (Y, N)
============================= ================================================================
