ACA diagnostic telemetry
-------------------------

The :mod:`mica.archive.aca_hdr3` module works with Header 3 data
(extended ACA diagnostic telemetry) available in 8x8 ACA
L0 image data.  The module provies an MSID class and MSIDset class to fetch
these data as "pseudo-MSIDs" and return masked array data structures.
See `ACA HDR3 Pseudo-MSIDS`_ for the list of available pseudo-MSIDs.

   >>> from mica.archive import aca_hdr3
   >>> ccd_temp = aca_hdr3.MSID('ccd_temp', '2012:001', '2012:020')
   >>> type(ccd_temp.vals)
   'numpy.ma.core.MaskedArray'
   >>> from Ska.Matplotlib import plot_cxctime
   >>> figure(figsize=(5, 3.5))
   >>> plot_cxctime(ccd_temp.times, ccd_temp.vals, '.')

   .. image:: plots/plot_cxctime_ccd_temp.png

   >>> perigee_data = aca_hdr3.MSIDset(['ccd_temp', 'aca_temp', 'dac'],
   ...                                 '2012:125', '2012:155')
   >>> figure(figsize=(5, 3.5))
   >>> plot(perigee_data['aca_temp'].vals - perigee_data['ccd_temp'].vals,
   ...      perigee_data['dac'].vals, '.')
   >>> subplots_adjust(bottom=0.15)
   >>> ylabel('TEC DAC Control Level')
   >>> xlabel('ACA temp - CCD temp (C)')

   .. image:: plots/dac_vs_tempdiff.png

Retrieving pseudo-MSIDs with this module will be slower than
Ska.engarchive fetches of similar telemetry, as the aca_hdr3 module
reads from each of the collection of original fits.gz files for a specified
time range.  Ska.engarchive, in contrast, reads from HDF5 files (per
MSID) optimized for fast reads.::


   In [3]: %time ccd_temp = aca_hdr3.MSID('ccd_temp', '2012:001', '2012:020')
   CPU times: user 5.18 s, sys: 0.12 s, total: 5.29 s
   Wall time: 7.46 s

   In [9]: %time quick_ccd = Ska.engarchive.fetch.MSID('AACCCDPT', '2012:001', '2012:020')
   CPU times: user 0.02 s, sys: 0.00 s, total: 0.03 s
   Wall time: 0.81 s



ACA HDR3 Pseudo-MSIDS
---------------------

.. raw:: html
   :file: hdr3_only_msids.html
