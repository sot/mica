ACA dark current
-------------------------

The :mod:`mica.archive.aca_dark` package provides modules related to the ACA dark current:

 * `Dark current calibrations`_
 * `Dark current modeling`_
 * `Processing`_

Dark current calibrations
^^^^^^^^^^^^^^^^^^^^^^^^^^

The :mod:`mica.archive.aca_dark.dark_cal` module provides functions for retrieving
data for the ACA full-frame dark current calibrations which occur about four
times per year (see the `ACA dark calibrations TWiki page <http://occweb.cfa.harvard.edu/twiki/Aspect/AcaDarkCal>`_).

The functions available are documented in the :ref:`api_aca_dark` section, but the most useful are:

 * :func:`~mica.archive.aca_dark.dark_cal.dark_temp_scale`: get the temperature scaling correction
 * :func:`~mica.archive.aca_dark.dark_cal.get_dark_cal_dirs`: get an ordered dict of dark cal identifer and directory
 * :func:`~mica.archive.aca_dark.dark_cal.get_dark_cal_image`: get a single dark cal image
 * :func:`~mica.archive.aca_dark.dark_cal.get_dark_cal_props`: get properties (e.g. date, temperature) of a dark cal
 * :func:`~mica.archive.aca_dark.dark_cal.get_dark_cal_props_table`: get properties of dark cals over a time range as a table

As an example, let's plot the raw and corrected warm pixel fraction over the mission.  The correction in this case is done to a reference temperature of -15 C::

  from mica.archive.aca_dark import dark_cal
  from Ska.Matplotlib import plot_cxctime
  from Chandra.Time import DateTime

  dark_cals = dark_cal.get_dark_cal_dirs()
  times = []
  n100 = []
  n100_m15 = []
  npix = 1024. * 1024.

  for dark_id in dark_cals:
      print('Reading {}'.format(dark_id))
      props = dark_cal.get_dark_cal_props(dark_id, include_image=True)
      scale = dark_cal.dark_temp_scale(props['ccd_temp'], -15.0)
      image = props['image']
      times.append(props['date'])
      n100.append(np.count_nonzero(image > 100.0) / npix)
      n100_m15.append(np.count_nonzero(image * scale > 100.0) / npix)

  times = DateTime(times).secs
  figure(figsize=(6, 4))
  plot_cxctime(times, n100, 'o-', color='red')
  plot_cxctime(times, n100_m15, 's-', color='cyan')
  grid(True)
  xlim(DateTime('2000:001').plotdate, DateTime().plotdate)
  ylim(0, None)
  title('Warm pixel fraction')

.. image:: plots/aca_dark_warm_pix_frac.png
   :width: 500

Note that the temperature assigned to a dark calibration is the mean of the temperature
for the invidivual dark replicas (typically 5).  These in turn use ACA hdr3 diagnostic
telemetry for high-resolution temperature readouts which are available before and after
(but not during) each replica.

Processing
^^^^^^^^^^^^^^

This module updates the MICA ACA dark current archive when new dark current calibrations
are completed.  For details see the API documentation at :ref:`api_update_aca_dark`.
