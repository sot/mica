Example: Temperatures Report
----------------------------

The temperatures topic report area would retrieve temperature
information from the Ska engineering archive and from either raw ACA
hdr 3 telemetry or files from the perigee_health processing.  The
resulting datas would include:

* DAC, CCD temp (fine), CCD temp, AACCCDPT (rough), and ACA temp statistics tables
* plots for each plus the DAC vs delta temp plot from perigee_health
* a test result plot for exceeded thresholds

Example: Collective Temperature & Star Tracking Report
------------------------------------

The test infrastructure we intend should allow the integration of data
types to ask new questions.  We are mostly interested in temperature
as it applies to camera capabilities, to answer questions such as:

* How much more likely is loss of track on a 10.3mag star with the
  current warm pixel distribution, if the CCD is at -17C?

Should this framework answer that question?
 
