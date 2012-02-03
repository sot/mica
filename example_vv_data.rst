Example: V&V Report
--------------------

The V&V tool of the mica suite takes as inputs:

* an obspar
* a complete set of aspect level 1 files
* the log from the aspect 1 run for each interval

The tool then calculates the aspect solution residuals for each slot
for evaluation.  Results include:

* Annotated top level tables including a slots table and annotated
  guide and fid props tables.  
* The obspar
* A table of the time ranges and aspect intervals of the observation
* plots of residuals for each slot
* plots of the aspect solution
* the data used to make those plots
* annotated log file

Tests will include the current V&V tests, including:

* on a per-slot-basis

   * was this star BAD, UNUSED, or GOOD?
   * Does the RMS in Y,Z, or radius exceed low or high thresholds?
   * Does the Mean or Median offset in Y,Z, or radius exceed
     thresholds?
   * was the slot present in all intervals?
   * did the slot spend a significant fraction with a large deviation?
   * were there any log errors 

* on a per-obsid basis:

   * what should the V&V disposition of the observation be?


Example: Image Reconstruction Checking
--------------------------------------

It is expected that the tool will have the capability to assist in
evaluating image reconstruction by inclusion of a path by which the
aspect solution will be recalculated after removing a slot, and the
V&V outputs of the with and without-slot solutions will be directly
compared.  The V&V report pieces should therefore also have mechanisms
for a collective / comparison report for two or more runs of the same
obsid.

