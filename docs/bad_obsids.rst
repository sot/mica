.. _bad_obsids-label:

Mica.bad_obsids
===============

The :mod:`mica.bad_obsids` module includes a single function `bad_obsids()` which
returns a list of the obsids that are reasonable to exclude from most trending applications.
The list includes observations of Venus and observations with multiple obis.
The observations of Venus are problematic for star trending.  The observations with multiple
obis may have more than one attitude for the same obsid.
