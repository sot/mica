#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import mica.centroid_dashboard

# Cheat. Needs entrypoint scripts
mica.centroid_dashboard.update_observed_metrics(save=True, make_plots=True)
