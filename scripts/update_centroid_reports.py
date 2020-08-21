#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import argparse

from mica import centroid_dashboard


# Cheat. Needs entrypoint scripts
centroid_dashboard.update_observed_metrics()
