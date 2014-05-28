"""Define bins for dark current distribution calculations"""
import numpy as np

x0 = 2.
x1 = 20000.
dx = 0.05
log_bins = np.arange(np.log10(x0), np.log10(x1), dx)
bins = 10 ** log_bins
bin_centers = 10 ** ((log_bins[1:] + log_bins[:-1]) / 2.0)
