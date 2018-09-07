#!/usr/bin/env python

# The cco2 44Ti Magkotsios plot has a high-abundance yellow region
# Report the particle IDs of everything inside that yellow region
# ...

# Last modified 9/7/18 by Greg Vance

import numpy as np

# Data file locations
DATADIR = "/home/gsvance/results/plotting/"
DATAFILE = DATADIR + "cco2_plotting.out"
COLSFILE = DATADIR + "columns"

# Boundaries (by eye) of the high-abundance "yellow region"
# See young/shared_files/magkotsios/cco2_magkotsios_hist.png
PEAKTEMPMIN = 8e9 # K
PEAKTEMPMAX = 10.5e9 # K
PEAKRHOMIN = 10.**6.5 # g cm^-3
PEAKRHOMAX = 10.**7.0 # g cm^-3

# Read columns file and extract indicies of relevant columns
with open(COLSFILE, "r") as colsfile:
	cols = colsfile.read().strip().split('\n')
i_pid = cols.index("id")
i_peaktemp = cols.index("peak temp")
i_peakrho = cols.index("peak density")

# Read the needed columns of the cco2 data file
datatype = {"names": ("pid", "peaktemp", "peakrho"), "formats": (np.int,
	np.float, np.float)}
pid, peaktemp, peakrho = np.loadtxt(DATAFILE, datatype, delimiter=", ",
	skiprows=1, usecols=(i_pid, i_peaktemp, i_peakrho), unpack=True)
print "cco2 data loaded from file"

# Make a list of the particle ID for every particle in the yellow region
yellow = []
for i in xrange(pid.size):
	if PEAKTEMPMIN < peaktemp[i] < PEAKTEMPMAX \
		and PEAKRHOMIN < peakrho[i] < PEAKRHOMAX:
		yellow.append(pid[i])
yellow = np.array(yellow)

# Print details
print yellow
print yellow.size

