#!/usr/bin/env python

# The cco2 44Ti Magkotsios plot has a high-abundance yellow region
# Save a list of particle IDs for everything within that yellow region
# Save a list of all the cco2 SDF files in order with tpos values
# These lists will be read by fetch_data.c when it does its job next

# Last modified 9/21/18 by Greg Vance

import numpy as np
import os

# Names of output list files
PID_FILE = "pid_list"
SDF_FILE = "sdf_list"

# Data file locations on saguaro for cco2
PLOT_DIR = "/home/gsvance/results/plotting/"
PLOT_FILE = PLOT_DIR + "cco2_plotting.out"
COLS_FILE = PLOT_DIR + "columns"
SDF_DIR = "/home/gsvance/supernova_data/cco2/cco2/"

# Boundaries (by eye) of the high-abundance "yellow region"
# See young/shared_files/magkotsios/cco2_magkotsios_hist.png
PEAK_TEMP_MIN = 8e9 # K
PEAK_TEMP_MAX = 10.5e9 # K
PEAK_RHO_MIN = 10.**6.5 # g/cm^3
PEAK_RHO_MAX = 10.**7. # g/cm^3

# Read columns file and extract indicies of relevant columns
with open(COLS_FILE, "r") as cols_file:
	cols = cols_file.read().strip().split('\n')
i_pid = cols.index("id")
i_peak_temp = cols.index("peak temp")
i_peak_rho = cols.index("peak density")

# Read the needed columns of data from the cco2 plotting file
datatype = {"names": ("pid", "peak temp", "peak rho"), "formats": (np.int,
	np.float, np.float)}
pid, peak_temp, peak_rho = np.loadtxt(PLOT_FILE, datatype, delimiter=", ",
	skiprows=1, usecols=(i_pid, i_peak_temp, i_peak_rho), unpack=True)
print "cco2 plotting data obtained"

# Make a list of the particle ID for every particle in the yellow region
yellow_pids = []
for i in xrange(pid.size):
	if PEAK_TEMP_MIN < peak_temp[i] < PEAK_TEMP_MAX \
		and PEAK_RHO_MIN < peak_rho[i] < PEAK_RHO_MAX:
		yellow_pids.append(pid[i])
yellow_pids.sort()

# Write out the list of particle IDs from the yellow region
with open(PID_FILE, "w") as pid_file:
	head = "n_id %d\n" % (len(yellow_pids))
	pid_file.write(head)
	for ypid in yellow_pids:
		line = "%d\n" % (ypid)
		pid_file.write(line)
print "list of yellow region ids produced"

# Function to extract tpos values from SDF files
def get_tpos(sdf_name):
	with open(sdf_name, "rb") as sdf_file:
		for line in sdf_file:
			if line.startswith("float tpos = "):
				return line.strip("float ps=;\n")
	msg = "SDF %s missing tpos value" % (sdf_name)
	raise ValueError(msg)

# Form a list of all the cco2 SDF file names
sdf_names = []
for maybe_sdf in os.listdir(SDF_DIR):
	extension = maybe_sdf.split('.')[-1]
	try:
		iter = int(extension)
	except ValueError:
		continue
	sdf_name = SDF_DIR + maybe_sdf
	tpos = get_tpos(sdf_name)
	sdf_names.append((iter, tpos, sdf_name))
sdf_names.sort(key=lambda s : s[0])

# Write out the list of SDF files with their data
with open(SDF_FILE, "w") as sdf_file:
	head = "n_sdf %d\n" % (len(sdf_names))
	head += "iter tpos name\n"
	sdf_file.write(head)
	for iter, tpos, name in sdf_names:
		line = "%d %s %s\n" % (iter, tpos, name)
		sdf_file.write(line)
print "list of cco2 sdf files produced"

