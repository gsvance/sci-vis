#!/usr/bin/env python

# The cco2 histogram of Fe/Ti vs radiation entropy displays 3 distinct regions:
#  - West region: bright blob with a few red bins at log(S_rad) ~ 22.5
#  - East region: extended shape with many bright bins at log(S_rad) ~ 22.6
#  - North region: small bright patch of bins up at log(Fe/Ti) ~ 3.5
# Save a list of particle IDs for everything within these 3 regions
# Save a list of all the cco2 SDF files in order with tpos values
# These lists will be read by fetch_data.c when it does its job next
# Also write a file listing particle IDs by region for the plots later

# Last modified 9/28/18 by Greg Vance

import numpy as np
import os

# Names of output list files
PID_FILE = "pid_list"
SDF_FILE = "sdf_list"

# Name of output file listing particle IDs by region
REGIONS_FILE = "regions.out"

# Data file locations on saguaro for cco2
PLOT_DIR = "/home/gsvance/results/plotting/"
PLOT_FILE = PLOT_DIR + "cco2_plotting.out"
COLS_FILE = PLOT_DIR + "columns"
SDF_DIR = "/home/gsvance/supernova_data/cco2/cco2/"

# Approximate boundaries of the regions from an interactive MPL session
# See young/shared_files/casA_plots/hist_fe_ti_vs_rad_ent.png
WEST_SR_MIN = 10.**22.47 # K^3 cm^3 g^-1
WEST_SR_MAX = 10.**22.56 # K^3 cm^3 g^-1
WEST_FETI_MIN = 10.**2.4
WEST_FETI_MAX = 10.**2.8
EAST_SR_MIN = 10.**22.56 # K^3 cm^3 g^-1
EAST_SR_MAX = 10.**22.69 # K^3 cm^3 g^-1
EAST_FETI_MIN = 10.**2.2
EAST_FETI_MAX = 10.**2.8
NORTH_SR_MIN = 10.**22.54 # K^3 cm^3 g^-1
NORTH_SR_MAX = 10.**22.62 # K^3 cm^3 g^-1
NORTH_FETI_MIN = 10.**3.5
NORTH_FETI_MAX = 10.**3.8

# Read columns file and extract indicies of relevant columns
with open(COLS_FILE, "r") as cols_file:
	cols = cols_file.read().strip().split('\n')
i_pid = cols.index("id")
i_peak_temp = cols.index("peak temp")
i_peak_rho = cols.index("peak density")
i_fe = cols.index("X_{Fe}")
i_56ni = cols.index("X_{56Ni}")
i_44ti = cols.index("X_{44Ti}")

# Read the needed columns of data from the cco2 plotting file
datatype = {"names": ("pid", "peak temp", "peak rho", "Fe", "56Ni", "44Ti"),
	"formats": (np.int, np.float, np.float, np.float, np.float, np.float)}
pid, peak_temp, peak_rho, fe, ni, ti = np.loadtxt(PLOT_FILE, datatype,
	delimiter=", ", skiprows=1, usecols=(i_pid, i_peak_temp, i_peak_rho, i_fe,
	i_56ni, i_44ti), unpack=True)
print "cco2 plotting data obtained"

# Calculate the Fe/Ti ratios and radiation entropies to identify regions
srad = peak_temp**3 / peak_rho
feti = (fe + ni) / ti
usable = np.isfinite(feti) & (feti != 0.)
pid, feti, srad = pid[usable], feti[usable], srad[usable]

# Make a list of the particle ID for every particle in each region
west_pids, east_pids, north_pids = [], [], []
for i in xrange(pid.size):
	if WEST_SR_MIN < srad[i] < WEST_SR_MAX \
		and WEST_FETI_MIN < feti[i] < WEST_FETI_MAX:
		west_pids.append(pid[i])
	elif EAST_SR_MIN < srad[i] < EAST_SR_MAX \
		and EAST_FETI_MIN < feti[i] < EAST_FETI_MAX:
		east_pids.append(pid[i])
	elif NORTH_SR_MIN < srad[i] < NORTH_SR_MAX \
		and NORTH_FETI_MIN < feti[i] < NORTH_FETI_MAX:
		north_pids.append(pid[i])
all_pids = []
all_pids.extend([(wpid, "W") for wpid in west_pids])
all_pids.extend([(epid, "E") for epid in east_pids])
all_pids.extend([(npid, "N") for npid in north_pids])
all_pids.sort(key=lambda x: x[0])

# Write out the list of particle IDs from all 3 regions
with open(PID_FILE, "w") as pid_file:
	head = "n_id %d\n" % (len(all_pids))
	pid_file.write(head)
	for apid, region in all_pids:
		line = "%d\n" % (apid)
		pid_file.write(line)
print "list of regions ids produced"

# Write out an alternate list of particle IDs tagged by region
with open(REGIONS_FILE, "w") as regions_file:
	head = "# pid region\n"
	regions_file.write(head)
	for apid, region in all_pids:
		line = "%d %s\n" % (apid, region)
		regions_file.write(line)
print "alternate list of regions ids produced"

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

