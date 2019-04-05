#!/usr/bin/env python

# The cco2 histogram of Fe/Ti vs radiation entropy displays 4 distinct regions:
#  - Region 1: bright yellow blob appearing at log(S_rad) ~ 22.5
#  - Region 2: smeared patch of reddish bins up at log(Fe/Ti) ~ 3.5
#  - Region 3: crescent shape with bright yellow bins at log(S_rad) ~ 22.6
#  - Region 4: extended region of high-entropy particles at log(S_rad) ~ 23.1
# Save a list of particle IDs for everything within these 4 regions
# Save a list of all the cco2 SDF files in order with tpos values
# These lists will be read by fetch_data.c when it does its job next
# Also write a file listing particle IDs by region for the plots later

# Last modified 4/5/19 by Greg Vance

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
REG1_SRAD_MIN = 10.**22.47 # K^3 cm^3 g^-1
REG1_SRAD_MAX = 10.**22.55 # K^3 cm^3 g^-1
REG1_FETI_MIN = 10.**2.4
REG1_FETI_MAX = 10.**2.7
REG2_SRAD_MIN = 10.**22.55 # K^3 cm^3 g^-1
REG2_SRAD_MAX = 10.**22.63 # K^3 cm^3 g^-1
REG2_FETI_MIN = 10.**3.5
REG2_FETI_MAX = 10.**4.1
REG3_SRAD_MIN = 10.**22.56 # K^3 cm^3 g^-1
REG3_SRAD_MAX = 10.**22.70 # K^3 cm^3 g^-1
REG3_FETI_MIN = 10.**2.2
REG3_FETI_MAX = 10.**2.8
REG4_SRAD_MIN = 10.**23.02 # K^3 cm^3 g^-1
REG4_SRAD_MAX = 10.**23.16 # K^3 cm^3 g^-1
REG4_FETI_MIN = 10.**0.9
REG4_FETI_MAX = 10.**4.1

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
reg1_pids, reg2_pids, reg3_pids, reg4_pids = [], [], [], []
for i in xrange(pid.size):
	if REG1_SRAD_MIN < srad[i] < REG1_SRAD_MAX \
		and REG1_FETI_MIN < feti[i] < REG1_FETI_MAX:
		reg1_pids.append(pid[i])
	elif REG2_SRAD_MIN < srad[i] < REG2_SRAD_MAX \
		and REG2_FETI_MIN < feti[i] < REG2_FETI_MAX:
		reg2_pids.append(pid[i])
	elif REG3_SRAD_MIN < srad[i] < REG3_SRAD_MAX \
		and REG3_FETI_MIN < feti[i] < REG3_FETI_MAX:
		reg3_pids.append(pid[i])
	elif REG4_SRAD_MIN < srad[i] < REG4_SRAD_MAX \
		and REG4_FETI_MIN < feti[i] < REG4_FETI_MAX:
		reg4_pids.append(pid[i])
all_pids = []
all_pids.extend([(pid1, "R1") for pid1 in reg1_pids])
all_pids.extend([(pid2, "R2") for pid2 in reg2_pids])
all_pids.extend([(pid3, "R3") for pid3 in reg3_pids])
all_pids.extend([(pid4, "R4") for pid4 in reg4_pids])
all_pids.sort(key=lambda x: x[0])

# Write out the list of particle IDs from all 4 regions
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

