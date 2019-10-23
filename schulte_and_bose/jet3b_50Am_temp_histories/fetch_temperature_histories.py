#!/usr/bin/env python
# Fetch the temperature histories for particles as described in README

# Last modified 10/22/19 by Greg Vance

import os
import sdfpy as sp
import numpy as np

# Relevant data directories
SDF_DIR_50AM = "/home/gsvance/supernova_data/50Am/run3g_50Am6/"
SDF_DIR_JET3B = "/home/gsvance/supernova_data/jet3b/jet3b/"
PLOTTING_DIR = "/home/gsvance/results/plotting/"

# Code units from SNSPH
SNSPH_MASS = 1e-6 * 1.9889e33 # g
SNSPH_LENGTH = 6.955e10 # cm
SNSPH_TIME = 100.0 # s
SNSPH_DENSITY = SNSPH_MASS / SNSPH_LENGTH**3
SNSPH_VELOCITY = SNSPH_LENGTH / SNSPH_TIME

# Which columns of the plotting files are PID, 12C, and 13C?
cols_file = os.path.join(PLOTTING_DIR, "columns")
with open(cols_file, "r") as cf:
	cols = [line.strip() for line in cf if line.strip() != ""]
i_pid = cols.index("id")
i_12c = cols.index("X_{12C}")
i_13c = cols.index("X_{13C}")
my_cols = (i_pid, i_12c, i_13c)

# Function to extract file extensions
def get_extension(filename):
	return filename.split('.')[-1]

# Function to get the tpos value from an SDF
def get_tpos(filename):
	sdf = sp.SDFRead(filename)
	tpos = sdf.parameters["tpos"]
	del sdf
	return tpos

# Organize the string info for each simulation
info = list()
info.append(("50Am", SDF_DIR_50AM, \
	os.path.join(PLOTTING_DIR, "50Am_plotting.out")))
info.append(("jet3b", SDF_DIR_JET3B, \
	os.path.join(PLOTTING_DIR, "jet3b_plotting.out")))

# Loop over the two simulations
for sim_name, sdf_dir, plotting_file in info:
	
	# Read the PID, 12C, and 13C data from sim's plotting file
	dtype = {"names": ("PID", "12C", "13C"), \
		"formats": (np.int, np.float, np.float)}
	plot_data = np.loadtxt(plotting_file, dtype, delimiter=", ", skiprows=1, \
		usecols=my_cols)
	
	# Identify the particles that have 12C/13C < 20
	select = (plot_data["13C"] != 0.) \
		& (plot_data["12C"] < 20. * plot_data["13C"])
	selected_pids = plot_data["PID"][select]
	n_pids = len(selected_pids)
	
	# Get a list of all the sim's SDFs sorted by tpos value
	all_files = os.listdir(sdf_dir)
	all_sdfs = [f for f in all_files if get_extension(f).isdigit()]
	sdfs_with_tpos = list()
	for filename in all_sdfs:
		fullname = os.path.join(sdf_dir, filename)
		tpos = get_tpos(fullname)
		sdfs_with_tpos.append((tpos, fullname))
	sdfs_with_tpos.sort(key=lambda x : x[0])
	n_sdfs = len(sdfs_with_tpos)
	
	# Set up arrays to store the temperature history data etc.
	iter_arr = np.empty(n_sdfs, np.int)
	tpos_arr = np.empty(n_sdfs, np.float)
	temp_arr = np.empty((n_pids, n_sdfs), np.float)
	dens_arr = np.empty((n_pids, n_sdfs), np.float)
	
	# Read all of the desired data from the sim's SDFs
	for j, sdfname in enumerate([x[1] for x in sdfs_with_tpos]):
		
		# Open the SDF and read the time information
		sdf = sp.SDFRead(sdfname)
		iter_arr[j] = sdf.parameters["iter"]
		tpos_arr[j] = sdf.parameters["tpos"]
		
		# Read the data for the selected particles
		i = 0
		for sdf_i, sdf_pid in enumerate(sdf["ident"]):
			if sdf_pid in selected_pids:
				assert sdf_pid == selected_pids[i]
				temp_arr[i,j] = sdf["temp"][sdf_i]
				dens_arr[i,j] = sdf["rho"][sdf_i]
				i += 1
		
		# Make sure the SDF file is closed
		del sdf
	
	# Convert all of the data from code units to CGS
	time = tpos_arr * SNSPH_TIME
	temp = temp_arr * 1.0
	dens = dens_arr * SNSPH_DENSITY
	
	# Write out all of the data to files
	outfile = "%s_histories_12C13C.dat" % (sim_name)
	with open(outfile, "w") as of:
		for j in xrange(n_sdfs):
			of.write("iteration %d, time %s\n" % (iter_arr[j], time[j]))
			of.write("id, temp, density\n")
			for i in xrange(n_pids):
				of.write("%d, %s, %s\n" % (selected_pids[i], temp[i,j], dens[i,j]))

