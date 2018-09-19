#!/usr/bin/env python

# The cco2 44Ti Magkotsios plot has a high-abundance yellow region
# Report the particle IDs of everything inside that yellow region
# ...

# Last modified 9/18/18 by Greg Vance

import numpy as np
import os
import sdfpy as sp
import matplotlib.pyplot as plt

# Data file locations on saguaro
DATADIR = "/home/gsvance/results/plotting/"
DATAFILE = DATADIR + "cco2_plotting.out"
COLSFILE = DATADIR + "columns"
SDFDIR = "/home/gsvance/supernova_data/cco2/cco2/"

# Boundaries (by eye) of the high-abundance "yellow region"
# See young/shared_files/magkotsios/cco2_magkotsios_hist.png
PEAKTEMPMIN = 8e9 # K
PEAKTEMPMAX = 10.5e9 # K
PEAKRHOMIN = 10.**6.5 # g cm^-3
PEAKRHOMAX = 10.**7.0 # g cm^-3

# Time unit for SNSPH tpos values
SNSPHTIME = 100. # s

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
yellow_pids = list()
for i in xrange(pid.size):
	if PEAKTEMPMIN < peaktemp[i] < PEAKTEMPMAX \
		and PEAKRHOMIN < peakrho[i] < PEAKRHOMAX:
		yellow_pids.append(pid[i])
yellow_pids.sort()

# Print the details of the yellow region IDs
print "ids in the yellow region:", len(yellow_pids)

# Form a list of all the cco2 SDF file names
sdf_names = list()
for maybe_sdf in os.listdir(SDFDIR):
	extension = maybe_sdf.split('.')[-1]
	try:
		ext = int(extension)
	except ValueError:
		continue
	sdf_name = SDFDIR + maybe_sdf
	sdf_names.append(sdf_name)
print "sdf files identified:", len(sdf_names)

# Binary search
def bin_search_in(x, l):
	lo, hi = 0, len(l) - 1
	while hi - lo >= 0:
		mid = (lo + hi) / 2
		if x < l[mid]:
			hi = mid - 1
		elif x > l[mid]:
			lo = mid + 1
		else:
			return True
	return False

# Compile the necessary trajectory data from the SDFs
sdfs_data = list()
for sdf_name in sdf_names:
	print sdf_name
	sdf = sp.SDFRead(sdf_name)
	dat = dict()
	dat["time"] = float(sdf.parameters["tpos"] / SNSPHTIME)
	for ident, rho, temp in zip(sdf["ident"], sdf["rho"], sdf["temp"]):
		if bin_search_in(ident, yellow_pids):
			id = int(ident)
			dat[id] = dict()
			dat[id]["rho"] = float(rho)
			dat[id]["temp"] = float(temp)
	sdfs_data.append(dat)
	del sdf
print "sdf data acquired"

# Set up arrays to organize the trajectory data for plotting
n_part = len(yellow_pids)
n_time = len(sdfs_data)
sdfs_data.sort(key=lambda d : d["time"])
times = np.array([d["time"] for d in sdfs_data])
rho_traj = np.empty((n_part, n_time))
temp_traj = np.empty((n_part, n_time))
print "sdf data organized for plotting"

# Fill the plotting arrays with data
for i_t in xrange(n_time):
	for i_p, ypid in enumerate(yellow_pids):
		rho_traj[i_p,i_t] = sdfs_data[i_t][ypid]["rho"]
		temp_traj[i_p,i_t] = sdfs_data[i_t][ypid]["temp"]

# Make the spaghetti plot for density
plt.figure()
for i_p in xrange(n_part):
	plt.plot(times, rho_traj[i_p], "g-")
plt.savefig("spaghetti_rho.png")
plt.close()

# Make the spaghetti plot for temperature
plt.figure()
for i_p in xrange(n_part):
	plt.plot(times, temp_traj[i_p], "b-")
plt.savefig("spaghetti_temp.png")
plt.close()

# 
rho_mean = np.mean(rho_traj, axis=0)
rho_sigma = np.std(rho_traj, axis=0)
temp_mean = np.mean(temp_traj, axis=0)
temp_sigma = np.std(temp_traj, axis=0)


