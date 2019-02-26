#!/usr/bin/env python
# plot_radial_profiles.py

# Read the SDF data retrieved by the program get_sdf_data.c
# Make a series of radial plots for studying cco2 the shock breakout

# Last modified 2/25/19 by Greg Vance

# Import packagess (these can take time on saguaro...)
print "starting imports"
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
print "imports completed"

# Directory with the SDF data to be used
DATA_DIR = "data_from_sdfs/"

# Output directory for completed images
PLOT_DIR = "radial_profiles/"

# Units from SNSPH for converting to CGS
SNSPH_TIME = 100. # s
SNSPH_MASS = 1e-6 * 1.9889e33 # g
SNSPH_LENGTH = 6.955e10 # cm
SNSPH_DENSITY = SNSPH_MASS / SNSPH_LENGTH**3
SNSPH_VELOCITY = SNSPH_LENGTH / SNSPH_TIME

# Useful units for axis scales
MSUN = 1.9889e33 # g
ASTU = 1.495978707e13 # cm
KMPS = 1e5 # cm/s

rmin, rmax = 8.0, 13.0
vmin, vmax = 6.0, 8.0
dmin, dmax = -5.0, 2.0

print "beginning to read data files"
data_file_list = sorted(os.listdir(DATA_DIR))
n_files = len(data_file_list)
for i, data_file_name in enumerate(data_file_list):
	
	print "processing data file %d/%d" % (i + 1, n_files)
	data_file_name = DATA_DIR + data_file_name
	
	# Read the iter and tpos values from this file
	with open(data_file_name, "r") as data_file:
		iter = int(data_file.readline().strip().split()[1])
		tpos = float(data_file.readline().strip().split()[1])
	
	# Instruct numpy to load the float data from this file
	rad, vrad, rho, mass = np.loadtxt(data_file_name, np.float32, skiprows=3,
		unpack=True)
	
	# Skip the few "impostor" SDF files at the start of cco2
	if all(rad == 0.0):
		print "data file %d was skipped" % (i + 1)
		continue
	
	# Convert the SDF data to CGS units
	tpos *= SNSPH_TIME
	rad *= SNSPH_LENGTH
	vrad *= SNSPH_VELOCITY
	rho *= SNSPH_DENSITY
	mass *= SNSPH_MASS
	
	# Convert the CGS data to units convenient for the plot axes
	rad = np.log10(rad)
	vrad = np.log10(vrad)
	rho = np.log10(rho)
	mass /= MSUN
	
	if min(rad) < rmin:
		rmin = min(rad)
	if max(rad) > rmax:
		rmax = max(rad)
	if min(vrad) < vmin:
		vmin = min(vrad)
	if max(vrad) > vmax:
		vmax = max(vrad)
	if min(rho) < dmin:
		dmin = min(rho)
	if max(rho) > dmax:
		dmax = max(rho)
	
	# Set up a big figure to plot two radial profiles together
	# The figsize here is twice the height of the MPL default
	plt.figure(figsize=(6.4, 9.6))
	
	# Upper subplot of vr vs r
	plt.subplot(211)
	plt.hist2d(rad, vrad, 100, weights=mass, norm=clr.LogNorm(1e-6, 1e0))
	plt.xlim(6.5, 14.0)
	plt.ylim(3.5, 10.0)
	plt.xlabel("log Radius [cm]")
	plt.ylabel("log Radial Velocity [cm/s]")
	plt.colorbar(label="Total Bin Mass (M$_\\odot$)")
	plt.title("Radial Profiles for Simulation cco2\n" +
		"Iteration %d (t = %g sec)" % (iter, tpos))
	
	# Lower subplot of rho vs r
	plt.subplot(212)
	plt.hist2d(rad, rho, 100, weights=mass, norm=clr.LogNorm(1e-6, 1e0))
	plt.xlim(6.5, 14.0)
	plt.ylim(-10.5, 7.0)
	plt.xlabel("log Radius [cm]")
	plt.ylabel("log Density [g/cm$^3$]")
	plt.colorbar(label="Total Bin Mass (M$_\\odot$)")
	
	plt.savefig(PLOT_DIR + "profile_iter%06d.png" % (iter), dpi=120)
	plt.close()

print "rad", rmin, rmax
print "vrad", vmin, vmax
print "rho", dmin, dmax

