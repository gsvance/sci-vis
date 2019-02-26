#!/usr/bin/env python
# plot_radial_profiles.py

# Read the SDF data retrieved by the program get_sdf_data.c
# Make a series of radial plots for studying cco2 the shock breakout

# Last modified 2/26/19 by Greg Vance

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

# Iteration where shock breakout occurs -- make a nicer plot
BREAKOUT_ITER = 18400

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
	vrad /= 1e3 * KMPS
	rho = np.log10(rho)
	mass /= MSUN
	
	# Set up a big figure to plot two radial profiles together
	# The figsize here is twice the height of the MPL default
	plt.figure(figsize=(6.4, 9.6))
	
	# Upper subplot of vr vs r
	plt.subplot(211)
	plt.hist2d(rad, vrad, 100, ((6.5, 14.0), (-40.0, 55.0)), weights=mass,
		norm=clr.LogNorm(10.0**-5.5, 10.0**0.5))
	plt.xlim(6.5, 14.0)
	plt.ylim(-40.0, 55.0)
	plt.xlabel("log Radius [cm]")
	plt.ylabel("Radial Velocity (10$^3$ km/s)")
	plt.colorbar(label="Total Bin Mass (M$_\\odot$)")
	plt.title("Radial Profiles for Simulation cco2\n" +
		"Iteration %d (t = %g sec)" % (iter, tpos))
	
	# Lower subplot of rho vs r
	plt.subplot(212)
	plt.hist2d(rad, rho, 100, ((6.5, 14.0), (-10.5, 7.0)), weights=mass,
		norm=clr.LogNorm(10.0**-5.5, 10.0**0.5))
	plt.xlim(6.5, 14.0)
	plt.ylim(-10.5, 7.0)
	plt.xlabel("log Radius [cm]")
	plt.ylabel("log Density [g/cm$^3$]")
	plt.colorbar(label="Total Bin Mass (M$_\\odot$)")
	
	plt.savefig(PLOT_DIR + "profile_iter%06d.png" % (iter), dpi=120)
	plt.close()
	
	# Make some nicer plots at the shock breakout point
	if iter == BREAKOUT_ITER:
		print "data file %d matches breakout iter" % (i + 1)
		
		plt.figure(figsize=(6.4, 9.6))
		
		plt.subplot(211)
		plt.hist2d(rad, vrad, 100, weights=mass, norm=clr.LogNorm())
		plt.xlabel("log Radius [cm]")
		plt.ylabel("Radial Velocity (10$^3$ km/s)")
		plt.colorbar(label="Total Bin Mass (M$_\\odot$)")
		plt.title("Radial Profiles for Simulation cco2\nat Shock Breakout")
		
		plt.subplot(212)
		plt.hist2d(rad, rho, 100, weights=mass, norm=clr.LogNorm())
		plt.xlabel("log Radius [cm]")
		plt.ylabel("log Density [g/cm$^3$]")
		plt.colorbar(label="Total Bin Mass (M$_\\odot$)")
		
		plt.savefig(PLOT_DIR + "profile_breakout.png", dpi=120)
		plt.close()

