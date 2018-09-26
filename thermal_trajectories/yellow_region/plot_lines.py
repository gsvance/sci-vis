#!/usr/bin/env python

# Read the fetched data files and do the actual plotting work
# Make spaghetti plots of the yellow region particle trajectories
# Also produce cleaner plots with the visual spread statistics

# Last modified 9/21/18 by Greg Vance

import numpy as np
import matplotlib.pyplot as plt

# Fetched data files for this script to read
TIME_FILE = "time_fetched.dat"
RHO_FILE = "rho_fetched.dat"
TEMP_FILE = "temp_fetched.dat"

# Subdirectory to store all plot images
PLOT_DIR = "images/"

# Units from the SNSPH code
SNSPH_TIME = 100. # s
SNSPH_MASS = 1e-6 * 1.9889e33 # g
SNSPH_LENGTH = 6.955e10 # cm
SNSPH_DENSITY = SNSPH_MASS / SNSPH_LENGTH**3

# Read the iter and tpos data points from file
iter = np.loadtxt(TIME_FILE, np.int, usecols=[0])
tpos = np.loadtxt(TIME_FILE, np.float, usecols=[1])

# Sanity check on the time data
assert len(iter) == len(tpos)
n_time = len(tpos)

# Convert times to seconds
time = tpos * SNSPH_TIME

# Read the density data from file
rho_id = np.loadtxt(RHO_FILE, np.int, usecols=[0])
rho = np.loadtxt(RHO_FILE, np.float, usecols=range(1, n_time + 1))

# Read the temperature data from file
temp_id = np.loadtxt(TEMP_FILE, np.int, usecols=[0])
temp = np.loadtxt(TEMP_FILE, np.float, usecols=range(1, n_time + 1))

# Sanity checks on the lists of IDs
assert len(rho_id) == len(temp_id)
n_id = len(rho_id)
assert all(rho_id == temp_id)
id = rho_id

# Sanity checks on the rho and temp data
assert rho.shape == (n_id, n_time)
assert temp.shape == (n_id, n_time)

# Convert densities to g/cm^3
dens = rho * SNSPH_DENSITY

# Print the results of the sanity checks
print "fetched data read in and sanity-checked"
print "number of points in time:", n_time
print "number of particle ids:", n_id

############################################################
# SPAGHETTI PLOTS FOR DENS AND TEMP VS ITER
############################################################

# Plot of iteration vs density with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(iter, dens[i], "b-", alpha=0.2)
plt.yscale("log")
plt.title("Densities of Yellow Region Particles vs. Iteration")
plt.xlabel("SNSPH Iteration")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "iter_vs_dens_lines.png", dpi=150)
plt.close()

# Plot of iteration vs temperature with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(iter, temp[i], "r-", alpha=0.2)
plt.yscale("log")
plt.title("Temperatures of Yellow Region Particles vs. Iteration")
plt.xlabel("SNSPH Iteration")
plt.ylabel("Temperature (K)")
plt.savefig(PLOT_DIR + "iter_vs_temp_lines.png", dpi=150)
plt.close()

############################################################
# SPAGHETTI PLOTS FOR DENS AND TEMP VS TIME
############################################################

# Plot of time vs density with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(time, dens[i], "b-", alpha=0.2)
plt.xscale("log")
plt.yscale("log")
plt.title("Densities of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "time_vs_dens_lines.png", dpi=150)
plt.close()

# Plot of time vs temperature with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(time, temp[i], "r-", alpha=0.2)
plt.xscale("log")
plt.yscale("log")
plt.title("Temperatures of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.savefig(PLOT_DIR + "time_vs_temp_lines.png", dpi=150)
plt.close()

