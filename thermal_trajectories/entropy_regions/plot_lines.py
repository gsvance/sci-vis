#!/usr/bin/env python

# Read all fetched data files and do some trajectory plotting work
# Make spaghetti plots of the three entropy regions in different colors
# Create cleaner plots showing the mean and spread of the trajectories
# Create plots showing best-fit Magkotsios trajectories for each region

# Last modified 10/10/18 by Greg Vance

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop

# Fetched data files for this script to read
TIME_FILE = "time_fetched.dat"
RHO_FILE = "rho_fetched.dat"
TEMP_FILE = "temp_fetched.dat"
VRAD_FILE = "vrad_fetched.dat"

# File with an alternate particle IDs list tagging them by region
REGIONS_FILE = "regions.out"

# Subdirectory to store all plot images
PLOT_DIR = "images/"

# Units from the SNSPH code
SNSPH_TIME = 100. # s
SNSPH_MASS = 1e-6 * 1.9889e33 # g
SNSPH_LENGTH = 6.955e10 # cm
SNSPH_DENSITY = SNSPH_MASS / SNSPH_LENGTH**3
SNSPH_VELOCITY = SNSPH_LENGTH / SNSPH_TIME

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

# Read the radial velocity data from file
vrad_id = np.loadtxt(VRAD_FILE, np.int, usecols=[0])
vrad = np.loadtxt(VRAD_FILE, np.float, usecols=range(1, n_time + 1))

# Read the region-tagged ID list from file
region_id = np.loadtxt(REGIONS_FILE, np.int, usecols=[0])
region = np.loadtxt(REGIONS_FILE, "S1", usecols=[1])

# Sanity checks on the lists of IDs
assert len(rho_id) == len(temp_id)
assert len(rho_id) == len(vrad_id)
assert len(rho_id) == len(region_id)
n_id = len(rho_id)
assert all(rho_id == temp_id)
assert all(rho_id == vrad_id)
assert all(rho_id == region_id)
id = np.copy(rho_id)
del rho_id, temp_id, vrad_id, region_id
assert len(region) == n_id

# Sanity checks on the rho, temp, and vrad data
assert rho.shape == (n_id, n_time)
assert temp.shape == (n_id, n_time)
assert vrad.shape == (n_id, n_time)

# Convert densities to g/cm^3
dens = rho * SNSPH_DENSITY

# Convert radial velocities to cm/s
rvel = vrad * SNSPH_VELOCITY

# Print the results of the sanity checks
print "fetched data read in and sanity-checked"
print "number of points in time:", n_time
print "number of particle ids:", n_id

############################################################
# SPAGHETTI PLOTS FOR DENS, TEMP, AND RVEL VS ITER
############################################################

# Plot of iteration vs density with all lines
plt.figure()
labeled = [False, False, False]
for i in xrange(n_id):
	r = {"W": 0, "E": 1, "N": 2}[region[i]]
	col = ["red", "yellowgreen", "blue"][r]
	if not labeled[r]:
		lab = ["West", "East", "North"][r]
		plt.plot(iter, dens[i], color=col, label=lab)
		labeled[r] = True
	else:
		plt.plot(iter, dens[i], color=col, alpha=0.05)
plt.yscale("log")
plt.legend()
plt.title("Densities of 3 Regions of Particles vs. Iteration")
plt.xlabel("SNSPH Iteration")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "iter_vs_dens_lines.png", dpi=150)
plt.close()

# Plot of iteration vs temperature with all lines
plt.figure()
labeled = [False, False, False]
for i in xrange(n_id):
	r = {"W": 0, "E": 1, "N": 2}[region[i]]
	col = ["red", "yellowgreen", "blue"][r]
	if not labeled[r]:
		lab = ["West", "East", "North"][r]
		plt.plot(iter, temp[i], color=col, label=lab)
		labeled[r] = True
	else:
		plt.plot(iter, temp[i], color=col, alpha=0.05)
plt.yscale("log")
plt.legend()
plt.title("Temperatures of 3 Regions of Particles vs. Iteration")
plt.xlabel("SNSPH Iteration")
plt.ylabel("Temperature (K)")
plt.savefig(PLOT_DIR + "iter_vs_temp_lines.png", dpi=150)
plt.close()

# Plot of iteration vs radial velocity with all lines
plt.figure()
labeled = [False, False, False]
for i in xrange(n_id):
	r = {"W": 0, "E": 1, "N": 2}[region[i]]
	col = ["red", "yellowgreen", "blue"][r]
	if not labeled[r]:
		lab = ["West", "East", "North"][r]
		plt.plot(iter, rvel[i], color=col, label=lab)
		labeled[r] = True
	else:
		plt.plot(iter, rvel[i], color=col, alpha=0.05)
plt.yscale("log")
plt.legend()
plt.title("Velocities of 3 Regions of Particles vs. Iteration")
plt.xlabel("SNSPH Iteration")
plt.ylabel("Radial Velocity (cm/s)")
plt.savefig(PLOT_DIR + "iter_vs_rvel_lines.png", dpi=150)
plt.close()

############################################################
# SPAGHETTI PLOTS FOR DENS, TEMP, AND RVEL VS TIME
############################################################

# Plot of time vs density with all lines
plt.figure()
labeled = [False, False, False]
for i in xrange(n_id):
	r = {"W": 0, "E": 1, "N": 2}[region[i]]
	col = ["red", "yellowgreen", "blue"][r]
	if not labeled[r]:
		lab = ["West", "East", "North"][r]
		plt.plot(time, dens[i], color=col, label=lab)
		labeled[r] = True
	else:
		plt.plot(time, dens[i], color=col, alpha=0.05)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.title("Densities of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "time_vs_dens_lines.png", dpi=150)
plt.close()

# Plot of time vs temperature with all lines
plt.figure()
labeled = [False, False, False]
for i in xrange(n_id):
	r = {"W": 0, "E": 1, "N": 2}[region[i]]
	col = ["red", "yellowgreen", "blue"][r]
	if not labeled[r]:
		lab = ["West", "East", "North"][r]
		plt.plot(time, temp[i], color=col, label=lab)
		labeled[r] = True
	else:
		plt.plot(time, temp[i], color=col, alpha=0.05)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.title("Temperatures of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.savefig(PLOT_DIR + "time_vs_temp_lines.png", dpi=150)
plt.close()

# Plot of time vs radial velocity with all lines
plt.figure()
labeled = [False, False, False]
for i in xrange(n_id):
	r = {"W": 0, "E": 1, "N": 2}[region[i]]
	col = ["red", "yellowgreen", "blue"][r]
	if not labeled[r]:
		lab = ["West", "East", "North"][r]
		plt.plot(time, rvel[i], color=col, label=lab)
		labeled[r] = True
	else:
		plt.plot(time, rvel[i], color=col, alpha=0.05)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.title("Velocities of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Radial Velocity (cm/s)")
plt.savefig(PLOT_DIR + "time_vs_rvel_lines.png", dpi=150)
plt.close()

############################################################
# CALCULATING MEAN AND SPREAD STATISTICS
############################################################

# Mean and std. dev. of density in log space at every time
mean_dens = 10.**np.mean(np.log10(dens), axis=0)
sigma_dens = 10.**np.std(np.log10(dens), axis=0)

# Mean and std. dev. of temperature in log space at every time
mean_temp = 10.**np.mean(np.log10(temp), axis=0)
sigma_temp = 10.**np.std(np.log10(temp), axis=0)

# Mean and std. dev. of radial velocity in log space at every time
mean_rvel = 10.**np.mean(np.log10(rvel), axis=0)
sigma_rvel = 10.**np.std(np.log10(rvel), axis=0)

# Make boolean arrays for each of the 3 regions and check sanity
west = (region == "W")
east = (region == "E")
north = (region == "N")
assert np.sum(west) + np.sum(east) + np.sum(north) == n_id

# Calculate density mean and std. dev. for the west region
mean_dens_west = 10.**np.mean(np.log10(dens[west]), axis=0)
sigma_dens_west = 10.**np.std(np.log10(dens[west]), axis=0)

# Calculate density mean and std. dev. for the east region
mean_dens_east = 10.**np.mean(np.log10(dens[east]), axis=0)
sigma_dens_east = 10.**np.std(np.log10(dens[east]), axis=0)

# Calculate density mean and std. dev. for the north region
mean_dens_north = 10.**np.mean(np.log10(dens[north]), axis=0)
sigma_dens_north = 10.**np.std(np.log10(dens[north]), axis=0)

# Calculate temperature mean and std. dev. for the west region
mean_temp_west = 10.**np.mean(np.log10(temp[west]), axis=0)
sigma_temp_west = 10.**np.std(np.log10(temp[west]), axis=0)

# Calculate temperature mean and std. dev. for the east region
mean_temp_east = 10.**np.mean(np.log10(temp[east]), axis=0)
sigma_temp_east = 10.**np.std(np.log10(temp[east]), axis=0)

# Calculate temperature mean and std. dev. for the north region
mean_temp_north = 10.**np.mean(np.log10(temp[north]), axis=0)
sigma_temp_north = 10.**np.std(np.log10(temp[north]), axis=0)

# Calculate radial velocity mean and std. dev. for the west region
mean_rvel_west = 10.**np.mean(np.log10(rvel[west]), axis=0)
sigma_rvel_west = 10.**np.std(np.log10(rvel[west]), axis=0)

# Calculate radial velocity mean and std. dev. for the east region
mean_rvel_east = 10.**np.mean(np.log10(rvel[east]), axis=0)
sigma_rvel_east = 10.**np.std(np.log10(rvel[east]), axis=0)

# Calculate radial velocity mean and std. dev. for the north region
mean_rvel_north = 10.**np.mean(np.log10(rvel[north]), axis=0)
sigma_rvel_north = 10.**np.std(np.log10(rvel[north]), axis=0)

############################################################
# PLOTTING MEAN TRAJECTORIES OF EACH REGION WITH SPREADS
############################################################

# Plot mean density trajectories of each region with spread
plt.figure()
for r in xrange(3):
	mu = [mean_dens_west, mean_dens_east, mean_dens_north][r]
	sig = [sigma_dens_west, sigma_dens_east, sigma_dens_north][r]
	lab = ["West", "East", "North"][r] + " $\\pm$ 1$\\sigma$"
	col = ["red", "yellowgreen", "blue"][r]
	plus = 10.**(np.log10(mu) + np.log10(sig))
	minus = 10.**(np.log10(mu) - np.log10(sig))
	plt.fill_between(time, plus, minus, color=col, alpha=0.2)
	plt.plot(time, mu, color=col, label=lab)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.title("Densities of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "time_vs_dens_sigma.png", dpi=150)
plt.close()

# Plot mean temperature trajectories of each region with spread
plt.figure()
for r in xrange(3):
	mu = [mean_temp_west, mean_temp_east, mean_temp_north][r]
	sig = [sigma_temp_west, sigma_temp_east, sigma_temp_north][r]
	lab = ["West", "East", "North"][r] + " $\\pm$ 1$\\sigma$"
	col = ["red", "yellowgreen", "blue"][r]
	plus = 10.**(np.log10(mu) + np.log10(sig))
	minus = 10.**(np.log10(mu) - np.log10(sig))
	plt.fill_between(time, plus, minus, color=col, alpha=0.2)
	plt.plot(time, mu, color=col, label=lab)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.title("Temperatures of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.savefig(PLOT_DIR + "time_vs_temp_sigma.png", dpi=150)
plt.close()

# Plot mean radial velocity trajectories of each region with spread
plt.figure()
for r in xrange(3):
	mu = [mean_rvel_west, mean_rvel_east, mean_rvel_north][r]
	sig = [sigma_rvel_west, sigma_rvel_east, sigma_rvel_north][r]
	lab = ["West", "East", "North"][r] + " $\\pm$ 1$\\sigma$"
	col = ["red", "yellowgreen", "blue"][r]
	plus = 10.**(np.log10(mu) + np.log10(sig))
	minus = 10.**(np.log10(mu) - np.log10(sig))
	plt.fill_between(time, plus, minus, color=col, alpha=0.2)
	plt.plot(time, mu, color=col, label=lab)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.title("Velocities of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Radial Velocity (cm/s)")
plt.savefig(PLOT_DIR + "time_vs_rvel_sigma.png", dpi=150)
plt.close()

############################################################
# DEFINING AND FITTING MAGKOTSIOS TRAJECTORIES
############################################################

# Exponential trajectories from Magkotsios Eq. 2
def temp_exp(t, T0, tau):
	return T0 * np.exp(-t / (3. * tau))
def dens_exp(t, rho0, tau):
	return rho0 * np.exp(-t / tau)

# Power law trajectories from Magkotsios Eq. 5
def temp_pow(t, T0):
	return T0 / (2. * t + 1.)
def dens_pow(t, rho0):
	return rho0 / (2. * t + 1.)**3

# Power law trajectories with slope parameter alpha
def temp_pow2(t, T0, alpha):
	return T0 * (2. * t + 1.)**(-alpha)
def dens_pow2(t, rho0, alpha):
	return rho0 * (2. * t + 1.)**(-3. * alpha)

# Define a general function to do the exponential trajectory fitting
def exp_fitter(dens_data, temp_data):
	fixed_rho0 = dens_data[0]
	fixed_T0 = temp_data[0]
	def dens_exp_fit(t, tau):
		return np.log10(dens_exp(t, fixed_rho0, tau))
	def temp_exp_fit(t, tau):
		return np.log10(temp_exp(t, fixed_T0, tau))
	dopt, dcov = spop.curve_fit(dens_exp_fit, time, np.log10(dens_data), [1e3])
	topt, tcov = spop.curve_fit(temp_exp_fit, time, np.log10(temp_data), [1e3])
	dtup = tuple([fixed_rho0] + list(dopt))
	ttup = tuple([fixed_T0] + list(topt))
	return dtup, ttup

# Fit exponential trajectories in log space with only the tau parameter
dens_exp_best_west, temp_exp_best_west = exp_fitter(mean_dens_west,
	mean_temp_west)
dens_exp_best_east, temp_exp_best_east = exp_fitter(mean_dens_east,
	mean_temp_east)
dens_exp_best_north, temp_exp_best_north = exp_fitter(mean_dens_north,
	mean_temp_north)

# "Fit" the power law trajectories with no real fitting parameters
dens_pow_best_west = tuple([mean_dens_west[0]])
temp_pow_best_west = tuple([mean_temp_west[0]])
dens_pow_best_east = tuple([mean_dens_east[0]])
temp_pow_best_east = tuple([mean_temp_east[0]])
dens_pow_best_north = tuple([mean_dens_north[0]])
temp_pow_best_north = tuple([mean_temp_north[0]])

# Define a general function to do the sloped power law trajectory fitting
def pow2_fitter(dens_data, temp_data):
	fixed_rho0 = dens_data[0]
	fixed_T0 = temp_data[0]
	def dens_pow2_fit(t, alpha):
		return np.log10(dens_pow2(t, fixed_rho0, alpha))
	def temp_pow2_fit(t, alpha):
		return np.log10(temp_pow2(t, fixed_T0, alpha))
	dopt, dcov = spop.curve_fit(dens_pow2_fit, time, np.log10(dens_data), [1.])
	topt, tcov = spop.curve_fit(temp_pow2_fit, time, np.log10(temp_data), [1.])
	dtup = tuple([fixed_rho0] + list(dopt))
	ttup = tuple([fixed_T0] + list(topt))
	return dtup, ttup

# Fit power law trajectories in log space with added alpha slope parameter
dens_pow2_best_west, temp_pow2_best_west = pow2_fitter(mean_dens_west,
	mean_temp_west)
dens_pow2_best_east, temp_pow2_best_east = pow2_fitter(mean_dens_east,
	mean_temp_east)
dens_pow2_best_north, temp_pow2_best_north = pow2_fitter(mean_dens_north,
	mean_temp_north)

############################################################
# REGIONAL SPAGHETTI PLOTS FOR DENS, TEMP, AND RVEL VS TIME
############################################################

# Plots of density vs time in each region with fitted trajectories
for r in xrange(3):
	tag = ["W", "E", "N"][r]
	lab = ["West", "East", "North"][r]
	col = ["red", "yellowgreen", "blue"][r]
	exp_best = [dens_exp_best_west, dens_exp_best_east, dens_exp_best_north][r]
	pow_best = [dens_pow_best_west, dens_pow_best_east, dens_pow_best_north][r]
	pow2_best = [dens_pow2_best_west, dens_pow2_best_east,
		dens_pow2_best_north][r]
	plt.figure()
	for i in xrange(n_id):
		if region[i] == tag:
			plt.plot(time, dens[i], color=col, alpha=0.05)
	plt.plot(time, dens_exp(time, *exp_best), "k--", label="Exponential")
	plt.plot(time, dens_pow(time, *pow_best), "k:", label="Power Law")
	plt.plot(time, dens_pow2(time, *pow2_best), "k-.", label="Power Law 2")
	plt.xscale("log")
	plt.yscale("log")
	plt.legend()
	plt.title("Densities of %s Region Particles vs. Time" % (lab))
	plt.xlabel("Time (s)")
	plt.ylabel("Density (g/cm$^3$)")
	plotfile = PLOT_DIR + "time_vs_dens_%s.png" % (lab.lower())
	plt.savefig(plotfile, dpi=150)
	plt.close()

# Plots of temperature vs time in each region with fitted trajectories
for r in xrange(3):
	tag = ["W", "E", "N"][r]
	lab = ["West", "East", "North"][r]
	col = ["red", "yellowgreen", "blue"][r]
	exp_best = [temp_exp_best_west, temp_exp_best_east, temp_exp_best_north][r]
	pow_best = [temp_pow_best_west, temp_pow_best_east, temp_pow_best_north][r]
	pow2_best = [temp_pow2_best_west, temp_pow2_best_east,
		temp_pow2_best_north][r]
	plt.figure()
	for i in xrange(n_id):
		if region[i] == tag:
			plt.plot(time, temp[i], color=col, alpha=0.05)
	plt.plot(time, temp_exp(time, *exp_best), "k--", label="Exponential")
	plt.plot(time, temp_pow(time, *pow_best), "k:", label="Power Law")
	plt.plot(time, temp_pow2(time, *pow2_best), "k-.", label="Power Law 2")
	plt.xscale("log")
	plt.yscale("log")
	plt.legend()
	plt.title("Temperatures of %s Region Particles vs. Time" % (lab))
	plt.xlabel("Time (s)")
	plt.ylabel("Temperature (K)")
	plotfile = PLOT_DIR + "time_vs_temp_%s.png" % (lab.lower())
	plt.savefig(plotfile, dpi=150)
	plt.close()

# Plots of radial velocity vs time in each region without fitted trajectories
for r in xrange(3):
	tag = ["W", "E", "N"][r]
	lab = ["West", "East", "North"][r]
	col = ["red", "yellowgreen", "blue"][r]
	#exp_best = [dens_exp_best_west, dens_exp_best_east, dens_exp_best_north][r]
	#pow_best = [dens_pow_best_west, dens_pow_best_east, dens_pow_best_north][r]
	#pow2_best = [dens_pow2_best_west, dens_pow2_best_east,
	#	dens_pow2_best_north][r]
	plt.figure()
	for i in xrange(n_id):
		if region[i] == tag:
			plt.plot(time, rvel[i], color=col, alpha=0.05)
	#plt.plot(time, dens_exp(time, *exp_best), "k--", label="Exponential")
	#plt.plot(time, dens_pow(time, *pow_best), "k:", label="Power Law")
	#plt.plot(time, dens_pow2(time, *pow2_best), "k-.", label="Power Law 2")
	plt.xscale("log")
	plt.yscale("log")
	#plt.legend()
	plt.title("Velocities of %s Region Particles vs. Time" % (lab))
	plt.xlabel("Time (s)")
	plt.ylabel("Radial Velocity (cm/s)")
	plotfile = PLOT_DIR + "time_vs_rvel_%s.png" % (lab.lower())
	plt.savefig(plotfile, dpi=150)
	plt.close()

