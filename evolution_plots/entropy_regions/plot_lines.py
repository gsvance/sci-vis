#!/usr/bin/env python

# Read all fetched data files and do some trajectory plotting work
# Make spaghetti plots of the three entropy regions in different colors
# Create cleaner plots showing the mean and spread of the trajectories
# Create plots showing best-fit Magkotsios trajectories for each region

# Last modified 11/8/18 by Greg Vance

# Run the various import statements, which can sometimes take a while
print "starting imports"
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop
print "imports complete"

# Fetched data files for this script to read
TIME_FILE = "time_fetched.dat"
RHO_FILE = "rho_fetched.dat"
TEMP_FILE = "temp_fetched.dat"
VRAD_FILE = "vrad_fetched.dat"
TI44_FILE = "ti44_fetched.dat"
NI56_FILE = "ni56_fetched.dat"

# File with an alternate particle IDs list tagging them by region
REGIONS_FILE = "regions.out"

# Subdirectory to store all plot images
PLOT_DIR = "images/"

# Subdirectory to store the mean trajectories
MEAN_DIR = "mean_traj/"

# Units from the SNSPH code
SNSPH_TIME = 100. # s
SNSPH_MASS = 1e-6 * 1.9889e33 # g
SNSPH_LENGTH = 6.955e10 # cm
SNSPH_DENSITY = SNSPH_MASS / SNSPH_LENGTH**3
SNSPH_VELOCITY = SNSPH_LENGTH / SNSPH_TIME

# Alert for start of data read process
print "beginning to read data files"

# Read the iter and tpos data points from file
iter = np.loadtxt(TIME_FILE, np.int, usecols=[0])
tpos = np.loadtxt(TIME_FILE, np.float, usecols=[1])

# Sanity check on the time data
assert len(iter) == len(tpos)
n_time = len(tpos)

# Convert times to seconds
time = tpos * SNSPH_TIME

# Progress printing
print "time data read"

# Read the density data from file
rho_id = np.loadtxt(RHO_FILE, np.int, usecols=[0])
rho = np.loadtxt(RHO_FILE, np.float, usecols=range(1, n_time + 1))
print "density data read"

# Read the temperature data from file
temp_id = np.loadtxt(TEMP_FILE, np.int, usecols=[0])
temp = np.loadtxt(TEMP_FILE, np.float, usecols=range(1, n_time + 1))
print "temperature data read"

# Read the radial velocity data from file
vrad_id = np.loadtxt(VRAD_FILE, np.int, usecols=[0])
vrad = np.loadtxt(VRAD_FILE, np.float, usecols=range(1, n_time + 1))
print "radial velocity data read"

# Read the 44Ti abundance data from file
ti44_id = np.loadtxt(TI44_FILE, np.int, usecols=[0])
ti44 = np.loadtxt(TI44_FILE, np.float, usecols=range(1, n_time + 1))
print "44Ti data read"

# Read the 56Ni abundance data from file
ni56_id = np.loadtxt(NI56_FILE, np.int, usecols=[0])
ni56 = np.loadtxt(NI56_FILE, np.float, usecols=range(1, n_time + 1))
print "56Ni data read"

# Read the region-tagged ID list from file
region_id = np.loadtxt(REGIONS_FILE, np.int, usecols=[0])
region = np.loadtxt(REGIONS_FILE, "S1", usecols=[1])
print "region tag data read"

# Sanity checks on the lists of IDs
assert len(rho_id) == len(temp_id)
assert len(rho_id) == len(vrad_id)
assert len(rho_id) == len(ti44_id)
assert len(rho_id) == len(ni56_id)
assert len(rho_id) == len(region_id)
n_id = len(rho_id)
assert all(rho_id == temp_id)
assert all(rho_id == vrad_id)
assert all(rho_id == ti44_id)
assert all(rho_id == ni56_id)
assert all(rho_id == region_id)
id = np.copy(rho_id)
del rho_id, temp_id, vrad_id, ti44_id, ni56_id, region_id
assert len(region) == n_id

# Sanity checks on the rho, temp, and vrad data
assert rho.shape == (n_id, n_time)
assert temp.shape == (n_id, n_time)
assert vrad.shape == (n_id, n_time)
assert ti44.shape == (n_id, n_time)
assert ni56.shape == (n_id, n_time)

# Convert densities to g/cm^3
dens = rho * SNSPH_DENSITY

# Convert radial velocities to cm/s
rvel = vrad * SNSPH_VELOCITY

# Print the results of the sanity checks
print "all data sanity-checked"
print "number of points in time:", n_time
print "number of particle ids:", n_id

############################################################
# SPAGHETTI PLOTS FOR DENS, TEMP, AND RVEL VS ITER
############################################################

# Plot of iteration vs density with all lines
print "plotting dens vs iter spaghetti plot"
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
print "plotting temp vs iter spaghetti plot"
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
print "plotting vrad vs iter spaghetti plot"
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
print "plotting dens vs time spaghetti plot"
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
print "plotting temp vs time spaghetti plot"
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
print "plotting vrad vs time spaghetti plot"
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
plt.xlim(1e0, 1e5)
plt.legend()
plt.title("Velocities of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Radial Velocity (cm/s)")
plt.savefig(PLOT_DIR + "time_vs_rvel_lines.png", dpi=150)
plt.close()

############################################################
# CALCULATING MEAN AND SPREAD STATISTICS
############################################################

# Alert for user
print "calculating mean and spread stats"

# Mean and std. dev. of density in log space at every time
#mean_dens = 10.**np.mean(np.log10(dens), axis=0)
#sigma_dens = 10.**np.std(np.log10(dens), axis=0)

# Mean and std. dev. of temperature in log space at every time
#mean_temp = 10.**np.mean(np.log10(temp), axis=0)
#sigma_temp = 10.**np.std(np.log10(temp), axis=0)

# Mean and std. dev. of radial velocity in log space at every time
#mean_rvel = 10.**np.mean(np.log10(rvel), axis=0)
#sigma_rvel = 10.**np.std(np.log10(rvel), axis=0)

# Mean and std. dev. of 44Ti abundance in log space at every time
#mean_ti44 = 10.**np.mean(np.log10(ti44), axis=0)
#sigma_ti44 = 10.**np.std(np.log10(ti44), axis=0)

# Mean and std. dev. of 56Ni abundance in log space at every time
#mean_ni56 = 10.**np.mean(np.log10(ni56), axis=0)
#sigma_ni56 = 10.**np.std(np.log10(ni56), axis=0)

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

# Write mean density and temperature trajectories to file
np.savetxt(MEAN_DIR + "time_dens_temp_west.txt", np.column_stack((time,
	mean_dens_west, mean_temp_west)), header="time dens temp")
np.savetxt(MEAN_DIR + "time_dens_temp_east.txt", np.column_stack((time,
	mean_dens_east, mean_temp_east)), header="time dens temp")
np.savetxt(MEAN_DIR + "time_dens_temp_north.txt", np.column_stack((time,
	mean_dens_north, mean_temp_north)), header="time dens temp")

# Calculate radial velocity mean and std. dev. for the west region
mean_rvel_west = 10.**np.mean(np.log10(rvel[west]), axis=0)
sigma_rvel_west = 10.**np.std(np.log10(rvel[west]), axis=0)

# Calculate radial velocity mean and std. dev. for the east region
mean_rvel_east = 10.**np.mean(np.log10(rvel[east]), axis=0)
sigma_rvel_east = 10.**np.std(np.log10(rvel[east]), axis=0)

# Calculate radial velocity mean and std. dev. for the north region
mean_rvel_north = 10.**np.mean(np.log10(rvel[north]), axis=0)
sigma_rvel_north = 10.**np.std(np.log10(rvel[north]), axis=0)

# Calculate 44Ti abundance mean and std. dev. for the west region
mean_ti44_west = 10.**np.mean(np.log10(ti44[west]), axis=0)
sigma_ti44_west = 10.**np.std(np.log10(ti44[west]), axis=0)

# Calculate 44Ti abundance mean and std. dev. for the east region
mean_ti44_east = 10.**np.mean(np.log10(ti44[east]), axis=0)
sigma_ti44_east = 10.**np.std(np.log10(ti44[east]), axis=0)

# Calculate 44Ti abundance mean and std. dev. for the north region
mean_ti44_north = 10.**np.mean(np.log10(ti44[north]), axis=0)
sigma_ti44_north = 10.**np.std(np.log10(ti44[north]), axis=0)

# Calculate 56Ni abundance mean and std. dev. for the west region
mean_ni56_west = 10.**np.mean(np.log10(ni56[west]), axis=0)
sigma_ni56_west = 10.**np.std(np.log10(ni56[west]), axis=0)

# Calculate 56Ni abundance mean and std. dev. for the east region
mean_ni56_east = 10.**np.mean(np.log10(ni56[east]), axis=0)
sigma_ni56_east = 10.**np.std(np.log10(ni56[east]), axis=0)

# Calculate 56Ni abundance mean and std. dev. for the north region
mean_ni56_north = 10.**np.mean(np.log10(ni56[north]), axis=0)
sigma_ni56_north = 10.**np.std(np.log10(ni56[north]), axis=0)

############################################################
# PLOTTING MEAN TRAJECTORIES OF EACH REGION WITH SPREADS
############################################################

# Plot mean density trajectories of each region with spread
print "plotting mean densities"
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
plt.xlim(3e-2, 10.)
plt.ylim(1e3, 1e7)
plt.legend(loc="upper right")
plt.savefig(PLOT_DIR + "time_vs_dens_sigma_zoom.png", dpi=150)
plt.close()

# Plot mean temperature trajectories of each region with spread
print "plotting mean temperatures"
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
plt.xlim(3e-2, 10.)
plt.ylim(1e8, 1e10)
plt.legend(loc="upper right")
plt.savefig(PLOT_DIR + "time_vs_temp_sigma_zoom.png", dpi=150)
plt.close()

# Plot mean radial velocity trajectories of each region with spread
print "plotting mean radial velocities"
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
plt.xlim(1e0, 1e5)
plt.legend()
plt.title("Velocities of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Radial Velocity (cm/s)")
plt.savefig(PLOT_DIR + "time_vs_rvel_sigma.png", dpi=150)
plt.xlim(1e0, 10.)
plt.legend(loc="lower right")
plt.savefig(PLOT_DIR + "time_vs_rvel_sigma_zoom.png", dpi=150)
plt.close()

# Plot mean 44Ti abundance trajectories of each region with spread
print "plotting mean 44Ti abundances"
plt.figure()
for r in xrange(3):
	mu = [mean_ti44_west, mean_ti44_east, mean_ti44_north][r]
	sig = [sigma_ti44_west, sigma_ti44_east, sigma_ti44_north][r]
	lab = ["West", "East", "North"][r] + " $\\pm$ 1$\\sigma$"
	col = ["red", "yellowgreen", "blue"][r]
	plus = 10.**(np.log10(mu) + np.log10(sig))
	minus = 10.**(np.log10(mu) - np.log10(sig))
	plt.fill_between(time, plus, minus, color=col, alpha=0.2)
	plt.plot(time, mu, color=col, label=lab)
plt.xscale("log")
plt.yscale("log")
plt.xlim(1e0, 1e5)
plt.legend(loc="lower left")
plt.title("${}^{44}$Ti Abundances of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("${}^{44}$Ti Mass Fraction")
plt.savefig(PLOT_DIR + "time_vs_ti44_sigma.png", dpi=150)
plt.xlim(1e0, 10.)
#plt.legend(loc="lower right")
plt.savefig(PLOT_DIR + "time_vs_ti44_sigma_zoom.png", dpi=150)
plt.close()

# Plot mean 56Ni abundance trajectories of each region with spread
print "plotting mean 56Ni abundances"
plt.figure()
for r in xrange(3):
	mu = [mean_ni56_west, mean_ni56_east, mean_ni56_north][r]
	sig = [sigma_ni56_west, sigma_ni56_east, sigma_ni56_north][r]
	lab = ["West", "East", "North"][r] + " $\\pm$ 1$\\sigma$"
	col = ["red", "yellowgreen", "blue"][r]
	plus = 10.**(np.log10(mu) + np.log10(sig))
	minus = 10.**(np.log10(mu) - np.log10(sig))
	plt.fill_between(time, plus, minus, color=col, alpha=0.2)
	plt.plot(time, mu, color=col, label=lab)
plt.xscale("log")
plt.yscale("log")
plt.xlim(1e0, 1e5)
plt.legend(loc="lower left")
plt.title("${}^{56}$Ni Abundances of 3 Regions of Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("${}^{56}$Ni Mass Fraction")
plt.savefig(PLOT_DIR + "time_vs_ni56_sigma.png", dpi=150)
plt.xlim(1e0, 10.)
#plt.legend(loc="lower right")
plt.savefig(PLOT_DIR + "time_vs_ni56_sigma_zoom.png", dpi=150)
plt.close()

############################################################
# DEFINING AND FITTING MAGKOTSIOS TRAJECTORIES
############################################################

# Progress alert for the user
print "fitting the magkotsios-inspired trajectories"

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
print "plotting densities by region with function fits"
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
	plt.plot(time, dens_exp(time, *exp_best), "k--", label="Exponential" \
		+ "\n($\\rho_0$ = %.2e, $\\tau$ = %.2e)" % (exp_best))
	plt.plot(time, dens_pow(time, *pow_best), "k:", label="Power Law" \
		+ "\n($\\rho_0$ = %.2e)" % (pow_best))
	plt.plot(time, dens_pow2(time, *pow2_best), "k-.", label="Power Law 2" \
		+ "\n($\\rho_0$ = %.2e, $\\alpha$ = %.2e)" % (pow2_best))
	plt.xscale("log")
	plt.yscale("log")
	plt.legend(loc="lower left", fontsize=9)#, title="$\\rho_0$ = %.2e g/cc" % (pow_best))
	plt.title("Densities of %s Region Particles vs. Time" % (lab))
	plt.xlabel("Time (s)")
	plt.ylabel("Density (g/cm$^3$)")
	plotfile = PLOT_DIR + "time_vs_dens_%s.png" % (lab.lower())
	plt.savefig(plotfile, dpi=150)
	plt.close()

# Plots of temperature vs time in each region with fitted trajectories
print "plotting temperatures by region with function fits"
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
	plt.plot(time, temp_exp(time, *exp_best), "k--", label="Exponential" \
		+ "\n($T_0$ = %.2e, $\\tau$ = %.2e)" % (exp_best))
	plt.plot(time, temp_pow(time, *pow_best), "k:", label="Power Law" \
		+ "\n($T_0$ = %.2e)" % (pow_best))
	plt.plot(time, temp_pow2(time, *pow2_best), "k-.", label="Power Law 2" \
		+ "\n($T_0$ = %.2e, $\\alpha$ = %.2e)" % (pow2_best))
	plt.xscale("log")
	plt.yscale("log")
	plt.legend(loc="lower left", fontsize=9)#, title="$T_0$ = %.2e K" % (pow_best))
	plt.title("Temperatures of %s Region Particles vs. Time" % (lab))
	plt.xlabel("Time (s)")
	plt.ylabel("Temperature (K)")
	plotfile = PLOT_DIR + "time_vs_temp_%s.png" % (lab.lower())
	plt.savefig(plotfile, dpi=150)
	plt.close()

# Plots of radial velocity vs time in each region without fitted trajectories
print "plotting radial velocities by region"
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
	plt.xlim(1e0, 1e5)
	#plt.legend()
	plt.title("Velocities of %s Region Particles vs. Time" % (lab))
	plt.xlabel("Time (s)")
	plt.ylabel("Radial Velocity (cm/s)")
	plotfile = PLOT_DIR + "time_vs_rvel_%s.png" % (lab.lower())
	plt.savefig(plotfile, dpi=150)
	plt.close()

