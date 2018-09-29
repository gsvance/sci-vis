#!/usr/bin/env python

## Read the fetched data files and do the actual plotting work
## Make spaghetti plots of the yellow region particle trajectories
## Overlay exponential and power law trjectories a la Magkotsios
## Also produce cleaner plots with the visual spread statistics

# Last modified 9/28/18 by Greg Vance

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop

# Fetched data files for this script to read
TIME_FILE = "time_fetched.dat"
RHO_FILE = "rho_fetched.dat"
TEMP_FILE = "temp_fetched.dat"

# File with an alternate particle IDs list tagging them by region
REGIONS_FILE = "regions.out"

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

# Read the region-tagged ID list from file
region_id = np.loadtxt(REGIONS_FILE, np.int, usecols=[0])
region = np.loadtxt(REGIONS_FILE, "S1", usecols=[1])

# Sanity checks on the lists of IDs
assert len(rho_id) == len(temp_id)
assert len(rho_id) == len(region_id)
n_id = len(rho_id)
assert all(rho_id == temp_id)
assert all(rho_id == region_id)
id = rho_id
assert len(region) == n_id

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

raise SystemExit

# Plot of iteration vs temperature with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(iter, temp[i], "r-", alpha=0.05)
plt.yscale("log")
plt.legend()
plt.title("Temperatures of Yellow Region Particles vs. Iteration")
plt.xlabel("SNSPH Iteration")
plt.ylabel("Temperature (K)")
plt.savefig(PLOT_DIR + "iter_vs_temp_lines.png", dpi=150)
plt.close()

############################################################
# CALCULATING MEAN AND SPREAD STATISTICS
############################################################

# Mean and std. dev. of density in log space at every time
mean_dens = 10.**np.mean(np.log10(dens), axis=0)
assert mean_dens.size == n_time
sigma_dens = 10.**np.std(np.log10(dens), axis=0)
assert sigma_dens.size == n_time

# Mean and std. dev. of temperature in log space at every time
mean_temp = 10.**np.mean(np.log10(temp), axis=0)
assert mean_temp.size == n_time
sigma_temp = 10.**np.std(np.log10(temp), axis=0)
assert sigma_temp.size == n_time

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

# Fixed points density and temperature fits must start at
fixed_rho0 = mean_dens[0]
fixed_T0 = mean_temp[0]

# Fit exponential trajectories using scipy.optimize.curve_fit
# We are fitting only the tau parameter, and it is being fit in log space
def dens_exp_fit(t, tau):
	global fixed_rho0
	return np.log10(dens_exp(t, fixed_rho0, tau))
opt, cov = spop.curve_fit(dens_exp_fit, time, np.log10(mean_dens), [1e3])
dens_exp_best = tuple([fixed_rho0] + list(opt))
print "dens exp best: rho0, tau = %.2e, %.2e" % (dens_exp_best)
def temp_exp_fit(t, tau):
	global fixed_T0
	return np.log10(temp_exp(t, fixed_T0, tau))
opt, cov = spop.curve_fit(temp_exp_fit, time, np.log10(mean_temp), [1e3])
temp_exp_best = tuple([fixed_T0] + list(opt))
print "temp exp best: T0, tau = %.2e, %.2e" % (temp_exp_best)

# "Fit" the power law trajectories with no real fitting parameters
dens_pow_best = tuple([fixed_rho0])
print "dens pow best: rho0 = %.2e" % (dens_pow_best)
temp_pow_best = tuple([fixed_T0])
print "temp pow best: T0 = %.2e" % (temp_pow_best)

# Fit power law trajectories with slopes using scipy.optimize.curve_fit
# We are fitting only the alpha parameter, and it is being fit in log space
def dens_pow2_fit(t, alpha):
	global fixed_rho0
	return np.log10(dens_pow2(t, fixed_rho0, alpha))
opt, cov = spop.curve_fit(dens_pow2_fit, time, np.log10(mean_dens), [1.])
dens_pow2_best = tuple([fixed_rho0] + list(opt))
print "dens pow2 best: rho0, alpha = %.2e, %.2e" % (dens_pow2_best)
def temp_pow2_fit(t, alpha):
	global fixed_T0
	return np.log10(temp_pow2(t, fixed_T0, alpha))
opt, cov = spop.curve_fit(temp_pow2_fit, time, np.log10(mean_temp), [1.])
temp_pow2_best = tuple([fixed_T0] + list(opt))
print "temp pow2 best: T0, alpha = %.2e, %.2e" % (temp_pow2_best)

############################################################
# SPAGHETTI PLOTS FOR DENS AND TEMP VS TIME
############################################################

# Plot of time vs density with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(time, dens[i], "b-", alpha=0.05)
plt.plot(time, dens_exp(time, *dens_exp_best), "k--", label="Exponential")
plt.plot(time, dens_pow(time, *dens_pow_best), "k:", label="Power Law")
plt.plot(time, dens_pow2(time, *dens_pow2_best), "k-.", label="Power Law 2")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="lower left")
plt.title("Densities of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "time_vs_dens_lines.png", dpi=150)
plt.close()

# Plot of time vs temperature with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(time, temp[i], "r-", alpha=0.05)
plt.plot(time, temp_exp(time, *temp_exp_best), "k--", label="Exponential")
plt.plot(time, temp_pow(time, *temp_pow_best), "k:", label="Power Law")
plt.plot(time, temp_pow2(time, *temp_pow2_best), "k-.", label="Power Law 2")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="lower left")
plt.title("Temperatures of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.savefig(PLOT_DIR + "time_vs_temp_lines.png", dpi=150)
plt.close()

############################################################
# CLEANER PLOTS OF DENS AND TEMP VS TIME WITH SIGMAS
############################################################

# Plot of time vs density with mean and sigmas
plt.figure()
for s in xrange(4):
	plus = 10.**(np.log10(mean_dens) + s * np.log10(sigma_dens))
	col = ["black", "darkblue", "blue", "lightblue"][s]
	sty = ["solid", "dashdot", "dashed", "dotted"][s]
	lab = ["Mean", "1$\\sigma$", "2$\\sigma$", "3$\\sigma$"][s]
	plt.plot(time, plus, color=col, linestyle=sty, label=lab)
	if s > 0:
		minus = 10.**(np.log10(mean_dens) - s * np.log10(sigma_dens))
		plt.plot(time, minus, color=col, linestyle=sty)
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="lower left")
plt.title("Densities of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "time_vs_dens_sigma.png", dpi=150)
plt.close()

# Plot of time vs temperature with mean and sigmas
plt.figure()
for s in xrange(4):
	plus = 10.**(np.log10(mean_temp) + s * np.log10(sigma_temp))
	col = ["black", "darkred", "red", "pink"][s]
	sty = ["solid", "dashdot", "dashed", "dotted"][s]
	lab = ["Mean", "1$\\sigma$", "2$\\sigma$", "3$\\sigma$"][s]
	plt.plot(time, plus, color=col, linestyle=sty, label=lab)
	if s > 0:
		minus = 10.**(np.log10(mean_temp) - s * np.log10(sigma_temp))
		plt.plot(time, minus, color=col, linestyle=sty)
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="lower left")
plt.title("Temperatures of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperatures (K)")
plt.savefig(PLOT_DIR + "time_vs_temp_sigma.png", dpi=150)
plt.close()

