#!/usr/bin/env python

# Read the fetched data files and do the actual plotting work
# Make spaghetti plots of the yellow region particle trajectories
# Overlay exponential and power law trjectories a la Magkotsios
# Also produce cleaner plots with the visual spread statistics

# Last modified 9/21/18 by Greg Vance

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop

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
# CALCULATING MEAN & SPREAD STATS AND FITTING TRAJECTORIES
############################################################

# Mean and std. dev. of density at every time
mean_dens = np.mean(dens, axis=0)
assert mean_dens.size == n_time
sigma_dens = np.std(dens, axis=0)
assert sigma_dens.size == n_time

# Mean and std. dev. of temperature at every time
mean_temp = np.mean(temp, axis=0)
assert mean_temp.size == n_time
sigma_temp = np.std(temp, axis=0)
assert sigma_temp.size == n_time

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

# Fit exponential trajectories using scipy.optimize.curve_fit
opt, cov = spop.curve_fit(dens_exp, time, mean_dens, [1e4, 1e4], sigma_dens)
dens_exp_best = tuple(opt)
print "dens exp best: rho0, tau = %.2e, %.2e" % (dens_exp_best)
opt, cov = spop.curve_fit(temp_exp, time, mean_temp, [1e7, 100.], sigma_temp)
temp_exp_best = tuple(opt)
print "temp exp best: T0, tau = %.2e, %.2e" % (temp_exp_best)

# Fit power law trajectories using scipy.optimize.curve_fit
opt, cov = spop.curve_fit(dens_pow, time, mean_dens, [1e6], sigma_dens)
dens_pow_best = tuple(opt)
print "dens pow best: rho0 = %.2e" % (dens_pow_best)
opt, cov = spop.curve_fit(temp_pow, time, mean_temp, [1e10], sigma_temp)
temp_pow_best = tuple(opt)
print "temp pow best: T0 = %.2e" % (temp_pow_best)

############################################################
# SPAGHETTI PLOTS FOR DENS AND TEMP VS TIME
############################################################

# Plot of time vs density with all lines
plt.figure()
for i in xrange(n_id):
	plt.plot(time, dens[i], "b-", alpha=0.2)
plt.plot(time, dens_exp(time, *dens_exp_best), "k--", label="Exponential")
plt.plot(time, dens_pow(time, *dens_pow_best), "k:", label="Power Law")
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
	plt.plot(time, temp[i], "r-", alpha=0.2)
plt.plot(time, temp_exp(time, *temp_exp_best), "k--", label="Exponential")
plt.plot(time, temp_pow(time, *temp_pow_best), "k:", label="Power Law")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="lower left")
plt.title("Temperatures of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.ylim(1e-5, 1e12)
plt.savefig(PLOT_DIR + "time_vs_temp_lines.png", dpi=150)
plt.close()

############################################################
# CLEANER PLOTS OF DENS AND TEMP VS TIME WITH SIGMAS
############################################################

# Plot of time vs density with mean and sigmas
plt.figure()
for s in xrange(4):
	plus = mean_dens + s * sigma_dens
	col = ["black", "navy", "blue", "cyan"][s]
	sty = ["solid", "dashdot", "dashed", "dotted"][s]
	lab = ["Mean", "1$\\sigma$", "2$\\sigma$", "3$\\sigma$"][s]
	plt.plot(time, plus, color=col, linestyle=sty, label=lab)
	if s > 0:
		minus = mean_dens - s * sigma_dens
		plt.plot(time, minus, color=col, linestyle=sty)
#plt.plot(time, temp_exp(time, *temp_exp_best), "k--", label="Exponential")
#plt.plot(time, temp_pow(time, *temp_pow_best), "k:", label="Power Law")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="lower left")
plt.title("Densities of Yellow Region Particles vs. Time")
plt.xlabel("Time (s)")
plt.ylabel("Density (g/cm$^3$)")
plt.savefig(PLOT_DIR + "time_vs_dens_sigma.png", dpi=150)
plt.close()

