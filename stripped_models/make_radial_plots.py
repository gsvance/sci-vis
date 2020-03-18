#!/usr/bin/env python

# Make comparative radial plots of the 1D input files
# Also work out and plot differences for the paper

# Last modified 18 Mar 2020 by Greg Vance

import numpy as np
#from scipy import integrate
import matplotlib.pyplot as plt

from calctemp import calctemp

dtype = {"names": ["radius", "density", "internal energy", "pressure",
	"velocity"], "formats": 5 * [np.float]}

c16r4 = np.loadtxt("c16r4-inputmodel.dat", dtype)
run3g = np.loadtxt("run3g-inputmodel.dat", dtype)

# The c16r4 file is in code units, so convert those numbers to cgs
# The run3g file is already in cgs units, so don't worry about those numbers

solmass = 1.9889e33 # g

codemass = 1e-6 * solmass # g
codelength = 6.955e10 # cm
codetime = 100.0 # s
codetemp = 1.0 # K

codeforce = codemass * codelength / codetime ** 2 # dyn
codeenergy = codeforce * codelength # erg

codedensity = codemass / codelength ** 3 # g/cm^3
codeintenergy = codeenergy / codemass # erg/g
codepressure = codeforce / codelength ** 2 # dyn/cm^2
codevelocity = codelength / codetime # cm/s

c16r4["radius"] *= codelength
c16r4["density"] *= codedensity
c16r4["internal energy"] *= codeintenergy
c16r4["pressure"] *= codepressure
c16r4["velocity"] *= codevelocity

# Work out the enclosed mass coordinate the hard (but correct) way

c16r4_massenc = np.empty(c16r4.size, np.float)
c16r4_massenc[0] = 1.35 * solmass  # mass of compact object cut out of center
for index in xrange(1, c16r4.size):
	r_outer = c16r4["radius"][index]
	r_inner = c16r4["radius"][index-1]
	rho_outer = c16r4["density"][index]
	rho_inner = c16r4["density"][index-1]
	shell_volume = (4./3.) * np.pi * (r_outer**3 - r_inner**3)
	rho_average = np.sqrt(rho_outer * rho_inner)  # log-space average
	shell_mass = shell_volume * rho_average
	c16r4_massenc[index] = c16r4_massenc[index-1] + shell_mass
c16r4_massenc /= solmass  # store enclosed mass in units of Msun

run3g_massenc = np.empty(run3g.size, np.float)
run3g_massenc[0] = 1.35 * solmass  # mass of compact object cut out of center
for index in xrange(1, run3g.size):
	r_outer = run3g["radius"][index]
	r_inner = run3g["radius"][index-1]
	rho_outer = run3g["density"][index]
	rho_inner = run3g["density"][index-1]
	shell_volume = (4./3.) * np.pi * (r_outer**3 - r_inner**3)
	rho_average = np.sqrt(rho_outer * rho_inner)  # log-space average
	shell_mass = shell_volume * rho_average
	run3g_massenc[index] = run3g_massenc[index-1] + shell_mass
run3g_massenc /= solmass  # store enclosed masses in usints of Msun

# Plot the enclosed mass as a function of radius for both models
plt.figure()
plt.plot(run3g["radius"], run3g_massenc, color="red", linestyle="solid",
	label="cco2 model (run3g)")
plt.plot(c16r4["radius"], c16r4_massenc, color="blue", linestyle="dotted",
	label="stripped model (c16r4)")
plt.xscale("log")
plt.yscale("linear")
plt.xlabel("Radius (cm)")
plt.ylabel("Enclosed mass (M$_\\odot$)")
plt.legend(loc="lower right")
plt.savefig("plots/radius_vs_massenc.png", dpi=180)
plt.close()

# Don't do this, it's wrong  |
#                            v
## Approximate the enclosed mass by trapezoidal rule using rho and r arrays
#c16r4_dm = 4. * np.pi * c16r4["radius"] ** 2 * c16r4["density"]
#c16r4_massenc = integrate.cumtrapz(c16r4_dm, c16r4["radius"], initial=0.)
#c16r4_massenc /= 1e6 * codemass # enclosed mass in units of Msun
#run3g_dm = 4. * np.pi * run3g["radius"] ** 2 * run3g["density"]
#run3g_massenc = integrate.cumtrapz(run3g_dm, run3g["radius"], initial=0.)
#run3g_massenc /= 1e6 * codemass # enclosed mass in units of Msun

# Read in the composition data for c16r4 and sanity-check it
with open("c16r4_abun.dat", "r") as f:
	c16r4_header = [l.strip() for l in list(f)[:5]]
c16r4_nnet = int(c16r4_header[0])
c16r4_nzone = int(c16r4_header[1]) - 1 # forget outer boundary condition
c16r4_nz = np.array([int(nz) for nz in c16r4_header[2].split()])
assert c16r4_nz.size == c16r4_nnet
c16r4_nn = np.array([int(nn) for nn in c16r4_header[3].split()])
assert c16r4_nn.size == c16r4_nnet
c16r4_iname = np.array([s for s in c16r4_header[4].split()])
assert c16r4_iname.size == c16r4_nnet + 1
assert c16r4_iname[0] == "r"
c16r4_iname = c16r4_iname[1:] # cut off the first column, which is radius
c16r4_irad = np.loadtxt("c16r4_abun.dat", usecols=[0], skiprows=5)
assert c16r4.size == c16r4_nzone
c16r4_abun = np.loadtxt("c16r4_abun.dat", usecols=range(1, c16r4_nnet + 1),
	skiprows=5)
assert c16r4_abun.shape == (c16r4_nzone + 1, c16r4_nnet)
c16r4_abun = c16r4_abun[:-1, :] # drop the outer boundary condition

# Read in the composition data for run3g and sanity-check it
with open("run3g_abun.dat", "r") as f:
	run3g_header = [l.strip() for l in list(f)[:5]]
run3g_nnet = int(run3g_header[0])
run3g_nzone = int(run3g_header[1])
run3g_nz = np.array([int(nz) for nz in run3g_header[2].split()])
assert run3g_nz.size == run3g_nnet
run3g_nn = np.array([int(nn) for nn in run3g_header[3].split()])
assert run3g_nn.size == run3g_nnet
run3g_iname = np.array([s for s in run3g_header[4].split()])
assert run3g_iname.size == run3g_nnet
run3g_abun = np.loadtxt("run3g_abun.dat", skiprows=5)
assert run3g_abun.shape == (run3g_nzone, run3g_nnet)

# Check the composition data against the other data
assert c16r4_nzone == c16r4.size
assert run3g_nzone == run3g.size

# Calculate all the zone temperatures for c16r4
c16r4_temp = np.empty(c16r4_nzone, np.float)
for k in xrange(c16r4_nzone):
	c16r4_temp[k] = calctemp(c16r4["internal energy"][k], c16r4["density"][k],
		c16r4_nz, c16r4_nn, c16r4_abun[k, :])

# Calculate all the zone temperatures for run3g
run3g_temp = np.empty(run3g_nzone, np.float)
for k in xrange(run3g_nzone):
	run3g_temp[k] = calctemp(run3g["internal energy"][k], run3g["density"][k],
		run3g_nz, run3g_nn, run3g_abun[k, :])

# Make graphs of the various interesting quantities

c16r4_ydata = [c16r4["density"], c16r4["internal energy"], c16r4["pressure"],
	c16r4["velocity"], c16r4_temp]
run3g_ydata = [run3g["density"], run3g["internal energy"], run3g["pressure"],
	run3g["velocity"], run3g_temp]
ylabels = ["Density", "Internal energy", "Pressure", "Radial velocity",
	"Temperature"]
yunits = [" (g cm$^{-3}$)", " (erg g$^{-1}$)", " (dyn cm$^{-2}$)",
	" (cm s$^{-1}$)", " (K)"]
short = ["rho", "u", "pr", "vr", "T"]

for i in xrange(len(short)):
	
	# Print first checkpoint for the loop
	print short[i], 1
	
	# Plot ydata against radius for both models
	plt.figure()
	plt.plot(run3g["radius"], run3g_ydata[i], color="red", linestyle="solid",
		label="cco2 model (run3g)")
	plt.plot(c16r4["radius"], c16r4_ydata[i], color="blue", linestyle="dotted",
		label="stripped model (c16r4)")
	plt.xscale("log")
	if short[i] != "vr":
		plt.yscale("log")
	else:
		plt.yscale("linear")
	plt.xlabel("Radius (cm)")
	plt.ylabel(ylabels[i] + yunits[i])
	plt.legend(loc="upper right")
	plt.savefig("plots/radius_%s_profiles.png" % (short[i]), dpi=180)
	plt.close()
	
	# Plot ydata against enclosed mass coordinate for both models
	plt.figure()
	plt.plot(run3g_massenc, run3g_ydata[i], color="red", linestyle="solid",
		label="cco2 model (run3g)")
	plt.plot(c16r4_massenc, c16r4_ydata[i], color="blue", linestyle="dotted",
		label="stripped model (c16r4)")
	if short[i] != "vr":
		plt.yscale("log")
	else:
		plt.yscale("linear")
	plt.xlabel("Enclosed mass (M$_\\odot$)")
	plt.ylabel(ylabels[i] + yunits[i])
	plt.legend(loc="upper right")
	plt.savefig("plots/massenc_%s_profiles.png" % (short[i]), dpi=180)
	plt.close()
	
	# Print second checkpoint for the loop
	print short[i], 2
	
	# Interpolate each model's ydata values at the other's radii values
	# Then calculate relative ratios with the run3g model as baseline
	
	all_radii = np.union1d(run3g["radius"], c16r4["radius"])
	radius_diffs = list()
	for radius in all_radii:
		
		if radius in run3g["radius"]:
			run3g_yvalue = run3g_ydata[i][run3g["radius"] == radius]
			assert run3g_yvalue.size == 1
			run3g_yvalue = run3g_yvalue[0]
		else:
			# Interpolate ydata in log space except for radial velocity
			# Always take the log of the radii
			if short[i] != "vr":
				run3g_yvalue = 10. ** np.interp(np.log10(radius),
					np.log10(run3g["radius"]), np.log10(run3g_ydata[i]),
					np.nan, np.nan)
			else:
				run3g_yvalue = np.interp(np.log10(radius),
					np.log10(run3g["radius"]), run3g_ydata[i], np.nan, np.nan)
		
		if radius in c16r4["radius"]:
			c16r4_yvalue = c16r4_ydata[i][c16r4["radius"] == radius]
			assert c16r4_yvalue.size == 1
			c16r4_yvalue = c16r4_yvalue[0]
		else:
			# Interpolate ydata in log space except for radial velocity
			# Always take the log of the radii
			if short[i] != "vr":
				c16r4_yvalue = 10. ** np.interp(np.log10(radius),
					np.log10(c16r4["radius"]), np.log10(c16r4_ydata[i]),
					np.nan, np.nan)
			else:
				c16r4_yvalue = np.interp(np.log10(radius),
					np.log10(c16r4["radius"]), c16r4_ydata[i], np.nan, np.nan)
		
		radius_diffs.append(c16r4_yvalue / run3g_yvalue)
	radius_diffs = np.array(radius_diffs)
	
	# Interpolate each model's ydata values at the other's mass coord values
	# Then calculate relative ratios with the run3g model as baseline
	
	all_massenc = np.union1d(run3g_massenc, c16r4_massenc)
	massenc_diffs = list()
	for emass in all_massenc:
		
		if emass in run3g_massenc:
			run3g_yvalue = run3g_ydata[i][run3g_massenc == emass]
			assert run3g_yvalue.size == 1
			run3g_yvalue = run3g_yvalue[0]
		else:
			# Interpolate ydata in log space except for radial velocity
			# Never take the log of the mass coords
			if short[i] != "vr":
				run3g_yvalue = 10. ** np.interp(emass, run3g_massenc,
					np.log10(run3g_ydata[i]), np.nan, np.nan)
			else:
				run3g_yvalue = np.interp(emass, run3g_massenc, run3g_ydata[i],
					np.nan, np.nan)
		
		if emass in c16r4_massenc:
			c16r4_yvalue = c16r4_ydata[i][c16r4_massenc == emass]
			assert c16r4_yvalue.size == 1
			c16r4_yvalue = c16r4_yvalue[0]
		else:
			# Interpolate ydata in log space except for radial velocity
			# Never take the log of the mass coords
			if short[i] != "vr":
				c16r4_yvalue = 10. ** np.interp(emass, c16r4_massenc,
					np.log10(c16r4_ydata[i]), np.nan, np.nan)
			else:
				run3g_yvalue = np.interp(emass, c16r4_massenc, c16r4_ydata[i],
					np.nan, np.nan)
		
		massenc_diffs.append(c16r4_yvalue / run3g_yvalue)
	massenc_diffs = np.array(massenc_diffs)
	
	## Point-by-point difference calculations, which require interpolation
	#r_diffs, m_diffs = list(), list()
	#r_interp = np.interp(run3g["radius"], c16r4["radius"], c16r4_ydata[i],
	#	np.nan, np.nan)
	#m_interp = np.interp(run3g_massenc, c16r4_massenc, c16r4_ydata[i],
	#	np.nan, np.nan)
	#for j in xrange(run3g.size):
	#	r_d = 100. * (r_interp[j] - run3g_ydata[i][j]) / (run3g_ydata[i][j])
	#	m_d = 100. * (m_interp[j] - run3g_ydata[i][j]) / (run3g_ydata[i][j])
	#	if not np.isnan(r_d):
	#		r_diffs.append((run3g["radius"][j], r_d))
	#	if not np.isnan(m_d):
	#		m_diffs.append((run3g_massenc[j], m_d))
	#r_interp = np.interp(c16r4["radius"], run3g["radius"], run3g_ydata[i],
	#	np.nan, np.nan)
	#m_interp = np.interp(c16r4_massenc, run3g_massenc, run3g_ydata[i],
	#	np.nan, np.nan)
	#for j in xrange(c16r4.size):
	#	r_d = 100. * (c16r4_ydata[i][j] - r_interp[j]) / (r_interp[j])
	#	m_d = 100. * (c16r4_ydata[i][j] - m_interp[j]) / (m_interp[j])
	#	if not np.isnan(r_d):
	#		r_diffs.append((c16r4["radius"][j], r_d))
	#	if not np.isnan(m_d):
	#		m_diffs.append((c16r4_massenc[j], m_d))
	#r_diffs.sort(key=lambda x: x[0])
	#m_diffs.sort(key=lambda x: x[0])
	#
	#r_r, r_y = np.array(r_diffs).T
	#m_m, m_y = np.array(m_diffs).T
	
	# Print third checkpoint for the loop
	print short[i], 3
	
	plt.figure()
	plt.plot(all_radii, radius_diffs, "k-")
	plt.xscale("log")
	plt.yscale("log")
	plt.xlabel("Radius (cm)")
	plt.ylabel("%s ratio (stripped model / cco2 model)" % (ylabels[i]))
	plt.savefig("plots/radius_%s_ratio.png" % (short[i]), dpi=180)
	plt.close()
	
	plt.figure()
	plt.plot(all_massenc, massenc_diffs, "k-")
	plt.yscale("log")
	plt.xlabel("Enclosed mass (M$_\\odot$)")
	plt.ylabel("%s ratio (stripped model / cco2 model)" % (ylabels[i]))
	plt.savefig("plots/massenc_%s_ratio.png" % (short[i]), dpi=180)
	plt.close()

