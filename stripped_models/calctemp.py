# Function for calculating the temperature of a zone of material in the models
# This code relies on knowing the density, internal energy, and composition
# Based on code repurposed from Carola's snsphtools/addabundance2v3.c
# None of the cooling table/ion fraction stuff from there is included here
# We're assuming everything is ionized, which is true above a few MK

# Last modified 1/24/20 by Greg Vance

import numpy as np

# Relevant physical constants from snsphtools/consts.h
A_RAD = 7.565700e-15 # erg cm^-3 K^-4
K_BOLTZ = 1.3806503e-16 # erg K^-1
N_AVOG = 6.02214179e23

# This next number (AMU) is the *approximate* mass of a single 1H atom
# It's also *approximately* 1/12 the mass of a single 12C atom
# For any isotope, (Z+N)*AMU is *approximately* the mass of a single atom
# The reason is something along the lines of "this is how moles work"
AMU = 1.0 / N_AVOG # grams

# Most important function in this file... see line 499 of addabundance2v3.c
# This does the inverse of that block of code, getting T from u instead
def calctemp(u_spec, rho, nz_arr, nn_arr, abun_arr):
	
	# Calculate the temperature of the zone of material, all inputs in CGS
	# u_spec is a scalar for the specific internal energy
	# rho is a scalar for the mass density
	# nz_arr is an array of ints which are the proton numbers for all isotopes
	# nn_arr is an array of ints which are the neutron numbers for all isotopes
	# abun_arr is an array of floats which are the isotope mass fractions X
	# These three arrays should all be 1D with the same length
	
	global A_RAD, K_BOLTZ, N_AVOG, AMU
	
	# Sanity-check the input arrays
	N_iso = nz_arr.size
	assert nn_arr.size == N_iso
	assert abun_arr.size == N_iso
	
	u_tot = u_spec * rho  # convert out of specific internal energy
	n_tot = 0.0  # total number density n
	
	# Add number densities for each isotope nucleus
	for i in xrange(N_iso):
		n_tot += (rho * abun_arr[i]) / (float(nz_arr[i] + nn_arr[i]) * AMU)
	
	# Add in the electron number density
	n_e = find_ne_ionized(rho, nz_arr, nn_arr, abun_arr)
	n_tot += n_e
	
	# Now we have a polynomial in T that needs to be solved for the temp
	# u_tot = (3/2) n_tot K_BOLTZ temp + A_RAD temp**4
	# We need to find roots of a simple quartic polynomial to do this
	# There will always be four roots (accounting for multiplicity)
	# There should always be exactly one root that is both real and positive
	# Negative and complex values for temperature are obviously wrong
	poly_coeff = np.array([A_RAD, 0.0, 0.0, 1.5 * n_tot * K_BOLTZ, -u_tot])
	roots = np.roots(poly_coeff)
	good = np.logical_and(np.isreal(roots), np.real(roots) > 0.0)
	temp = roots[good]
	assert temp.size == 1, "more than one temp root available!!!"
	
	return np.real(temp[0])

# Electron number density subroutine for use by calctemp... see line 735 of
# addabundance2v3.c... this is much simplified from that original function
# The simplifying assumption is that everything is ionized, so we always have
# the ion fraction equation to 1.0 for every isotope
# Patrick says this will always be true for temps > a few MK
def find_ne_ionized(rho, nz_arr, nn_arr, abun_arr):
	
	# Calculate the electron number density of fully ionized material
	# rho is a scalar for the mass density
	# nz_arr is an array of ints which are the proton numbers for all isotopes
	# nn_arr is an array of ints which are the neutron numbers for all isotopes
	# abun_arr is an array of floats which are the isotope mass fractions X
	# These three arrays should all be 1D with the same length
	
	global A_RAD, K_BOLTZ, N_AVOG, AMU
	
	# Sanity check the input arrays
	N_iso = nz_arr.size
	assert nn_arr.size == N_iso
	assert abun_arr.size == N_iso
	
	# We need to know how many different elements we have in the abundances
	# This is going to intentionally exclude free neutrons
	# The answer should just be max(nz_arr) assuming the network doesn't skip
	# any elements along the way
	N_el = np.amax(nz_arr)
	
	# Set up the all-important X_el array, which stores element fraction
	# number densities indexed by Z-1 (so H is X_el[0], He is X_el[1], ...)
	X_el = np.zeros(N_el, np.float)
	for i in xrange(N_iso):
		if nz_arr[i] > 0:
			X_el[ nz_arr[i]-1 ] += (rho * abun_arr[i]) \
				/ (float(nz_arr[i] + nn_arr[i]) * AMU)
	
	# Set up the electron fraction at zero and then accumulate each element
	n_e = 0.0
	for n_el in xrange(N_el):
		nz = n_el + 1
		
		# Don't bother looping over all levels of ionization for each element
		# Everything is assumed fully ionized, so we only consider that
		# In the original code, this is the m = N (or m = nz) case
		fracn = 1.0  # assume ion fraction is 1.0, always fully ionized
		n_e += X_el[n_el] * float(nz) * fracn
	
	return n_e

