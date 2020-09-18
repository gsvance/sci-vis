#!/usr/bin/env python2
# Fetch thermodynamic histories (temp & rho) for all particles in these sims:
#   50Am, jet3b, cco2
# This task was needed for the referee report on Jack & Maitrayee's paper
# Kind of a messy, hacky implementation here...
# Watch out if you want to run this again in the future...

# THIS DOESN'T WORK AT ALL FOR THE CCO2 TRAJECTORIES--YOU WILL GET AN ERROR
# CCO2 SPLITS AND REPARTITIONS PARTICLES, SO THERE ARE NEW PARTICLES APPEARING
# THE FINAL SDF CONTAINS SOME PARTICLE IDS WHICH ARE NOT IN THE INITIAL SDF
# LUCKILY FOR ME, JACK ONLY NEEDED THE TRAJECTORIES FOR 50AM AND JET3B...

# Last modified 18 Sep 2020 by Greg Vance



print "preparing..."

import glob
import sdfpy
import numpy as np

SNDATA = "/home/gsvance/supernova_data/"

# Lists of things to process
n = 6
sims = ["50Am", "50Am", "jet3b", "jet3b", "cco2", "cco2"]
values = ["temp", "rho", "temp", "rho", "temp", "rho"]
sdfs = [
	glob.glob(SNDATA + "50Am/run3g_50Am6/run3g_50Am6_sph.????"),
	glob.glob(SNDATA + "50Am/run3g_50Am6/run3g_50Am6_sph.????"),
	glob.glob(SNDATA + "jet3b/jet3b/run3g_50Am_jet3b_sph.????"),
	glob.glob(SNDATA + "jet3b/jet3b/run3g_50Am_jet3b_sph.????"),
	glob.glob(SNDATA + "cco2/cco2/r3g_1M_cco2_sph.*00"),
	glob.glob(SNDATA + "cco2/cco2/r3g_1M_cco2_sph.*00")
]

assert len(sims) == n
assert len(values) == n
assert len(sdfs) == n

print "ready to go!\n"



for i in xrange(n):
	
	# Open files, sort by tpos
	readers = [sdfpy.SDFRead(s) for s in sdfs[i]]
	readers.sort(key = lambda read : read.parameters['tpos'])
	
	print sims[i] + " sdfs opened"
	
	# Set up arrays to hold data
	tpos = np.empty(len(readers), "float32")
	pid = np.array(readers[-1]['ident'], "uint32")
	data = np.empty((pid.size, tpos.size), "float32")
	
	print sims[i] + " " + values[i] + " arrays ready"
	
	# Read data into arrays carefully
	# Ignore any PIDs that aren't in the final SDF file
	for j, read in enumerate(readers):
		tpos[j] = read.parameters['tpos']
		final_index = 0
		for k in xrange(read['ident'].size):
			if read['ident'][k] < pid[final_index]:
				continue
			elif read['ident'][k] > pid[final_index]:
				raise Exception("%d %d %d %d %d" % (read.parameters['iter'],
					k, read['ident'][k], final_index, pid[final_index]))
			data[final_index, j] = read[values[i]][k]
			final_index += 1
		assert final_index == pid.size
	
	print sims[i] + " " + values[i] + " data read"
	
	# Convert data to standard units
	tpos *= 100.0 # sec / SNSPH_TIME
	if values[i] == "rho":
		msun = 1.9889e33 # grams
		rsun = 6.955e10 # cm
		data *= 1e-6 * msun / rsun**3 # grams per cm**3 / SNSPH_DENISTY
	elif values[i] == "temp":
		data *= 1.0 # kelvin / SNSPH_TEMP
	
	print sims[i] + " " + values[i] + " data converted"
	
	# Write data out to file
	fname = sims[i] + "_" + values[i] + ".dat"
	with open(fname, 'w') as outf:
		times = [str(time) for time in tpos]
		outf.write("time ")
		outf.write(" ".join(times))
		del times
		for k in xrange(pid.size):
			outf.write("\n" + str(pid[k]) + " ")
			row = [str(data[k, j]) for j in xrange(tpos.size)]
			outf.write(" ".join(row))
			del row
	
	print "wrote file: " + fname
	del fname
	
	# Clear data from RAM
	del tpos, pid, data
	while len(readers) > 0:
		del readers[-1]
	del readers
	
	print



