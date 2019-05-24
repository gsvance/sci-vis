#!/usr/bin/env python

# Extract total yields for each of the four regions from the cco2 HDF5 files

# Last modified 5/23/19 by Greg Vance

print "importing modules..."
import numpy as np
import glob
import h5py
print "imports complete!"

# File with particle IDs and region tags produced by make_lists.py
REGIONS_FILE = "regions.out"

# Directory path containing all the cco2 HDF5 data
HDF5_DIR = "/home/gsvance/supernova_data/cco2/cco2.dir/"

# Directory for storing this script's output data
OUT_DIR = "region_total_yields/"

# The number of isotopes in the nuclear network
NET_SIZE = 524

# Read all of the particle IDs and region tags into a dict
print "reading region tags..."
tags = dict()
with open(REGIONS_FILE, "r") as rf:
	rf.readline()  # Read and discard the header line
	for line in rf:
		id, rtag = line.strip().split()
		tags[int(id)] = rtag
print "region tags read!"

# Set up data structures to hold the summed yields from the HDF5 files
nn, nz = None, None
yields = dict()
for reg in xrange(1, 5):
	yields["R%d" % (reg)] = np.zeros(NET_SIZE, np.float64)

# Scour the HDF5 files for all the particle IDs in the region tag file
print "beginning to read hdf5 files..."
for hdf5 in sorted(glob.glob(HDF5_DIR + "outc2?.h5")):
	
	# Open the next HDF5 file
	print " opening %s" % (hdf5)
	f = h5py.File(hdf5, "r")
	
	# Run some checks with the file's nn and nz arrays
	fnn = f["nn"]["data"]
	fnz = f["nz"]["data"]
	if nn is None and nz is None:
		nn, nz = np.copy(fnn), np.copy(fnz)
	else:
		assert all(nn == fnn)
		assert all(nz == fnz)
	
	# Interate over all particle data ("cycles") stored in this file
	cycles = [k for k in f.keys() if (k != "nn" and k != "nz")]
	for cycle in cycles:
		
		# Extract the particle ID from the cycle string and check if we want it
		pid = int(cycle.strip("cycle"))
		if pid in tags:
			reg = tags[pid]
		else:
			continue
		
		# Extract the particle's total mass and fractional isotope masses
		mass = f[cycle].attrs['mass']
		fmass = f[cycle]['fmass']['data']
		assert mass.size == 1
		assert fmass.size == NET_SIZE
		
		# Add the mass of each isotope to the running totals for the correct region
		yields[reg] += fmass * mass
	
	# Close this HDF5 file
	f.close()
print "done reading hdf5 files!"

# Write the collected data out to text files
print "writing out data..."
for reg in xrange(1, 5):
	t = "R%d" % (reg)
	outfile = OUT_DIR + "cco2_reg%d_yields.txt" % (reg)
	with open(outfile, "w") as of:
		of.write("# nn nz mass")
		for inn, inz, imass in zip(nn, nz, yields[t]):
			of.write("\n%d %d %.6e" % (inn, inz, imass))
print "data write out complete!"

