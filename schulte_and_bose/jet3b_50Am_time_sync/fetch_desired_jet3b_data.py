#!/usr/bin/env python
# Fetch the desired jet3b SDF data at tpos ~ 43 hrs as described in README

# Last modified 10/22/19 by Greg Vance

import os
import sdfpy as sp
import numpy as np

# Paths to the 50Am and jet3b SDF directories
SDF_DIR_50AM = "/home/gsvance/supernova_data/50Am/run3g_50Am6/"
SDF_DIR_JET3B = "/home/gsvance/supernova_data/jet3b/jet3b/"

# Code units for SNSPH
SNSPH_MASS = 1e-6 * 1.9889e33 # g
SNSPH_LENGTH = 6.955e10 # cm
SNSPH_TIME = 100.0 # s
SNSPH_DENSITY = SNSPH_MASS / SNSPH_LENGTH**3
SNSPH_VELOCITY = SNSPH_LENGTH / SNSPH_TIME

# Quick and dirty extension splitter function
def get_ext(filename):
	return filename.split('.')[-1]

# Function to return the tpos value from an SDF
def get_tpos(filename):
	sdf = sp.SDFRead(filename)
	tpos = sdf.parameters["tpos"]
	del sdf
	return tpos

# List all the 50Am SDFs and sort them by tpos value
all_files = os.listdir(SDF_DIR_50AM)
all_sdfs = [f for f in all_files if get_ext(f).isdigit()]
sdfs_with_tpos = list()
for filename in all_sdfs:
	fullname = os.path.join(SDF_DIR_50AM, filename)
	tpos = get_tpos(fullname)
	sdfs_with_tpos.append((tpos, fullname))
sdfs_with_tpos.sort(key=lambda x: x[0])

# Pull out the final tpos value from the end of 50Am (~43 hrs)
final_50Am_tpos = sdfs_with_tpos[-1][0]
print "50Am final tpos value:", final_50Am_tpos / 36., "hrs"

# Now list all the jet3b SDFs and sort them by tpos value
del all_files, all_sdfs, sdfs_with_tpos
all_files = os.listdir(SDF_DIR_JET3B)
all_sdfs = [f for f in all_files if get_ext(f).isdigit()]
sdfs_with_tpos = list()
for filename in all_sdfs:
	fullname = os.path.join(SDF_DIR_JET3B, filename)
	tpos = get_tpos(fullname)
	sdfs_with_tpos.append((tpos, fullname))
sdfs_with_tpos.sort(key=lambda x: x[0])

# Which jet3b SDF has the nearest tpos value to the final 50Am tpos value?
bestdiff, bestfile = None, None
for jet3b_tpos, jet3b_sdf in sdfs_with_tpos:
	diff = abs(jet3b_tpos - final_50Am_tpos)
	if bestdiff is None or diff < bestdiff:
		bestdiff = diff
		bestfile = jet3b_sdf
best_jet3b_sdf = bestfile
print "nearest jet3b sdf:", best_jet3b_sdf
best_jet3b_tpos = get_tpos(best_jet3b_sdf)
print "nearest jet3b tpos value:", best_jet3b_tpos / 36., "hrs"

# Identify all desired data columns in the jet3b SDF
sdf = sp.SDFRead(best_jet3b_sdf)
iter = sdf.parameters["iter"]
id = sdf["ident"]
x, y, z = sdf["x"], sdf["y"], sdf["z"]
mass = sdf["mass"]
h = sdf["h"]
rho = sdf["rho"]

# Convert all values from code units to CGS
x2 = x * SNSPH_LENGTH
y2 = y * SNSPH_LENGTH
z2 = z * SNSPH_LENGTH
mass2 = mass * SNSPH_MASS
h2 = h * SNSPH_LENGTH
rho2 = rho * SNSPH_DENSITY

# Organize and write out the data
data = np.column_stack([id, x2, y2, z2, mass2, h2, rho2])
form = ["%d"] + 6 * ["%s"]
head = "id, x, y, z, mass, h, density"
outfile = "jet3b_sdf_%d.dat" % (iter)
np.savetxt(outfile, data, fmt=form, delimiter=", ", header=head, comments="")
del sdf

