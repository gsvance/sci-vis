#!/usr/bin/env python
# list_sdfs.py

# Save a list file of all the cco2 SDF names in order with tpos values
# The list will be read by get_sdf_data.c when it does its job next

# Last modified 2/25/19 by Greg Vance

import os

# Name of SDF list file to output
SDF_FILE = "sdf_list"

# SDF location on Saguaro for cco2
SDF_DIR = "/home/gsvance/supernova_data/cco2/cco2/"

# Function to extract a tpos value from an SDF file
def get_tpos(sdf_name):
	with open(sdf_name, "rb") as sdf_file:
		for line in sdf_file:
			if line.startswith("float tpos = "):
				return line.strip("float ps=;\n")  # Strip these characters
	msg = "SDF %s missing tpos value" % (sdf_name)
	raise ValueError(msg)

# Form a list of all the cco2 SDF file names
sdf_names = []
for maybe_sdf in os.listdir(SDF_DIR):
	extension = maybe_sdf.split('.')[-1]
	try:
		iter = int(extension)
	except ValueError:
		continue
	sdf_name = SDF_DIR + maybe_sdf
	tpos = get_tpos(sdf_name)
	sdf_names.append((iter, tpos, sdf_name))

# Sort the SDF files list by SNSPH iteration
sdf_names.sort(key=lambda s : s[0])

# Write out the list of SDF files with their time data
with open(SDF_FILE, "w") as sdf_file:
	head = "n_sdf %d\n" % (len(sdf_names))
	head += "iter tpos name\n"
	sdf_file.write(head)
	for iter, tpos, name in sdf_names:
		line = "%d %s %s\n" % (iter, tpos, name)
		sdf_file.write(line)

