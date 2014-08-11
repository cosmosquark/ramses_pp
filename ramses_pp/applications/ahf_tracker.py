#!/usr/bin/python
from ramses_pp.modules import Simulation
import numpy as np
from ramses_pp import config
import sys, os
import json, glob, uuid

sim = Simulation.load("grid_n08_halo_1")
snapno = sim.num_snapshots()
ahf_files = []
mfiles = []
mtree_files = np.zeros(snapno-1, dtype=str)
for isnap in range(snapno,0,-1):
	tipsy_dir = str("%s/output_%05d/output_%05d_tipsy/" % (sim.path(), isnap, isnap))
	#print tipsy_dir
#	print glob.glob("%s*.AHF_particles*" %tipsy_dir)
	ahf_files.append(glob.glob('%s*.AHF_halos'%tipsy_dir))

print ahf_files

for i in range (0, len(ahf_files)):
	word = ahf_files[i]
	word[0] = word[0][:-6]
	ahf_files[i] = word


#for m in range(snapno,1,-1):
#	mfiles.append("output_%05d_ahf_mtree"%m)

f = open("applications/ahf_track_infile", 'w')

## write the track file

for i in range(0,snapno):
#	if i == (snapno - 1):
	f.write(str(ahf_files[i][0]) + "\n")
#	else:
	#	f.write(str(ahf_files[i]) + "," + str(mfiles[i]) + "\n")
f.close()

# write the redshift file

red = sim.avail_redshifts()
if red[(sim.num_snapshots()-1)] < 0.00:
	red[(sim.num_snapshots()-1)] = 0

z = open("applications/zfile",'w')
for i in range(snapno-1,-1,-1):
	z.write(str(red[i]) + "\n")

z.close()

