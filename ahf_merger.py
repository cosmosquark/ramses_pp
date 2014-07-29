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
for isnap in range(1,snapno+1):
	tipsy_dir = str("%s/output_%05d/output_%05d_tipsy/" % (sim.path(), isnap, isnap))
	#print tipsy_dir
#	print glob.glob("%s*.AHF_particles*" %tipsy_dir)
	ahf_files.append(glob.glob('%s*.AHF_particles'%tipsy_dir)[0])

print ahf_files

for m in range(2,snapno+1):
	mfiles.append("output_%05d_ahf_mtree"%m)

f = open("applications/ahf_mtree_infile", 'w')
f.write(str(snapno))
f.write("\n")
for line in ahf_files:
	f.write(line)
	f.write("\n")
for line in mfiles:
	f.write(line)
	f.write("\n")
f.close()

execommand = "./applications/AHF-MergerTree < applications/ahf_mtree_infile"
os.system(execommand)

