#!/usr/bin/python
from ramses_pp.modules import Simulation
import numpy as np
from ramses_pp import config
import sys, os
import json, glob, uuid
import re

sim = Simulation.load("grid_n08_halo_1")
snapno = sim.num_snapshots()
ahf_files = []
red = []
mfiles = []
mtree_files = np.zeros(snapno-1, dtype=str)
for isnap in range(snapno,0,-1):
	tipsy_dir = str("%s/output_%05d/output_%05d_tipsy/" % (sim.path(), isnap, isnap))
	#print tipsy_dir
#	print glob.glob("%s*.AHF_particles*" %tipsy_dir)
	ahf_files.append(glob.glob('%s*.AHF_particles'%tipsy_dir)[0])
	red.append( re.split(r'\.(?!\d)', glob.glob('%s*.AHF_particles'%tipsy_dir)[0])[2] ) # this appends the string value of Z to the list by using regex to delimit the filename by full stops and grabbing the redshift string from that

print ahf_files

#red = sim.avail_redshifts()
#if red[(sim.num_snapshots()-1)] < 0.00:
#        red[(sim.num_snapshots()-1)] = 0

# read z safely from the files




for m in range(snapno,1,-1):
#	mfiles.append("output_%05d_ahf_mtree"%m)
	mfiles.append("%s/output_%05d/output_%05d_tipsy/output_%05d_fullbox.tipsy.%s.AHF" % (sim.path(), m, m, m,str(red[snapno-m])))

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

