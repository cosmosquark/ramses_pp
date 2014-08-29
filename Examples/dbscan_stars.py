from ramses_pp.modules import Simulation
import numpy as np
import sys

ioutput = int(sys.argv[1])
sim = Simulation.load('4Mpc_aton_256')
snap = sim.snapshot(ioutput, module='yt')
halos = snap.halos()
h_num = int(sys.argv[2])
h1 = halos[h_num]

db, X = h1.dbscan()
print X
labels = db.labels_

#import sys
#for label in labels:
#	print label

#sys.exit(1)

core_samples = db.core_sample_indices_

groups = []
#Find stars that are in groups (omit noise)
for k in set(labels):
	if k == -1: continue # noise

	class_members = [index[0] for index in np.argwhere(labels == k)]
	core_samples = [index for index in core_samples if labels[index] == k]

	#Grab the position of each star in the group k and check if they are a core sample\
	for index in class_members:
		pos = X[index]
		is_core_sample = (index in core_samples)

		#Do some science...
		groups.append([pos[0], pos[1], pos[2], k])

np.savetxt('./grps_%05d_%d.dat'%(ioutput, h_num), groups)