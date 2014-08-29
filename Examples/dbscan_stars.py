from ramses_pp.modules import Simulation
import numpy as np

sim = Simulation.load('4Mpc_aton_256')
snap = sim.snapshot(214, module='yt')
halos = snap.halos()
h1 = halos[1]

db, X = h1.dbscan()
labels = db.labels_
core_samples = db.core_sample_indices_

#Find stars that are in groups (omit noise)
for k in set(labels):
	if k == -1: continue # noise

	class_members = [index[0] for index in np.argwhere(labels == k)]
	core_samples = [index for index in core_samples if labels[index] == k]

	#Grab the position of each star in the group k and check if they are a core sample\
	group = []
	for index in class_members:
		pos = X[index]
		is_core_sample = (index in core_samples)

		#Do some science...
		group.append(pos)

	print group