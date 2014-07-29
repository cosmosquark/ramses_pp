#!/usr/bin/python
from ramses_pp.modules import Simulation
sim = Simulation.load("grid_n08_halo_1")
snapno = sim.num_snapshots()
for i in range(1,snapno+1):
	shot = sim.snapshot(i)
	shot.tipsy(gas=False)
	shot.halos()

