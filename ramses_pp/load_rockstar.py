#!/usr/bin/env python
from ramses_pp.modules import Simulation
import sys

sim = Simulation.load(sys.argv[1])
snap = sim.snapshot(int(sys.argv[2]))

halos = snap.halos()

for halo in halos:
	print halo.sphere()
