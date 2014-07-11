from ramses_pp.modules import Simulation
from ramses_pp.modules.pymses import Pymses
from ramses_pp.modules.pymses import PymsesProjection
from pymses.analysis.visualization import *
import matplotlib.pyplot as plt
from pymses.utils import constants as C

import sys, os

sim_name = "grid_n08_source"
sim = Simulation.load(sim_name)
num_snapshots = sim.num_snapshots()


for i in range(0,2000):
	for j in range(1,num_snapshots):
		snapshot = sim.snapshot(j, module='pymses')
		field = "rho"
		source_type = Pymses.Type.AMR
		method='RayTracer'
		cmap='jet'

		info = snapshot.info()
		scale = info["unit_length"].express(C.Mpc)
		#Temp hard-coding
		if field == 'rho':
			#Hydrogens per cubic centimetre
			factor = info["unit_density"].express(C.H_cc)
		else: factor = 1

		#Set our camera
		cam = PymsesProjection.camera(center=[0.5, 0.5, 0.5], region_size=[1., 1.])
		op = PymsesProjection.scalar_operator(field, factor)

		#cam = PymsesProjection.default_camera()

		#Create a PymsesProjection object which points to everything required for projections
		proj = PymsesProjection.load(snapshot, source_type, [field], camera=cam, operator=op)
		if hasattr(proj, method):
			#Call the desired function, if it exists
			func = getattr(proj, method)
			map = func()

		else:
			raise Exception("Method %s not found"%method)		
