import loader
from ramses_pp.modules import Simulation, Snapshot
from ramses_pp import config
import os

def setup(Simulation, subvol=False, ncpu=1, ):
	is_simulation = hasattr(Simulation, 'initial_conditions')
	if is_simulation == False:
		print "this is not a simulation object"
		return
	
	halomaker_dir = Simulation.data_dir(self) + "/HaloMaker"

	# make the halomaker data dir if it does not exist
	if not os.path.exists(directory):
		os.makedirs(directory)

class HaloMaker():
	def __init__(self):
		print "loading halomaker location"
		application_dir = config.applications_dir
		halomaker_dir = application_dir + "/HaloMaker/f90
		self._halomaker_exe = halomaker_dir + "/HaloFinder"
		self._halomaker_dir = halomaker_dir
		self._halomaker_pp = 
		print "loading initial variables"
		

	def executable_loc(self)
		return self._executable_loc
