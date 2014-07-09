import loader
from ramses_pp.modules import Simulation, Snapshot

def setup(Simulation):
	test = hasattr(Simulation, 'Simulation')
	print test
	return
