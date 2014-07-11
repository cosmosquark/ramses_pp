import ramses_pp
from ramses_pp.applications import halomaker as hm
from ramses_pp.modules import Simulation
sim = Simulation.load('selene')
hm.setup(sim)
