import ramses_pp
from ramses_pp.modules import Simulation
sim = Simulation.new('selene')
#sim = Simulation.new('grid_n08_source')
print sim.halomaker_info()


print "initial conditions"
print sim.initial_conditions()
print "end"

halomaker = sim.halomaker()
print halomaker.sim_info()

print halomaker.output_dir()

halomaker.run_halomaker()
