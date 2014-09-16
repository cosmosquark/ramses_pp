from ramses_pp.modules import Simulation

sim = Simulation.load("grid_n08_halo_1")
#sim.run_ahf()
#sim.run_ahf_merger()
#sim.run_ahf_tracker(snaptime=125,halos=[1])
sim.run_ahf_tracker() # in this case, run the tracker on all halos
