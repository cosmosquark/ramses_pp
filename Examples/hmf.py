from ramses_pp.modules import Simulation
import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

	sim = Simulation.load(sys.argv[1])
	snap = sim.snapshot(int(sys.argv[2]), module='yt')
	ds = snap.raw_snapshot()
	halos = snap.halos()
	mbinmps, mhist, mbinsize = halos.mass_function(units='Msun/h', nbins=100)
	snap_pynbody = sim.snapshot(sim.num_snapshots(), module='pynbody')
	M, sigma, N = snap_pynbody.analytical_mass_function(kern='REEDU')

	boxsize = ds.arr(ds.parameters['unit_l'], 'cm').in_units('Mpccm/h')
	print 'boxsize= ', boxsize

	#Plot the mass function	
	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	ax.semilogy(np.log10(M),N,label="Reed et al 2007")
	ax.semilogy(mbinmps,mhist/(boxsize**3)/mbinsize,'o')
	plt.xlabel(r'log$_{10}$(M) [M$_{\odot}$/h]')
	plt.ylabel('dN / dlog$_{10}$(M) [Mpc$^{-3}$ h$^{3}$]')
	fig.savefig('hmf_test.png')

