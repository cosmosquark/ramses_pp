from ramses_pp.modules import Simulation
import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

	sim = Simulation.load(sys.argv[1])
	snap = sim.snapshot(int(sys.argv[2]), module='yt')
	ds = snap.raw_snapshot()
	halosahf = snap.halos(finder="AHF")
	haloshm = snap.halos(finder="halomaker_simple")
	mbinmps, mhist, mbinsize = halosahf.mass_function(units='Msun/h', nbins=100)
	mbinmps_2, mhist_2, mbinsize_2 = haloshm.mass_function(units='Msun/h', nbins=100)
	snap_pynbody = sim.snapshot(sim.num_snapshots(), module='pynbody')
	M, sigma, N = snap_pynbody.analytical_mass_function()

	boxsize = ds.arr(ds.parameters['unit_l'], 'cm').in_units('Mpccm/h')
	print 'boxsize= ', boxsize

	#Plot the mass function	
	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	ax.semilogy(np.log10(M),N,label="Reed et al 2007")
	ax.semilogy(mbinmps,mhist/(boxsize**3)/mbinsize,'o',label="AHF")
	ax.semilogy(mbinmps_2,mhist_2/(boxsize**3)/mbinsize_2,'o',label="AdaptaHop")
	plt.xlabel(r'log$_{10}$(M) [M$_{\odot}$/h]')
	plt.ylabel('dN / dlog$_{10}$(M) [Mpc$^{-3}$ h$^{3}$]')
	plt.legend()
	fig.savefig('hmf_test_' + str(sys.argv[1]) + "_" + str(sys.argv[2]) + '_.png')

