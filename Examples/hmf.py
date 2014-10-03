from ramses_pp.modules import Simulation
import sys
import numpy as np
import matplotlib.pyplot as plt

def mass_function(snap,finder,boxsize,ax):
	halos = snap.halos(finder=finder)
	print 'Got %d halos'%len(halos)
	print 'boxsize= ', boxsize
	mbinmps, mhist, mbinsize = halos.mass_function(units='Msun/h', nbins=120)
	ax.semilogy(mbinmps,mhist/(boxsize**3)/mbinsize,'o',label=finder)

if __name__ == "__main__":

	kern="REEDU"
	if kern == "REEDU":
		kernlabel = "Reed et al 2007"

	sim = Simulation.load(sys.argv[1])
	snap = sim.snapshot(sim.redshift(float(sys.argv[2])), module='yt') # input redshift
	ds = snap.raw_snapshot()
	snap_pynbody = sim.snapshot(sim.num_snapshots(), module='pynbody')
	M, sigma, N = snap_pynbody.analytical_mass_function(kern=kern)

	boxsize = ds.arr(ds.parameters['unit_l'], 'cm').in_units('Mpccm/h')
#	print 'boxsize= ', boxsize

	#Plot the mass function	
	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	#ax.semilogy(np.log10(M),N,label="Reed et al 2007")
	ax.semilogy(np.log10(M*(ds.parameters['H0']/100)), N, label=kernlabel)
	mass_function(snap, 'AHF', boxsize, ax)
	mass_function(snap, 'halomaker_simple', boxsize, ax)
	
	plt.xlabel(r'log$_{10}$(M) [M$_{\odot}$/h]')
	plt.ylabel('dN / dlog$_{10}$(M) [Mpc$^{-3}$ h$^{3}$]')
	plt.legend(loc='lower left')
	plt.title(sys.argv[1])
	
	fig.savefig('hmf_test_' + str(sys.argv[1]) + "_" + str(sys.argv[2]) + '_.png')
