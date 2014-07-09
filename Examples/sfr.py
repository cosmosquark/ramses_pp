from ramses_pp.modules import Simulation
from ramses_pp.fortran import friedman as fm
import numpy as np
import matplotlib.pyplot as plt
import sys, pynbody, yt

def find_sfh(ax, tform, mass, halo_mass=None, bins=100):
		#Expects Mass in Msun
		trange = [tform.min(),tform.max()]

		#1e-6 for megayears
		#binnorm = 1e-6*bins / (trange[1] - trange[0])
		binnorm = 1e-9*bins / (trange[1] - trange[0])
		weight = mass * binnorm
		#sfh,sfhbines = np.histogram(tform, weights=weight, bins=bins)

		if halo_mass is None:
			sfhist, thebins, patches = ax.hist(tform, weights=weight, bins=bins,
						histtype='step')
		else:
			sfhist, thebins, patches = ax.hist(tform, weights=weight, bins=bins,
						histtype='step', label=r'M$_{H}$ = %1.2e M$_{\odot}$'%halo_mass)

		#sfhtimes = 0.5*(sfhbines[1:]+sfhbines[:-1])
		#return sfh, sfhtimes
		return sfhist, thebins

def z_ticks(s, mass, age, ticks):
	z_ticks = []
	# Skip last tick to keep it tidy
	for tick in ticks[:-1]:
		idx = (np.abs(age-tick)).argmin()

		z_ticks.append(pynbody.analysis.cosmology.redshift(s, age[idx]))
	print 'z ticks: ', z_ticks[::-1]

	return np.around(z_ticks[::-1], decimals=1)

def find_nearest(array,value):
	idx = (np.abs(array-value)).argmin()
	return array[idx]

def make_twin(s, host, thebins, age):
	# Figure out where to put z ticks
	twin = host.twiny()
	twin.set_xlabel(r'z')

	xlim = host.get_xlim()
	xvalues=thebins 

	ticks = np.linspace(min(xlim), max(xlim), num=8)
	print 'ticks: ', ticks
	redshift_ticks = z_ticks(s, age, xvalues, ticks)

	tick_pos = []
	# Skip last tick to keep it tidy
	for tick in ticks:
		tick_pos.append(find_nearest(age, tick))

	twin.set_xticks(tick_pos)
	twin.set_xticklabels(redshift_ticks)

def main(name, z, rt=True):
	fig, ax0 = plt.subplots(nrows=1, sharex=True)
	plt.xlabel(r'Lookback Time [Gyr]')
	plt.ylabel(r'Star Formation Rate [M$_{\odot}$/yr]')

	sim = Simulation.load(name)
	idx = sim.redshift(z)
	snapshot = sim.snapshot(idx)
	pf = snapshot.raw_snapshot()
	z = pf.current_redshift
	aexp = 1/(1+z)

	cosmology = snapshot.cosmology()
	h0 = cosmology['h']*100

	friedman = snapshot.integrate_friedman(store=True)
	print friedman

	age_simu = friedman['age_simu']
	time_simu = friedman['time_simu']

	halos = snapshot.halos()
	sorted_halos = sorted(halos, key=lambda x: x.total_mass(), reverse=True)

	for i in range(0, 8):

		halo = sorted_halos[i**3]

		try:
			sphere = halo.get_sphere()
			ct = sphere.particles.source['particle_age']
			stars = (ct != 0)

			gas_mass, particle_mass = sphere.quantities["TotalQuantity"](
								["CellMassMsun", "ParticleMassMsun"])

			total_mass = gas_mass+particle_mass

			ct_stars = ct[stars]
			age = snapshot.tform(ct_stars) # Gyr

			s = sim.snapshot(idx, module='pynbody').raw_snapshot()

			mass = sphere.particles.source['ParticleMassMsun'][stars]
			print 'Total stellar mass %e'%np.sum(mass)
			print 'Total halo mass %e'%total_mass
			print 'Min age = %f, Max age = %f, Mean age = %f'%(min(age), max(age), np.mean(age))

			sfh, sfhtimes = find_sfh(ax0, age, mass, total_mass, bins=50)
			print 'Min SFR = %f, Max SFR = %f'%(min(sfh), max(sfh))
			#plt.xlim(0, age_simu)
			if i == 0:
				make_twin(s, ax0, age_simu-sfhtimes, age_simu-age)

			#plt.semilogx(sfhtimes, sfh)
		except yt.utilities.exceptions.YTSphereTooSmall as e:
			print 'Sphere too small: ', e
		except ValueError as ve:
			print 'ValueError: ', ve

	#t = plt.title(r'Star Formation Rate. z=%10.5f, M$_{H}$=%e'%(z, total_mass))
	#t.set_y(1.09)

	box = ax0.get_position()
	ax0.set_position([box.x0, box.y0 + box.height * 0.1,
					 box.width, box.height * 0.9])

	# Put a legend below current axis
	ax0.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
			  fancybox=True, shadow=True, ncol=5, prop={'size':7.5})

	#plt.ylim(0, 0.2)
	#plt.semilogy()
	plt.savefig('%s_sfr_%05d.png'%(name, idx))

if __name__ == "__main__":
	name = sys.argv[1]
	z = float(sys.argv[2])
	rt = 'rt' in name
	main(name, z, rt=rt)
