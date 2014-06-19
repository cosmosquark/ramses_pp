'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''
import yt
from .. import Profile
import numpy as np
from yt.analysis_modules.halo_profiler.api import *

verbose = False

def load(snapshot):
	return YTProfile(snapshot)

class YTProfile(Profile.Profile):
	def __init__(self, snapshot):
		Profile.Profile.__init__(self, snapshot, 'yt')
		'''
		TODO Verify snapshot is a yt snapshot
		'''
		self._snapshot = snapshot

	#Implement abstract methods from Profile.py

	def HaloProfiler(self, halos=None):
		#Pretty much returns the usual yt halo profiler,
		#but makes sure halo catalogue exists and points
		#the profiler at existing catalogue
		snapshot = self._snapshot

		if halos == None:
			halos = snapshot.Halos()

		if (verbose): print 'Found %d halos'%len(halos)
		hop_fname = '../hop_dir/out_%05d_hop.out'%snapshot.OutputNumber()

		if (verbose): print 'Profiling halos at %s'%hop_fname

		hp = HaloProfiler(snapshot.raw_snapshot(), halo_list_file=hop_fname)
		return hp

	def fgas_halomass(self, halos=None):
		snapshot = self._snapshot

		if (halos == None):
			halos = snapshot.halos()

		if (verbose): print 'Found %d halos'%len(halos)

		#Use lists as some halos won't be large enough to define a sphere
		fgas = []
		mhalo = []

		#Loop through the halos
		i = 0
		for halo in halos:
			try:
				if (verbose): print 'Processing halo %d'%i
				#Get a sphere defining the halo
				sphere = halo.get_sphere()

				gas_mass, particle_mass = sphere.quantities["TotalQuantity"](
					["CellMassMsun", "ParticleMassMsun"])
				total_mass = particle_mass + gas_mass

				if (verbose): print 'gass mass = %e, total mass = %e, fgas = %f'%(gas_mass, total_mass, gas_mass/total_mass)

				fgas.append(gas_mass/total_mass)
				mhalo.append(total_mass)

			except yt.utilities.exceptions.YTSphereTooSmall as e:
				print 'Sphere too small: ', e
			except ValueError as ve:
				print 'Unexpected ValueError for halo %d: '%i, ve
			i = i + 1

			if (i % 100) == 0: print 'Processed %d halos...'%i

		return np.array(fgas), np.array(mhalo)
