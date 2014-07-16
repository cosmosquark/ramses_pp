'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import pynbody
import Pynbody as PynbodySnapshot
from .. import Profile
import numpy as np
import sys, os

verbose = False

def load(snapshot):
	return PynbodyProfile(snapshot)

class PynbodyProfile(Profile.Profile):
	def __init__(self, snapshot):
		Profile.Profile.__init__(self, snapshot, "Pynbody")
		'''
		TODO Verify snapshot exists
		'''
		self._snapshot = snapshot

	#Implement abstract methods from Profile.py

	def profile_halo(self, halo, species, min_r=0, max_r=50):
		s = self._snapshot.raw_snapshot()
		#Centre on the disk
		pynbody.analysis.angmom.faceon(halo)

		#Convert to physical units
		s.physical_units()

		#Create profile object for the given species
		if (species == PynbodySnapshot.Species.GAS):
			if(verbose): print 'Profiling halo gas...'
			prof = pynbody.analysis.profile.Profile(halo.g, min=min_r, max=max_r)
		elif (species == PynbodySnapshot.Species.STAR):
			if(verbose): print 'Profiling halo stars...'
			prof = pynbody.analysis.profile.Profile(halo.s, min=min_r, max=max_r)
		elif (species == PynbodySnapshot.Species.DM):
			if(verbose): print 'Profiling halo DM...'
			prof = pynbody.analysis.profile.Profile(halo.d, min=min_r, max=max_r)
		else:
			raise Exception("Unknown species %s"%species)

		return prof

	#TODO Filter by property e.g temperature
	def fgas_halomass(self, halos=None, baryon_fraction=False):
		cosmo = self._snapshot.cosmology()
		h = cosmo['h']

		if (halos == None):
			halos = self._snapshot.halos()

		fgas = [] # len(halos) - 1 ?
		mhalo = []

		print 'Found %d halos'%len(halos)

		#Loop over each halo and determine the gas fraction
		j = 0
		try:
			for halo in halos:
				halo.physical_units()
				print 'Processing halo %d'%j
				gas_mass = np.sum(halo.g['mass']*h)
				stellar_mass = np.sum(halo.s['mass']*h)
				total_mass = np.sum(halo['mass']*h)
				if(verbose): print 'mg = %e, mtot = %e, fg = %d'%(gas_mass, total_mass, gas_mass/total_mass)

				#If baryon_fraction is true, calculate fraction of baryons in halo instead of just gas
				mhalo.append(total_mass)
				if baryon_fraction:
					fgas.append((gas_mass + stellar_mass) / total_mass)
				else:
					fgas.append(gas_mass / total_mass)
				j = j + 1
				#if j == 50: break
		except KeyError:
			print 'KeyError on halo %d'%j

		return np.array(fgas), np.array(mhalo)
