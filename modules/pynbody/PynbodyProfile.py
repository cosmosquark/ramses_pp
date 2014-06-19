'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import pynbody
from ...pynbody import Pynbody as PynbodySnapshot
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
		s = self._snapshot.raw_snapshot()

		s.physical_units()

		if (halos == None):
			halos = s.halos()

		fgas = np.zeros(len(halos)) # len(halos) - 1 ?
		mhalo = np.zeros(len(halos))

		#Loop over each halo and determine the gas fraction
		i = 0
		for halo in halos:
			print 'Processing halo %d'%i
			gas_mass = np.sum(halo.g['mass'])
			stellar_mass = np.sum(halo.s['mass'])
			total_mass = np.sum(halo['mass'])
			if(verbose): print 'mg = %e, mtot = %e, fg = %d'%(gas_mass, total_mass, gas_mass/total_mass)

			#If baryon_fraction is true, calculate fraction of baryons in halo instead of just gas
			mhalo[i] = total_mass
			if baryon_fraction:
				fgas[i] = (gas_mass + stellar_mass) / total_mass
			else:
				fgas[i] = gas_mass / total_mass
			i = i + 1
			#if i == 50: break

		return fgas, mhalo
