'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import pynbody
from .. import Snapshot
import sys, os

def load(folder, ioutput=None):
	return PynbodySnapshot(folder, ioutput)

class PynbodySnapshot(Snapshot.Snapshot):
	def __init__(self, folder, ioutput=None):
		Snapshot.Snapshot.__init__(self, "Pynbody")
		'''
		Load the snapshot using pymses.RamsesOutput
		TODO Verify snapshot exists
		'''
		self._path = folder
		self._ioutput = ioutput
		if ioutput is not None:
			self._snapshot = pynbody.load('%s/output_%05d'%(folder, ioutput))
		else:
			self._snapshot = pynbody.load(folder)

	#Implement abstract methods from Snapshot.py

	def output_number(self):
		'''
		Return the output number for this snapshot
		'''
		return self._ioutput

	def path(self):
		'''
		Return the path to this snapshot
		'''
		return self._path

	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return self._snapshot

	def ncpu(self):
		'''
		Return the number of CPUs used
		'''
		s = self.raw_snapshot()
		return s.ncpu

	def current_redshift(self):
		'''
		Return the current redshift
		'''
		s = self.raw_snapshot()
		return s.properties['z']

	def cosmology(self):
		'''
		Return an object with cosmological parameters
		'''
		s = self.raw_snapshot()
		info = s.properties
		omega_m_0 = info['omegaM0']
		omega_l_0 = info['omegaL0']
		h = info['h']

		cosmology = {
			'omega_m_0':omega_m_0,
			'omega_l_0':omega_l_0,
			'h':h
		}
		return cosmology

	def info(self):
		'''
		Return info object
		'''
		return self.raw_snapshot().properties

	def halos(self):
		return self._snapshot.halos()

class Species:
	GAS = 1
	STAR = 2
	DM = 3