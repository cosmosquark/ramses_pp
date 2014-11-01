'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import pymses
from .. import Snapshot
import sys, os

from pymses import RamsesOutput
from pymses.sources.ramses import output
RamsesOutput.amr_field_descrs_by_file = {"3D": {
	"hydro" : [ output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]), output.Scalar("P", 4), output.Scalar("xHII", 5), output.Scalar("xHeII", 6), output.Scalar("xHeIII", 7) ],
	"grav"  : [ output.Vector("g", [0, 1, 2]) ]
	} }

def load(folder, ioutput):
	return PymsesSnapshot(folder, ioutput)

class PymsesSnapshot(Snapshot.Snapshot):
	def __init__(self, folder, ioutput):
		Snapshot.Snapshot.__init__(self, "Pymses")
		'''
		Load the snapshot using pymses.RamsesOutput
		TODO Verify snapshot exists
		'''
		self._path = folder
		self._ioutput = ioutput
		self._snappath = os.path.join('%s/output_%05d/'%(self._path, ioutput))
		self._snapshot = pymses.RamsesOutput(folder, ioutput)

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

	def amr_source(self, fields):
		if isinstance(fields, list):
			return self._snapshot.amr_source(fields)
		else:
			return self._snapshot.amr_source([fields])

	def particle_source(self, fields):
		if isinstance(fields, list):
			return self._snapshot.particle_source(fields)
		else:
			return self._snapshot.particle_source([fields])

	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return self._snapshot

	def ncpu(self):
		'''
		Return the number of CPUs used
		'''
		ro = self.raw_snapshot()
		return ro.info['ncpu']

	def current_redshift(self):
		'''
		Return the current redshift
		'''
		ro = self.raw_snapshot()
		aexp = ro.info['aexp']
		return (1/aexp) - 1

	def cosmology(self):
		'''
		Return an object with cosmological parameters
		'''
		ro = self.raw_snapshot()
		info = ro.info
		omega_m_0 = info['omega_m']
		omega_l_0 = info['omega_l']
		omega_k_0 = info['omega_k']
		omega_b_0 = info['omega_b']
		aexp = info['aexp']
		h = info['H0']/100

		cosmology = {
			'omega_m_0':omega_m_0,
			'omega_l_0':omega_l_0,
			'omega_k_0':omega_k_0,
			'omega_b_0':omega_b_0,
			'h':h,
			'aexp':aexp
		}
		return cosmology

	def info(self):
		'''
		Return info object
		'''
		return self.raw_snapshot().info

class Type:
	AMR = 1
	PART = 2
