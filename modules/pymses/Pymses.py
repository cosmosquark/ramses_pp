'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import pymses
from .. import Snapshot
from ... import config
import sys, os
import numpy as np

from pymses import RamsesOutput
from pymses.sources.ramses import output
RamsesOutput.amr_field_descrs_by_file = {"3D": {
		# default
		"hydro" : [ output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]), output.Scalar("P", 4), output.Scalar("metal", 5)],
		# CH
#		"hydro" : [ output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]), output.Scalar("P", 4), output.Scalar("metal", 5), 
#					 output.Scalar("H", 6), output.Vector("O", 7), output.Scalar("Fe", 8),	
#					 output.Scalar("C", 9), output.Vector("N", 10), output.Scalar("Mg", 11), output.Scalar("Si", 12),
#					 output.Scalar("Si", 13), output.Vector("Delay", 14),
# RT	"hydro" : [ output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]), output.Scalar("P", 4), 
#			output.Scalar("metal", 5), output.Scalar("xHII", 6), output.Scalar("xHeII", 7), 
#			output.Scalar("xHeIII", 7) ],
		"grav"  : [ output.Vector("g", [0, 1, 2]) ],
	} }

from pymses.filters import RegionFilter
from pymses.utils.regions import Sphere
from pymses.filters import CellsToPoints

#RamsesOutput.amr_field_descrs_by_file = {"3D": {
#	"hydro" : [ output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]), output.Scalar("P", 4), output.Scalar("xHII", 5), output.Scalar("J21", 6) ],
#	"grav"  : [ output.Vector("g", [0, 1, 2]) ]
#	} }	

def info_dict(snapshot):
	'''
	Expose info API
	'''
	fname = snapshot.info_path()
	from pymses.sources.ramses import info
	return info.read_ramses_info_file(fname)

def hilbert_dom_decomp(snapshot):
	'''
	Expose Hilbert domain decomposition API
	'''
	attrs = snapshot._attributes
	if attrs.has_key('dom_decomp'): return attrs['dom_decomp']

	from pymses.sources.ramses.hilbert import HilbertDomainDecomp
	info = snapshot._info if snapshot._info is not None else info_dict(snapshot)
	keys = info['dom_decomp_Hilbert_keys']
	dom_decomp = HilbertDomainDecomp(info['ndim'], keys[:-1], keys[1:], (info['levelmin'], info['levelmax']))
	snapshot._attributes['dom_decomp'] = dom_decomp
	return dom_decomp

def load(path, ioutput, **kwargs):
	return PymsesSnapshot(path, ioutput, **kwargs)

class PymsesSnapshot(Snapshot.Snapshot):
	def __init__(self, path, ioutput, **kwargs):
		Snapshot.Snapshot.__init__(self, path, ioutput, "pymses", **kwargs)
		'''
		Load the snapshot using pymses.RamsesOutput
		TODO Verify snapshot exists
		'''
		self._snapshot = pymses.RamsesOutput(path, ioutput)

	#Implement abstract methods from Snapshot.py

	def __getitem__(self, item, flatten=True, cpu=None):
		'''
		Return an object containing this field.
		flatten - Return a numpy array containing all points.
		cpu - Return only the domain belonging to this CPU.
		If flatten == False and cpu=None, return an iterable

		possible items
		
		'''
		#Check if AMR field
		#info = self.info()
		ndim = '%dD'%self.info()['ndim']
		fields = [v.name for v in output.RamsesOutput.amr_field_descrs_by_file[ndim]['hydro']]
		if item in fields:
			#Load and return (will work with multiple fields too!)
			#amr_source = self.amr_source(item)
			amr_source = self.amr_source(item)
			cells = CellsToPoints(amr_source)
			if flatten: 
				return cells.flatten()
			elif cpu is None:
				return cells.iter_dsets()
			else:
				return cells.get_domain_dset(cpu)
		else:
			#Particle field
			#part_source = self.particle_source(item)
			part_source = self.particle_source(item)
			if flatten:
				return part_source.flatten()
			elif cpu is None:
				return part_source.iter_dsets()
			else:
				return part_source.get_domain_dset(cpu)

	def get(self, item, flatten=True, cpu=None):
		return self.__getitem__(item, flatten, cpu)

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

	def sphere(self, halo, source):
		'''
		Return a pymses sphere centred on this halo
		'''
		pos = halo['pos'].in_units('code_length').value
		r = halo['Rvir'].in_units('code_length').value
		region = Sphere(pos, r)
		filt_region = RegionFilter(region, source)

		if filt_region is None:
			raise Exception("Unable to create sphere: pos - ", pos + " r - ", r)

		#return filt_region
		return Region(filt_region)

	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return self._snapshot

 	def filt_cube(self, cube, dset_type, fields):
		'''
		Return flattened dset of points contained by this cube
		'''
		if dset_type == Type.AMR:
			source = self.amr_source(fields)
		elif dset_type == Type.PART:
			source = self.particle_source(fields)
		else:
			raise Exception("No such type: ", dset_type)
		from pymses.filters import CellsToPoints
		from pymses.filters import RegionFilter

		filt_source = RegionFilter(cube, source)
		return CellsToPoints(filt_source).flatten()

	def cubes(self, n):
		'''
		Return a list of pymses cubes for a cartesian
		domain decomp
		n - number of cubes per dimension
		'''
		from pymses.utils.regions import Cube
		cube_pos, dx = self.cube_positions(n)
		cubes = []
		for pos in cube_pos:
			cubes.append(Cube(pos, dx))

		return np.array(cubes)

	def ndim(self):
		'''
		Return number of dimensions
		'''
		return self.info()['ndim']

	def boxlen(self):
		'''
		Return boxlen in code units
		'''
		return self.info()['boxlen']

	def ncpu(self):
		'''
		Return the number of CPUs used
		'''
		ro = self.raw_snapshot()
		return ro.info['ncpu']

	def halos(self, finder=config.default_finder, filename=None):
		'''
		Load a generic halo catalogue - default to rockstar if not overridden
		Override the snapshot method, and force loading using yt for unit coherence
		'''
		simulation = self._attributes['simulation']
		yt_snap = self.swap_modules(simulation, 'yt')
		from ...analysis.halo_analysis import halos
		if finder=='rockstar':
			return halos.RockstarCatalogue(yt_snap)
		elif finder=="AHF":
			return halos.AHFCatalogue(yt_snap, filename=filename)
		else:
			raise Exception("Unimplemented finder: %s"%finder)


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

class Region():
	def __init__(self, source, **kwargs):
		'''
		Load the snapshot using pymses.RamsesOutput
		TODO Verify snapshot exists
		'''
		self._snapshot = source

	def __getitem__(self, item, flatten=True, cpu=None):
		'''
		Return an object containing this field.
		flatten - Return a numpy array containing all points.
		cpu - Return only the domain belonging to this CPU.
		If flatten == False and cpu=None, return an iterable
		'''
		source = self._snapshot.source
		ndim = '3D'
		fields = [v.name for v in output.RamsesOutput.amr_field_descrs_by_file[ndim]['hydro']]
		if item in fields:
			#Load and return (will work with multiple fields too!)
			#amr_source = self.amr_source(item)
			cells = CellsToPoints(source)
			if flatten: 
				return cells.flatten()
			elif cpu is None:
				return cells.iter_dsets()
			else:
				return cells.get_domain_dset(cpu)
		else:
			#Particle field
			#part_source = self.particle_source(item)
			if flatten:
				return source.flatten()
			elif cpu is None:
				return source.iter_dsets()
			else:
				return source.get_domain_dset(cpu)

	def get(self, item, flatten=True, cpu=None):
		return self.__getitem__(item, flatten, cpu)[item]
