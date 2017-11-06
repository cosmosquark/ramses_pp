'''
Heavily based on halos.py from pynbody, so credit to those guys
Rewritten to allow loading of Rockstar catalogues when using any module

Currently depends on yt-3.0 for unit handling. I'll try and code this
dependency out sometime in the future if needed.

@author dsullivan
'''

super_verbose=False
OFFSET=30

import numpy as np
import sys, os.path, glob
import weakref
from ramses_pp import config as pp_cfg
#from yt.units.yt_array import YTArray
from halos import Halo
#import logging

#Should adopt this throughout
#logger = logging.getLogger('ramses_pp.halo')

class SimpleField():

	def __init__(self, field, unit):
		self.value = field
		self.unit = unit

	def __str__(self):
		return str(np.array([self.value, self.unit]))

class MergerHalo(Halo):
	"""
	Generic class representing a halo.
	TODO - Units registry and return object with units
	"""

	def __init__(self, halo_id, halo_catalogue, *args):
		Halo.__init__(self, halo_id, halo_catalogue, *args)

	def is_subhalo(self, otherhalo):
		raise NotImplementedError("Not setup children array for merger tree")

	def __getitem__(self, item):
		'''
		Load the correct snapshot from the simulation to make sure
		we have the right cosmological params for this redshift
		'''
		unit = self._halo_catalogue.units[item]
		return SimpleField(self.properties[item], unit)

	def get_sphere(self):
		'''
		YT only function. Currently, there is a bug with ds.sphere which ignores / h units.
		So for the time being this function will take care of that
		'''
		centre = self['pos']
		Rvir = self['Rvir']
		sim = self._halo_catalogue._base()
		ioutput = self._halo_catalogue['snap_num']+self._halo_catalogue.get_offset()
		snapshot = sim.snapshot(ioutput, module='yt')
		if not (str(type(snapshot)) == "<class 'ramses_pp.modules.yt.YT.YTSnapshot'>"):
			raise NotImplementedError("sphere only implemented for YT")
		ds = snapshot.raw_snapshot()

		#if(centre.uq == YTArray(1.0, 'Mpc/h') or (Rvir.uq == YTArray(1.0, 'kpc/h'))):
		#	if super_verbose:
		#		print "Warning: YT ds.sphere() but does not recognise /h units. Correcting, but double check!"
		#	return ds.sphere(centre/ds.hubble_constant, Rvir/ds.hubble_constant)

		return ds.sphere(centre, Rvir)

class MergerTree(object):

	'''
	Generic halo catalogue bases on pynbody, but with non-essential functionality
	stripped out
	TODO - Return halos.Halo object
	'''
	def __init__(self):
		self._halos = {}

	def calc_item(self, i):
		if i in self._halos:
			return self._halos[i]
		else:
			h = self._get_halo(i)
			self._halos[i] = h
			return h

	def __len__(self):
		return len(self._halos)

	def __iter__(self):
		return self._halo_generator()

	def __getitem__(self, item):
		if isinstance(item, slice):
			indices = item.indices(len(self._halos))
			[self.calc_item(i) for i in range(*indices)]
			return self._halos[item]
		else:
			return self.calc_item(item)

	def _halo_generator(self):
		i = 0
		while True:
			try:
				yield self[i]
				i += 1
				if i >= len(self._halos):
					break
			except RuntimeError:
				break

	def __contains__(self, haloid):
		return self.contains(haloid)

	@staticmethod
	def can_load(self):
		return False

	#TODO - Move all running of halo cats into here
	@staticmethod
	def can_run(self):
		return False

#Rockstar
class RockstarMergerTree(MergerTree):
	'''
	Class to handle catalogues produced by Rockstar (by Peter Behroozi).
	'''	

	halo_type = np.dtype([('scale','f'),('id',np.int64),('desc_scale','f'),('desc_id',np.int64),
					('num_prog',np.int64),('pid',np.int64),('upid',np.int64),
					('desc_pid',np.int64),('phantom','f'),('sam_mvir','f'),
					('Mvir','f'),('Rvir','f'),('Rs','f'),('v_rms','f'),
					('mmp',np.int64),('scale_of_last_MM','f'),('v_max','f'),
					('pos','f',3),('vel','f',3),('J','f',3),
					('spin','f'),('breadth_first_id',np.int64),('depth_first_id',np.int64),
					('tree_root_id',np.int64),('orig_halo_id',np.int64),
					('snap_num',np.int64),('next_prog_depth_first_id',np.int64),
					('last_prog_depth_first_id',np.int64),('klypin_rs','f'),
					('Mvir_all','f'),('M200b','f'),('M200c','f'),('M500c','f'),
					('M2500c','f'),('Xoff','f'),('Voff','f'),('bullock_spin','f'),
					('b_to_a','f'),('c_to_a','f'),('A','f',3),('T/|U|','f'),
					('M_pe_Behroozi','f'),('M_pe_Diemer','f')
					])

	units = {'scale':'dimensionless',
			'id':'dimensionless',
			'desc_scale':'dimensionless',
			'desc_id':'dimensionless',
			'num_prog':'dimensionless',
			'pid':'dimensionless',
			'upid':'dimensionless',
			'desc_pid':'dimensionless',
			'phantom':'bool',
			'sam_mvir':'Msun / h',
			'Mvir':'Msun / h',
			'Rvir':'kpccm / h',
			'Rs':'kpccm / h',
			'v_rms':'km / s',
			'mmp':'bool',
			'scale_of_last_MM':'dimensionless',
			'pos':'Mpccm / h',
			'vel':'km / s',			
			'J':'(Msun/h) * (Mpc/h) * km/s',
			'spin':'dimensionless',
			'breadth_first_id':'dimensionless',
			'depth_first_id':'dimensionless',
			'tree_root_id':'dimensionless',
			'orig_halo_id':'dimensionless',
			'snap_num':'dimensionless',
			'next_prog_depth_first_id':'dimensionless',
			'last_prog_depth_first_id':'dimensionless',
			'klypin_rs':'kpccm / h',
			'Mvir_all':'Msun / h',
			'M200b':'Msun / h',
			'M200c':'Msun / h',
			'M500c':'Msun / h',
			'M2500c':'Msun / h',
			'Xoff':'kpccm / h',
			'Voff':'km / s',
			'bullock_spin':'J / (2*GMVR)**0.5',
			'b_to_a':'dimensionless',
			'c_to_a':'dimensionless',
			'A':'kpccm / h',
			'T/|U|':'dimensionless',
			'M_pe_Behroozi':'Msun / h',
			'M_pe_Diemer':'Msun / h'}

	def __init__(self, simulation, filename=None, make_grp=None):
		# TODO - Read/store header
		# TODO - Abstract out file reading so a RockstarMergerTree can be a sub-catalogue of itself
		if not self._can_load(simulation):
			raise Exception("Cannot locate/load consistent trees")

		self._base = weakref.ref(simulation)
		#self._base = snap
		MergerTree.__init__(self)
		
		if filename is not None: self._rsFilename = filename
		else:
			fname = '%s/%s/trees/tree_0_0_0.dat'%(simulation.path(), pp_cfg.rockstar_base)
			self._rsFilename = cutgz(glob.glob(fname)[0])
			if True == pp_cfg.verbose:
				print 'Loading from %s'%fname

		#Try and open the files
		try:
			f = open_(self._rsFilename)
		except IOError:
			raise IOError("Consistent tree not found -- check the file name of catalogue data or try specifying a catalogue using the filename keyword")

		#self._head = np.fromstring(f.read(self.head_type.itemsize),
		#	dtype=self.head_type)
		#unused = f.read(256 - self._head.itemsize)

		#self._nhalos = self._head['num_halos'][0]

		if True == pp_cfg.verbose:
			print "RockstarMergerTree: loading consistent trees...",
			sys.stdout.flush()

		self._load_consistent_tree(f, simulation)
		f.close()

	@property
	def base(self):
		return self._base()

	def _load_consistent_tree(self, f, simulation):
		ctree = np.loadtxt(f, dtype=self.halo_type, comments='#', skiprows=46)
		self._tree = ctree
		self._nhalo = len(self._tree)

		for h in xrange(self._nhalo):
			self._halos[h] = MergerHalo(self._tree[h]['id'], self)
			self._halos[h].properties = self._tree[h]
		print 'done!'

	def load_copy(self, i):
		'''
		Load a fresh SimSnap with only the particle in halo i
		'''
		raise NotImplementedError("Requires pynbody functionality")

	def _get_halo(self, i):
		if self.base is None:
			raise RuntimeError("Parent snapshot has been deleted")

		return self._halos[i]

	def _get_by_id(self, halo_id):
		#Nasty! But, allows lookup by id only (do we need this?)
		idx = np.where(self._tree[:]['id'] == halo_id)[0][0]
		#print 'hn=', hn
		halo = self._halos[idx]
		return halo

	def find_progs(self, halo_list, all_halos, aexp_end=None):
		'''
		Recursively find all progenitors of a given halo
		If start_halo is not none, speficy the index to start with
		'''
		for halo in halo_list:
			if aexp_end and (halo['scale'] < aexp_end):
				break
			all_halos.append(halo)
			idxs = np.where(self._tree[:]['desc_id'] == halo['id'])[0]
			self.find_progs(self._tree[idxs], all_halos)

	def sub_catalogue(self, field, condition):
		'''
		Creata a sub-cataogue of halos
		'''
		idxs = np.where(self._tree[:][field] == condition)
		return RockstarSubCatalogue(self, self._tree[idxs])

	def sub_catalogue_aexp(self, aexp):
		'''
		Find closest available expansion factor
		'''
		all_aexp = self._tree[:]['scale']
		idx = (np.abs(all_aexp-aexp)).argmin()
		aexp_to_use = self._tree[idx]['scale']
		idxs = np.where(self._tree[:]['scale'] == aexp_to_use)
		return RockstarSubCatalogue(self, self._tree[idxs])

	@staticmethod
	def _can_load(simulation, **kwargs):
		fname = '%s/%s/trees/tree_0_0_0.dat'%(simulation.path(), pp_cfg.rockstar_base)
		print fname
		for file in glob.glob(fname):
			if os.path.exists(file):
				return True
		return False

	def get_offset(self):
		'''
		Determine the offset between snap_num in merger tree and actual
		'''
		simulation = self._base()
		num_snapshots = simulation.num_snapshots()

		#Take the first entry in the tree to find last snapshot
		last_snap = self[0]['snap_num'].value
		return num_snapshots - last_snap

class RockstarSubCatalogue(MergerTree):

	halo_type = RockstarMergerTree.halo_type
	units = RockstarMergerTree.units

	def __init__(self, parent_catalogue, tree):
		'''
		Index a sub-cataogue of a RockstarMergerTree
		'''
		MergerTree.__init__(self)
		self._base = parent_catalogue._base
		self._load_consistent_tree(tree)
		self._parent = parent_catalogue

	def _load_consistent_tree(self, tree):
		self._tree = tree
		self._nhalo = len(tree)

		for h in xrange(self._nhalo):
			self._halos[h] = MergerHalo(self._tree[h]['id'], self)
			self._halos[h].properties = self._tree[h]
		print 'done!'

	def _get_halo(self, i):
		if self.base is None:
			raise RuntimeError("Parent snapshot has been deleted")

		return self._halos[i]

	def _get_by_id(self, halo_id):
		#Nasty! But, allows lookup by id only (do we need this?)
		idx = np.where(self._tree[:]['id'] == halo_id)[0][0]
		#print 'hn=', hn
		halo = self._halos[idx]
		return halo

	@property
	def base(self):
		return self._base()

def open_(filename, *args):
	import gzip
	"""Open a file, determining from the filename whether to use
	gzip decompression"""

	if (filename[-3:] == '.gz'):
		return gzip.open(filename, *args)
	try:
		return open(filename, *args)
	except IOError:
		return gzip.open(filename + ".gz", *args)

def cutgz(x):
	"""Strip the .gz ending off a string"""
	if x[-3:] == '.gz':
		return x[:-3]
	else:
		return x
