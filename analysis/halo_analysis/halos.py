'''
Heavily based on halos.py from pynbody, so credit to those guys
Rewritten to allow loading of Rockstar catalogues when using any module

Currently depends on yt-3.0 for unit handling. I'll try and code this
dependency out sometime in the future if needed.

@author dsullivan
'''

super_verbose=True

import numpy as np
import sys, os.path, glob
import weakref, copy
from ramses_pp import config as pp_cfg
#from yt.units.yt_array import YTArray
#import logging

#Should adopt this throughout
#logger = logging.getLogger('ramses_pp.halo')

class DummyHalo(object):

	def __init__(self):
		self.properties = {}

class Halo():

	"""
	Generic class representing a halo.
	TODO - Units registry and return object with units
	"""

	def __init__(self, halo_id, halo_catalogue, *args):
		self.properties = {}
		self._halo_catalogue = halo_catalogue
		self._halo_id = halo_id
		self._descriptor = "halo_" + str(halo_id)

	def __getitem__(self, item):
		snap = self._halo_catalogue._base()
		ds = snap.raw_snapshot()
		unit = self._halo_catalogue.units[item]
		#return np.array([{item:self.properties[item]}, {'unit':unit}])
		#return YTArray(self.properties[item], unit)
		return ds.arr(self.properties[item], unit)

	def field_list(self):
		return self._halo_catalogue.halo_type

	def is_subhalo(self, otherhalo):
		"""
		Convenience function that calls the corresponding function in
		a halo catalogue.
		"""

		return self._halo_catalogue.is_subhalo(self._halo_id, otherhalo._halo_id)

	def get_sphere(self):
		'''
		YT only function. Currently, there is a bug with ds.sphere which ignores / h units.
		So for the time being this function will take care of that
		'''
		centre = self['pos']
		Rvir = self['Rvir']
		snapshot = self._halo_catalogue._base()
		if not (str(type(snapshot)) == "<class 'ramses_pp.modules.yt.YT.YTSnapshot'>"):
			raise NotImplementedError("sphere only implemented for YT")
		ds = snapshot.raw_snapshot()

		#if(centre.uq == YTArray(1.0, 'Mpc/h') or (Rvir.uq == YTArray(1.0, 'kpc/h'))):
		#	if super_verbose:
		#		print "Warning: YT ds.sphere() but does not recognise /h units. Correcting, but double check!"
		#	return ds.sphere(centre/ds.hubble_constant, Rvir/ds.hubble_constant)

		return ds.sphere(centre, Rvir)

	def dbscan(self, eps=0.4, min_samples=20):
		'''
		Define a sphere centred on the halo, and identify groups of stars using DBSCAN (i.e to find galaxies)
		eps defines the linking length to use, min_samples defines the minimum number of stars to define a galaxyTODO - normalise positions AND eps by the virial radius of the halo
		'''
		#First, we will need a sphere
		sphere = self.get_sphere()

		#Now, identify only star particles
		stars = sphere.particles.source['particle_age'] != 0
		#Make sure there are some stars to analyze
		if not True in stars:
			raise Exception("No star particles found in sphere")

		#Get the positions in code units
		if super_verbose:
			import time
			t0 = time.time()
			print 'Loading positions...'
		x, y, z = sphere['particle_position_x'][stars].value, \
						sphere['particle_position_y'][stars].value, \
						sphere['particle_position_z'][stars].value

		#Normalise these based on the virial radius
		pos = self['pos'].convert_to_units('code_length').value
		Rvir = self['Rvir'].convert_to_units('code_length').value

		min_x, max_x = pos[0] - Rvir, pos[0] + Rvir
		min_y, max_y = pos[1] - Rvir, pos[1] + Rvir
		min_z, max_z = pos[2] - Rvir, pos[2] + Rvir

		x_norm = (x - min_x)/(max_x - min_x)
		y_norm = (y - min_y)/(max_y - min_y)
		z_norm = (z - min_z)/(max_z - min_z)
		#x_norm = x
		#y_norm = y
		#z_norm = z

		print min(x_norm), max(x_norm)
		print min(y_norm), max(y_norm)
		print min(z_norm), max(z_norm)

		X = np.dstack([x_norm, y_norm, z_norm])[0]
		if super_verbose:
			print 'Got positions'
			t1 = time.time()
			print 'Took %f s'%(t1 - t0)

		#We have the data to run dbscan
		from scipy.spatial import distance
		from sklearn.cluster import DBSCAN
		#from sklearn import metrics
		# Compute similarities
		if super_verbose:
			print 'Computing similarities'
			t0 = time.time()
		D = distance.squareform(distance.pdist(X))
		S = 1 - (D / np.max(D))
		if super_verbose:
			print 'Got similarities'
			t1 = time.time()
			print 'Took %f s'%(t1 - t0)

		if super_verbose:
			print 'Starting DBSCAN'
			t0 = time.time()
		db = DBSCAN(eps=eps, min_samples=min_samples).fit(S)
		if super_verbose:
			print 'Finished DBSCAN'
			t1 = time.time()
			print 'Took %f s'%(t1 - t0)
		del(D)
		del(S)
		#core_samples = db.core_sample_indices_
		labels = db.labels_

		# Number of clusters in labels, ignoring noise if present.
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
		print 'Found %d groups in halo %d'%(n_clusters_, self._halo_id)

		return db, X

class HaloCatalogue(object):

	'''
	Generic halo catalogue bases on pynbody, but with non-essential functionality
	stripped out
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
				if i > len(self._halos)-1:
					break
			except RuntimeError:
				break

	def is_subhalo(self, childid, parentid):
		"""
		Checks whether the specified 'childid' halo is a subhalo
		of 'parentid' halo.
		"""
		parent = self._get_by_id(parentid)
		#print parent['id'] == parentid
		return (childid in parent._children)

	def __contains__(self, haloid):
		return self.contains(haloid)

	@staticmethod
	def can_load(self):
		return False

	#TODO - Move all running of halo cats into here
	@staticmethod
	def can_run(self):
		return False

	def sort(self, field, reverse=True):
		'''
		Sort halos by a given field
		'''
		return sorted(self, key = lambda x: x[field], reverse=reverse)

	def mass_function(self, units='Msun/h', nbins=100):
		'''
		Compute the halo mass function for the given catalogue
		'''
		masses = []
		for halo in self:
			Mvir = halo['Mvir'].in_units(units)
			masses.append(Mvir)

		mhist, mbin_edges = np.histogram(np.log10(masses),bins=nbins)
		mbinmps = np.zeros(len(mhist))
		mbinsize = np.zeros(len(mhist))
		for i in np.arange(len(mhist)):
			mbinmps[i] = np.mean([mbin_edges[i],mbin_edges[i+1]])
			mbinsize[i] = mbin_edges[i+1] - mbin_edges[i]

		return mbinmps, mhist, mbinsize


#Rockstar
class RockstarCatalogue(HaloCatalogue):
	'''
	Class to handle catalogues produced by Rockstar (by Peter Behroozi).
	'''	

	halo_type = np.dtype([('id',np.int64),('hostHalo',np.int64),('Mvir','f'),
					  ('v_max','f'),('v_rms','f'),('Rvir','f'),
					  ('Rs','f'),('num_p',np.int64),
					  ('pos','f',3),('vel','f',3),
					  ('J','f',3),('spin','f'),('klypin_rs','f'),
					  ('Mvir_all','f'),
					  ('M200b','f'),('M200c','f'),('M500c','f'),
					  ('Xoff','f'),('Voff','f'),
					  ('bullock_spin','f'),('b_to_a','f'),('c_to_a','f'),
					  ('A','f',3),('T/|U|','f')])

	mtree_type = np.dtype([('scale','f'),('id',np.int64),('desc_scale','f'),('desc_id',np.int64),
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

	#cm = comoving
	#nocm = physical
	units = {'id':'dimensionless',
			'hostHalo':'dimensionless',
			'Mvir':'Msun / h',
			'v_max':'km / s',
			'v_rms':'km / s',
			'Rvir':'kpccm / h',
			'Rs':'kpccm / h',
			'num_p':'dimensionless',
			'pos':'Mpccm / h',
			'vel':'km / s',
			'J':'(Msun/h)**2 / (km/s)',
			'spin':'dimensionless',
			'klypin_rs':'dimensionless',
			'Mvir_all':'Msun / h',
			'M200b':'Msun / h',
			'M200c':'Msun / h',
			'M500c':'Msun / h',
			'Xoff':'kpccm / h',
			'Voff':'km / s',
			'bullock_spin':'dimensionless',
			'b_to_a':'kpccm / h',
			'c_to_a':'kpccm / h',
			'A':'?',
			'T/|U|':'dimensionless'}

	def __init__(self, snap, filename=None, make_grp=None):
		# TODO - Read/store header
		if not self._can_load(snap):
			raise Exception("Cannot locate/load rockstar catalogue")

		self._base = weakref.ref(snap)
		#self._base = snap
		HaloCatalogue.__init__(self)
		
		if filename is not None: self._rsFilename = filename
		else:
			fname = '%s/%s/out_%d.list'%(snap.path(), pp_cfg.rockstar_base, snap.output_number()-1)
			self._rsFilename = cutgz(glob.glob(fname)[0])
			if True == pp_cfg.verbose:
				print 'Loading from %s'%fname

		#Try and open the files
		try:
			f = open_(self._rsFilename)
		except IOError:
			raise IOError("Halo catalogue not found -- check the file name of catalogue data or try specifying a catalogue using the filename keyword")

		#self._head = np.fromstring(f.read(self.head_type.itemsize),
		#	dtype=self.head_type)
		#unused = f.read(256 - self._head.itemsize)

		#self._nhalos = self._head['num_halos'][0]

		if True == pp_cfg.verbose:
			print "RockstarCatalogue: loading halos...",
			sys.stdout.flush()

		self._load_rs_halos(f, snap)
		f.close()

		self._setup_children()

		if make_grp is None:
			make_grp = pp_cfg.rockstar_autogrp

		if make_grp:
			self.make_grp()

	def __str__(self):
		return 'HaloCatalogue - Snapshot %s'%self._base()

	def make_grp(self):
		'''
		Creates a 'grp' array which labels each particle according to
		its parent halo.
		'''
		#try:
		#	self.base['grp']
		#except:
		#	self.base['grp'] = np.zeros(len(self.base),dtype='i')

		#for halo in self._halos.values():
		#	halo[name][:] = halo._halo_id

		#if config['verbose']:  print "writing %s"%(self._base().filename+'.grp')
		#self._base().write_array('grp',overwrite=True,binary=False)
		raise NotImplementedError("Requires pynbody functionality")

	def _setup_children(self):
		'''
		Creates a 'children' array inside each halo's 'properties'
		listing the halo IDs of its children. Used in case the reading
		of substructure data from the AHF-supplied _substructure file
		fails for some reason.
		'''
		for i in xrange(self._nhalos):
			self._halos[i]._children = []

		for i in xrange(self._nhalos):
			host = self._halos[i].properties['hostHalo']
			if host > -1:
				try:
					self._halos[host]._children.append(i)
				except KeyError:
					pass

	def _get_halo(self, i):
		if self.base is None:
			raise RuntimeError("Parent snapshot has been deleted")

		return self._halos[i]

	def _get_by_id(self, halo_id):
		#Nasty! But, allows lookup by id only (do we need this?)
		idx = np.where(self._haloprops[:]['id'] == halo_id)[0][0]
		hn = np.where(self._num_p_rank==idx)[0][0]
		#print 'hn=', hn
		halo = self._halos[hn]
		return halo

	@property
	def base(self):
		return self._base()

	def _load_rs_halos(self, f, snap):
		haloprops = np.loadtxt(f, dtype=self.halo_type, comments='#')
		self._nhalos = len(haloprops)
		self._haloprops = np.array(haloprops)
		# sort by number of particles to make compatible with AHF
		self._num_p_rank = np.flipud(self._haloprops[:]['num_p'].argsort(axis=0))

		for h in xrange(self._nhalos): # self._nhalos + 1?
			hn = np.where(self._num_p_rank==h)[0][0]

			#Is this really memory inefficient?
			self._halos[hn] = Halo(self._haloprops[h]['id'], self)
			# properties are in Msun / h, Mpc / h
			self._halos[hn].properties = self._haloprops[h]

	def load_copy(self, i):
		'''
		Load a fresh SimSnap with only the particle in halo i
		'''
		raise NotImplementedError("Requires pynbody functionality")

	def _load_rs_particles(self, f, snap):
		NotImplementedError("Only halo loading implemented")

	@staticmethod
	def _can_load(snap, **kwargs):
		fname = '%s/%s/out_%d.list'%(snap.path(), pp_cfg.rockstar_base, snap.output_number()-1)
		print fname
		for file in glob.glob(fname):
			if os.path.exists(file):
				return True
		return False
		

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