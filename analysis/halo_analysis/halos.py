'''
Heavily based on halos.py from pynbody, so credit to those guys
Rewritten to allow loading of Rockstar catalogues when using any module

Currently depends on yt-3.0 for unit handling. I'll try and code this
dependency out sometime in the future if needed.

@author dsullivan, bthompson
'''

super_verbose=True
OFFSET=30 # This is to match rockstar outputs to snapshot numbers (in the event that cayalogue does not start at 1) - manual at the moment (I cba to fix it)

import yt
import numpy as np
import sys, os.path, glob
import weakref, copy
import re, uuid
from ramses_pp import config as pp_cfg
from yt.units.yt_array import YTArray
#from ..analysis import center_finder
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

	def get_sphere(self, Rvir_fact=1.0):
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
	#	print Rvir

#		if(centre.uq == YTArray(1.0, 'Mpc/h') or (Rvir.uq == YTArray(1.0, 'kpc/h'))):
#			if super_verbose:
#				print "Warning: YT ds.sphere() but does not recognise /h units. Correcting, but double check!"
#			return ds.sphere(centre/ds.hubble_constant, Rvir/ds.hubble_constant)

		return ds.sphere(centre, (Rvir* Rvir_fact))

	def projection_plot(self, Rvir_fact=1.0,fixed_length=None,stars=False,axis="x",units="kpc/h",quantity="density"):
		'''
		YT only function, produce and save a gas projection plot of the Halo of interest
		Rvir_fact = multiply Rvir or fixed_length by a factor of something
		fixed_length: if None, then use Rvir... otherwise whatever custom length (with units) which you choose
		axis = projection axis
		units = projection axis units
		quantity = whatever you feel like plotting
		'''
		halo_sphere = self.get_sphere(Rvir_fact=Rvir_fact)
		if fixed_length == None:
			fixed_length = self["Rvir"]
		region = halo_sphere.ds.region(center=self["pos"],left_edge=[(i.in_units('code_length').value - (self["Rvir"].in_units('code_length').value * Rvir_fact ) ) for i in self["pos"] ], right_edge=[(i.in_units('code_length').value + (self["Rvir"].in_units('code_length').value * Rvir_fact ) ) for i in self["pos"] ])
		
		plt = yt.ProjectionPlot(halo_sphere.ds,axis,quantity,center=self["pos"].in_units(units), width=(fixed_length.in_units('code_length').value*Rvir_fact), axes_unit=units, data_source=region )
		plt.save(("Halo_" + str(self["id"].value)) )

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
		x, y, z, ind = sphere['particle_position_x'][stars].value, \
						sphere['particle_position_y'][stars].value, \
						sphere['particle_position_z'][stars].value, \
						sphere['particle_index'][stars].value

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

		# add the particle indices to the origional array
		Y = np.insert(X,3,ind,axis=1)

		# Number of clusters in labels, ignoring noise if present.
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
		print 'Found %d groups in halo %d'%(n_clusters_, self._halo_id)
		return db, Y

	def iterative_com(self, stars=True, dark=False, gas=False):
		"""
		find the center of the galaxy in the halo in a quick and dirty method (credit Gareth Few for the idea)
		1) put a sphere around the halo at Rvir = 1 centered on centre of halo
		2) find stellar CoM
		3) put a sphere around 0.9 Rvir of the halo on the stellar center
		4) find stellar CoM
		5) repeat until 0.1Rvir
		"""

		center = self['pos']
		Rvir = self['Rvir']
		snapshot = self._halo_catalogue._base()
		if not (str(type(snapshot)) == "<class 'ramses_pp.modules.yt.YT.YTSnapshot'>"):
			raise NotImplementedError("sphere only implemented for YT")
		ds = snapshot.raw_snapshot()

		## from stars, dark, gas .. what are we trying to find the centre of?
		## cases
		## if gas = True, stars = False, Dark = False ... case = 1 = gas CoM
		## if gas = False, stars = True, Dark = False ... case = 2 = star CoM
		## if gas = False, stars = False, Dark = True ... case = 3 = Dark CoM
		## if gas = True, stars = True, Dark = False ... case = 4 = Baryon CoM
		## if gas = False, stars = True, Dark = True ... case = 5 = Dark Star CoM
		## if gas = True, stars = False, Dark = True ... case = 6 = Dark Gas
		## if gas = True, stars = True, Dark = True ... case = 7 = All CoM
		## if gas = False, stars = False, Dark = False .. case = 8... nothing... return

		if stars == False and dark == False and gas == False:
			print "Finding the centre of nothing is boring.. returning nothing"
			return None
		

		for i in range(0,5):
			print i
			Rvir_fact = 1.0 - (i * 0.2) # shrink virial radius by 0.2 upon each run.. center on the new CoM
			try:
				sphere = ds.sphere(center, (Rvir*Rvir_fact))
			except yt.utilities.exceptions.YTSphereTooSmall:
				break
#			stars = None # potentially fixes a bug in YT
			stars_part = sphere.particles.source['particle_age'] != 0
			dark_part = sphere.particles.source['particle_age'] == 0
			print "stars ", len(stars_part)
#			center = sphere[stars].quantities.center_of_mass(use_gas = gas) ## find the CoM of the stars
#			print center

			# case 1

			if stars == False and gas == True and dark == False:
				center = sphere.quantities.center_of_mass(use_gas=True, use_particles=False)

			# case 2

			elif stars == True and gas == False and dark == False:
				com_x = (sphere["particle_position_x"][stars_part] * sphere["particle_mass"][stars_part]).sum(dtype=np.float64) / (sphere["particle_mass"][stars_part]).sum(dtype=np.float64)
				com_y = (sphere["particle_position_y"][stars_part] * sphere["particle_mass"][stars_part]).sum(dtype=np.float64) / (sphere["particle_mass"])[stars_part].sum(dtype=np.float64)
				com_z = (sphere["particle_position_z"][stars_part]* sphere["particle_mass"][stars_part]).sum(dtype=np.float64) / (sphere["particle_mass"][stars_part]).sum(dtype=np.float64)
				center = [com_x,com_y,com_z]
	
			# case 3
			elif stars == False and gas == False and dark == True:
				com_x = (sphere["particle_position_x"][dark_part] * sphere["particle_mass"][dark_part]).sum(dtype=np.float64) / (sphere["particle_mass"][dark_part]).sum(dtype=np.float64)
				com_y = (sphere["particle_position_y"][dark_part] * sphere["particle_mass"][dark_part]).sum(dtype=np.float64) / (sphere["particle_mass"])[dark_part].sum(dtype=np.float64)
				com_z = (sphere["particle_position_z"][dark_part]* sphere["particle_mass"][dark_part]).sum(dtype=np.float64) / (sphere["particle_mass"][dark_part]).sum(dtype=np.float64)
				center = [com_x,com_y,com_z]
			# case 4 baryon COM

			elif stars == True and gas == True and dark == False:
				com_x = (  (sphere["particle_position_x"][stars_part] * sphere["particle_mass"][stars_part]).sum(dtype=np.float64) + \
					(sphere["index","x"] * sphere["gas","cell_mass"]).sum(dtype=np.float64) )  / ((sphere["particle_mass"][stars_part]).sum(dtype=np.float64) + (sphere["gas","cell_mass"]).sum(dtype=np.float64) )
				com_y = ( ( (sphere["particle_position_y"][stars_part] * sphere["particle_mass"][stars_part]).sum(dtype=np.float64) + \
					(sphere["index","y"] * sphere["gas","cell_mass"]).sum(dtype=np.float64) ) )  / ((sphere["particle_mass"][stars_part]).sum(dtype=np.float64) + (sphere["gas","cell_mass"]).sum(dtype=np.float64) )
				com_z = ( (sphere["particle_position_z"][stars_part] * sphere["particle_mass"][stars_part]).sum(dtype=np.float64) +  \
					(sphere["index","z"] * sphere["gas","cell_mass"]).sum(dtype=np.float64) )   / ( (sphere["particle_mass"][stars_part]).sum(dtype=np.float64) + (sphere["gas","cell_mass"]).sum(dtype=np.float64) )
				center = [com_x,com_y,com_z]

			# case 5 dark star COM
			elif stars == True and gas == False and dark == True:
				center = sphere.quantities.center_of_mass(use_gas=False, use_particles=True) 

			# case 6

			elif stars == False and gas == True and dark == True:
				com_x = (  (sphere["particle_position_x"][dark_part] * sphere["particle_mass"][dark_part]).sum(dtype=np.float64) + \
					(sphere["index","x"] * sphere["gas","cell_mass"]).sum(dtype=np.float64) )  / ((sphere["particle_mass"][dark_part]).sum(dtype=np.float64) + (sphere["gas","cell_mass"]).sum(dtype=np.float64) )
				com_y = ( ( (sphere["particle_position_y"][dark_part] * sphere["particle_mass"][dark_part]).sum(dtype=np.float64) + \
					(sphere["index","y"] * sphere["gas","cell_mass"]).sum(dtype=np.float64) ) )  / ((sphere["particle_mass"][dark_part]).sum(dtype=np.float64) + (sphere["gas","cell_mass"]).sum(dtype=np.float64) )
				com_z = ( (sphere["particle_position_z"][dark_part] * sphere["particle_mass"][dark_part]).sum(dtype=np.float64) +  \
					(sphere["index","z"] * sphere["gas","cell_mass"]).sum(dtype=np.float64) )   / ( (sphere["particle_mass"][dark_part]).sum(dtype=np.float64) + (sphere["gas","cell_mass"]).sum(dtype=np.float64) )
				center = [com_x,com_y,com_z]

			

			# case 7 All COM
			else:
				center = sphere.quantities.center_of_mass(use_gas=True, use_particles=True)

			print center
			

		print "centre found at ", center, " with stars " + str(stars) + " dark " + str(dark) + " gas " + str(gas)
		return center

	def galaxy_disk(self,n_disks=5,disk_h=(10,"kpc/h"),disk_w=(50,"kpc/h"),center=None):
		# n_disks need to be odd... you have 1 disk.. or 1 disk sandwiched inbetween 2 others (n_disk=3) etc
		if center == None:
			center = self.iterative_com(stars=True, dark=False, gas=True)
		if n_disks % 2 == 0: # check for even number
			print "n_disks must be odd"
			return

		print disk_h[0], disk_h[1]
		disk_h = YTArray(disk_h[0],disk_h[1])
		print disk_h
		# get angular momentum
		Rvir = self['Rvir']
		snapshot = self._halo_catalogue._base()
		if not (str(type(snapshot)) == "<class 'ramses_pp.modules.yt.YT.YTSnapshot'>"):
			raise NotImplementedError("sphere only implemented for YT")
		ds = snapshot.raw_snapshot()

		for i in range(0,10):
			Rvir_fact = 1.0 - (i * 0.1) # shrink virial radius by 0.1 upon each run.. center on the new CoM
			try:
				sphere = ds.sphere(center, (Rvir*Rvir_fact))
			except yt.utilities.exceptions.YTSphereTooSmall:
				break
			L = sphere.quantities.angular_momentum_vector(use_gas=True,use_particles=False)
			print L
		# simple galaxy disk for height and width work
	

		L_mag = np.sqrt(L[0].value**2 + L[1].value**2 + L[2].value**2)
		print L_mag
		L_norm = np.zeros(3)
		L_norm[0] = L[0].value/L_mag
		L_norm[1] = L[1].value/L_mag
		L_norm[2] = L[2].value/L_mag
		print L_norm

		# bin of data disks

		# how small can we get the height
                disks = []

		stacks = (n_disks + 1)/2

	#	print center[0].convert_to_units("kpc/h")
		center = [center[0].convert_to_units("kpc/h"),center[1].convert_to_units("kpc/h"),center[2].convert_to_units("kpc/h")]

		print center	
		print disk_h	
		for i in range(0,stacks):
			if i == 0: # central region first
				cent = center
				h_bin = YTArray(0,"kpc/h")
				d_obj = ds.disk(center,L,Rvir,disk_h)
				disks.insert(0,[h_bin,cent,d_obj])
			else:
				cent_up = [center[0] + ((L_norm[0] * i * disk_h)), center[1] + ((L_norm[1]* i * disk_h)) , center[2] + ((L_norm[2] * i * disk_h))]
				cent_down = [center[0] - ((L_norm[0] * i * disk_h)), center[1] - ((L_norm[1]* i * disk_h)) , center[2] - ((L_norm[2] * i * disk_h))]
				h_bin_up = i * disk_h
				h_bin_down = i *(-disk_h)
				d_obj_up = ds.disk(cent_up,L,Rvir,disk_h)
				d_obj_down = ds.disk(cent_down,L,Rvir,disk_h)
				disks.insert(0,[h_bin_up,cent_up,d_obj_up])
				disks.append([h_bin_down,cent_down,d_obj_down])

		print disks
		cylinder = ds.disk(center,L,disk_w,Rvir)
		return disks, cylinder

		
		
		



	def stellar_center(self):
		print "not functioning yet, will do something with dbscan results"
		return 

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
				if i > len(self._halos):
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
		hn = np.where(self._num_p_rank==idx)[0][0]+1
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

# AHF

class AHFCatalogue(HaloCatalogue):
	'''
	Class to handle catalogues produced by Rockstar (by Peter Behroozi).
	'''	
## assume file structure like this

#ID(1)  hostHalo(2)     numSubStruct(3) Mvir(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) Rvir(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)        lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24)  b(25)   c(26)   Eax(27) Eay(28) Eaz(29) Ebx(30) Eby(31) Ebz(32) Ecx(33) Ecy(34) Ecz(35) ovdens(36)      nbins(37)       fMhires(38)     Ekin(39)        Epot(40)        SurfP(41)       Phi0(42)        cNFW(43)        n_gas(44)       M_gas(45)       lambda_gas(46)  lambdaE_gas(47) Lx_gas(48)      Ly_gas(49)      Lz_gas(50)      b_gas(51)       c_gas(52)       Eax_gas(53)     Eay_gas(54)     Eaz_gas(55)     Ebx_gas(56)     Eby_gas(57)     Ebz_gas(58)     Ecx_gas(59)     Ecy_gas(60)     Ecz_gas(61)     Ekin_gas(62)    Epot_gas(63)    n_star(64)      M_star(65)      lambda_star(66) lambdaE_star(67)        Lx_star(68)     Ly_star(69)     Lz_star(70)     b_star(71)      c_star(72)      Eax_star(73)    Eay_star(74)    Eaz_star(75)    Ebx_star(76)    Eby_star(77)    Ebz_star(78)    Ecx_star(79)    Ecy_star(80)    Ecz_star(81)    Ekin_star(82)   Epot_star(83)
	halo_type =  np.dtype([('id',np.int64),('hostHalo',np.int64),('numSubStruct',np.int64),('Mvir','f'),
					  ('num_p',np.int64),('pos','f',3),('vel','f',3),
					  ('Rvir','f'),('Rmax','f'),('r2','f'),('mpb_offset','f'),
					  ('com_offset','f'),('v_max','f'),('v_esc','f'),('sigV','f'),
					  ('bullock_spin','f'),('spin','f'),('L','f',3),('b','f'),('c','f'),
					  ('Ea','f',3),('Eb','f',3),('Ec','f',3),('ovdens','f'),
					  ('nbins',np.int64),('fMhires','f'),('Ekin','f'),
					  ('Epot','f'),('SurfP','f'),('PhiO','f'),('cNFW','f'),('n_gas',np.int64),
					  ('M_gas','f'),('bullock_spin_gas','f'),('spin_gas','f'),
					  ('L_gas','f',3),('b_gas','f'),('c_gas','f'),
					  ('Ea_gas','f',3),('Eb_gas','f',3),('Ec_gas','f',3,),
					  ('Ekin_gas','f'),('Epot_gas','f'),('n_star',np.int64),
					  ('M_star','f'),('bullock_spin_star','f'),('spin_star','f'),
					  ('L_star','f',3),('b_star','f'),('c_star','f'),
					  ('Ea_star','f',3),('Eb_star','f',3),('Ec_star','f',3,),
					  ('Ekin_star','f'),('Epot_star','f')])

	units = {'id':'dimensionless',
			'hostHalo':'dimensionless',
			'numSubStruct':'dimensionless',
			'Mvir':'Msun / h',
			'num_p':'dimensionless',
			'pos':'kpc / h',
			'vel':'km / s',
			'Rvir':'kpc / h',
			'Rmax':'kpc / h',
			'r2':'kpc / h',
			'mpb_offset': 'kpc / h',
			'com_offset': 'kpc / h',
			'v_max':'km / s',
			'v_esc':'km / s',
			'sigV': 'km / s',
			'bullock_spin': 'dimensionless',
			'spin': 'dimensionless',
			'L':'dimensionless',
			'b_to_a':'kpc / h',
			'c_to_a':'kpc / h',
			'Ea':'dimensionless',
			'Eb':'dimensionless',
			'Ec':'dimensionless',
			'ovdens':'dimensionless',
			'nbins':'dimensionless',
			'fMhires':'dimensionless',
			'Ekin': 'Msun / h (km / sec)**2',
			'Epot': 'Msun / h (km / sec)**2',
			'SurfP': 'Msun / h (km / sec)**2',
			'PhiO': '(km / sec)**2',
			'cNFW': 'dimensionless',
			'M_gas': 'Msun / h',
			'bullock_spin_gas': 'dimensionless',
			'spin_gas': 'dimensionless',
			'L_gas':'dimensionless',
			'b_to_a_gas':'kpc / h',
			'c_to_a_gas':'kpc / h',
			'Ea_gas':'dimensionless',
			'Eb_gas':'dimensionless',
			'Ec_gas':'dimensionless',
			'n_gas':'dimensionless',
			'Ekin_gas': 'Msun / h (km / sec)**2',
			'Epot_gas': 'Msun / h (km / sec)**2',
			'M_star': 'Msun / h',
			'bullock_spin_star': 'dimensionless',
			'spin_star': 'dimensionless',
			'L_star':'dimensionless',
			'b_to_a_star':'kpc / h',
			'c_to_a_star':'kpc / h',
			'Ea_star':'dimensionless',
			'Eb_star':'dimensionless',
			'Ec_star':'dimensionless',
			'n_star':'dimensionless',
			'Ekin_star': 'Msun / h (km / sec)**2',
			'Epot_star': 'Msun / h (km / sec)**2',
		}

#	halo_type = np.dtype([('id',np.int64),('hostHalo',np.int64),('Mvir','f'),
#					  ('v_max','f'),('v_rms','f'),('Rvir','f'),
#					  ('Rs','f'),('num_p',np.int64),
#					  ('pos','f',3),('vel','f',3),
#					  ('J','f',3),('spin','f'),('klypin_rs','f'),
#					  ('Mvir_all','f'),
#					  ('M200b','f'),('M200c','f'),('M500c','f'),
#					  ('Xoff','f'),('Voff','f'),
#					  ('bullock_spin','f'),('b_to_a','f'),('c_to_a','f'),
#					  ('A','f',3),('T/|U|','f')])

#	units = {'id':'dimensionless','hostHalo':'dimensionless','Mvir':'Msun / h',
#					  'v_max':'km / s','v_rms':'km / s','Rvir':'kpc / h',
#					  'Rs':'kpc / h','num_p':'dimensionless',
#					  'pos':'Mpc / h','vel':'km / s',
#					  'J':'(Msun/h)**2 / (km/s)','spin':'dimensionless','klypin_rs':'dimensionless',
#					  'Mvir_all':'Msun / h',
#					  'M200b':'Msun / h','M200c':'Msun / h','M500c':'Msun / h',
#					  'Xoff':'kpc / h','Voff':'km / s',
#					  'bullock_spin':'dimensionless','b_to_a':'kpc / h','c_to_a':'kpc / h',
#					  'A':'?','T/|U|':'K(?)'}

	def __init__(self, snap, filename=None, make_grp=None):
		# TODO - Read/store header
		if not self._can_load(snap):
			raise Exception("Cannot locate/load AHF catalogue")

		self._base = weakref.ref(snap)
		HaloCatalogue.__init__(self)
		if filename is not None: self._AHFFilename = filename
		else:
			# get the file name
			tipsy_dir = str("%soutput_%05d_tipsy/" % (snap.path(), snap.output_number()))
			if not os.path.isdir(tipsy_dir):
				print "AHF not run on this snapshot.. aborting"
				raise Exception("AHF not run, run AHF plz")

			if len(glob.glob('%s*.AHF_particles'%tipsy_dir)) == 0:
				print "AHF not run on this snapshot.. aborting"
				raise Exception("AHF not run, run AHF plz")
				
			fname = glob.glob('%s*.AHF_halos'%tipsy_dir)[0]
			self._AHFFilename = fname
			print 'Loading from %s'%fname

		#Try and open the files
		try:
			f = open_(self._AHFFilename)
		except IOError:
			raise IOError("Halo catalogue not found -- check the file name of catalogue data or try specifying a catalogue using the filename keyword")

		#self._head = np.fromstring(f.read(self.head_type.itemsize),
		#	dtype=self.head_type)
		#unused = f.read(256 - self._head.itemsize)

		#self._nhalos = self._head['num_halos'][0]

		if super_verbose == True:
			print "AHFCatalogue: loading halos...",
			sys.stdout.flush()

		self._load_ahf_halos(self._AHFFilename, snap)
		f.close()

		self._setup_children()

#		if make_grp is None:
#			make_grp = pp_cfg.rockstar_autogrp

		if make_grp:
			print "make grp not implemented yet"
#			self.make_grp()

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


	def _load_ahf_halos(self, f, snap):
		haloprops = np.loadtxt(f, dtype=self.halo_type, comments='#')
		self._nhalos = len(haloprops)
		self._haloprops = np.array(haloprops)
		# sort by number of particles to make compatible with AHF
		self._num_p_rank = np.flipud(self._haloprops[:]['num_p'].argsort(axis=0))

		for h in xrange(self._nhalos): # self._nhalos + 1?
			# AHF halos index start from 0... python start from 0.. lets be consitant :P
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

	def _load_ahf_particles(self, f, snap):
		NotImplementedError("Only halo loading implemented")

	@staticmethod
	def _can_load(snap, **kwargs):
		tipsy_dir = str("%soutput_%05d_tipsy/" % (snap.path(), snap.output_number()))
		if not os.path.isdir(tipsy_dir):
			return False

		if len(glob.glob('%s*.AHF_particles'%tipsy_dir)) == 0:
			return False
		
		fname = glob.glob('%s*.AHF_halos'%tipsy_dir)[0]
		if os.path.exists(fname):
			return True
		return False

class AHFHaloTracker(HaloCatalogue):
	'''
	Class to handle catalogues produced by Rockstar (by Peter Behroozi).
	'''	
## assume file structure like this

	halo_type =  np.dtype([('z','f'),('id',np.int64),('hostHalo',np.int64),('numSubStruct',np.int64),('Mvir','f'),
					  ('num_p',np.int64),('pos','f',3),('vel','f',3),
					  ('Rvir','f'),('Rmax','f'),('r2','f'),('mpb_offset','f'),
					  ('com_offset','f'),('v_max','f'),('v_esc','f'),('sigV','f'),
					  ('bullock_spin','f'),('spin','f'),('L','f',3),('b','f'),('c','f'),
					  ('Ea','f',3),('Eb','f',3),('Ec','f',3),('ovdens','f'),
					  ('nbins',np.int64),('fMhires','f'),('Ekin','f'),
					  ('Epot','f'),('SurfP','f'),('PhiO','f'),('cNFW','f'),('n_gas',np.int64),
					  ('M_gas','f'),('bullock_spin_gas','f'),('spin_gas','f'),
					  ('L_gas','f',3),('b_gas','f'),('c_gas','f'),
					  ('Ea_gas','f',3),('Eb_gas','f',3),('Ec_gas','f',3,),
					  ('Ekin_gas','f'),('Epot_gas','f'),('n_star',np.int64),
					  ('M_star','f'),('bullock_spin_star','f'),('spin_star','f'),
					  ('L_star','f',3),('b_star','f'),('c_star','f'),
					  ('Ea_star','f',3),('Eb_star','f',3),('Ec_star','f',3,),
					  ('Ekin_star','f'),('Epot_star','f')])

	units = {'z':'dimensionless',
			'id':'dimensionless',
			'hostHalo':'dimensionless',
			'numSubStruct':'dimensionless',
			'Mvir':'Msun / h',
			'num_p':'dimensionless',
			'pos':'kpc / h',
			'vel':'km / s',
			'Rvir':'kpc / h',
			'Rmax':'kpc / h',
			'r2':'kpc / h',
			'mpb_offset': 'kpc / h',
			'com_offset': 'kpc / h',
			'v_max':'km / s',
			'v_esc':'km / s',
			'sigV': 'km / s',
			'bullock_spin': 'dimensionless',
			'spin': 'dimensionless',
			'L':'dimensionless',
			'b_to_a':'kpc / h',
			'c_to_a':'kpc / h',
			'Ea':'dimensionless',
			'Eb':'dimensionless',
			'Ec':'dimensionless',
			'ovdens':'dimensionless',
			'nbins':'dimensionless',
			'fMhires':'dimensionless',
			'Ekin': 'Msun / h (km / sec)**2',
			'Epot': 'Msun / h (km / sec)**2',
			'SurfP': 'Msun / h (km / sec)**2',
			'PhiO': '(km / sec)**2',
			'cNFW': 'dimensionless',
			'M_gas': 'Msun / h',
			'bullock_spin_gas': 'dimensionless',
			'spin_gas': 'dimensionless',
			'L_gas':'dimensionless',
			'b_to_a_gas':'kpc / h',
			'c_to_a_gas':'kpc / h',
			'Ea_gas':'dimensionless',
			'Eb_gas':'dimensionless',
			'Ec_gas':'dimensionless',
			'n_gas':'dimensionless',
			'Ekin_gas': 'Msun / h (km / sec)**2',
			'Epot_gas': 'Msun / h (km / sec)**2',
			'M_star': 'Msun / h',
			'bullock_spin_star': 'dimensionless',
			'spin_star': 'dimensionless',
			'L_star':'dimensionless',
			'b_to_a_star':'kpc / h',
			'c_to_a_star':'kpc / h',
			'Ea_star':'dimensionless',
			'Eb_star':'dimensionless',
			'Ec_star':'dimensionless',
			'n_star':'dimensionless',
			'Ekin_star': 'Msun / h (km / sec)**2',
			'Epot_star': 'Msun / h (km / sec)**2',
		}

#	halo_type = np.dtype([('id',np.int64),('hostHalo',np.int64),('Mvir','f'),
#					  ('v_max','f'),('v_rms','f'),('Rvir','f'),
#					  ('Rs','f'),('num_p',np.int64),
#					  ('pos','f',3),('vel','f',3),
#					  ('J','f',3),('spin','f'),('klypin_rs','f'),
#					  ('Mvir_all','f'),
#					  ('M200b','f'),('M200c','f'),('M500c','f'),
#					  ('Xoff','f'),('Voff','f'),
#					  ('bullock_spin','f'),('b_to_a','f'),('c_to_a','f'),
#					  ('A','f',3),('T/|U|','f')])

#	units = {'id':'dimensionless','hostHalo':'dimensionless','Mvir':'Msun / h',
#					  'v_max':'km / s','v_rms':'km / s','Rvir':'kpc / h',
#					  'Rs':'kpc / h','num_p':'dimensionless',
#					  'pos':'Mpc / h','vel':'km / s',
#					  'J':'(Msun/h)**2 / (km/s)','spin':'dimensionless','klypin_rs':'dimensionless',
#					  'Mvir_all':'Msun / h',
#					  'M200b':'Msun / h','M200c':'Msun / h','M500c':'Msun / h',
#					  'Xoff':'kpc / h','Voff':'km / s',
#					  'bullock_spin':'dimensionless','b_to_a':'kpc / h','c_to_a':'kpc / h',
#					  'A':'?','T/|U|':'K(?)'}

	def __init__(self, snap, halo=None, filename=None, make_grp=None):
		# TODO - Read/store header
		if not self._can_load(snap):
			raise Exception("Cannot locate/load AHF catalogue")

		self._base = weakref.ref(snap)
		HaloCatalogue.__init__(self)
		if filename is not None: self._AHFFilename = filename
		else:
			# get the file name
			tipsy_dir = str("%soutput_%05d_tipsy/" % (snap.path(), snap.output_number()))
			if not os.path.isdir(tipsy_dir):
				print "AHF not run on this snapshot.. aborting"
				raise Exception("AHF not run, run AHF plz")

			if len(glob.glob('%s*.AHF_particles'%tipsy_dir)) == 0:
				print "AHF not run on this snapshot.. aborting"
				raise Exception("AHF not run, run AHF plz")
				
			fname = glob.glob('%s*.AHF_halos'%tipsy_dir)[0]
			self._AHFFilename = fname
			print 'Loading from %s'%fname

		#Try and open the files
		try:
			f = open_(self._AHFFilename)
		except IOError:
			raise IOError("Halo catalogue not found -- check the file name of catalogue data or try specifying a catalogue using the filename keyword")

		#self._head = np.fromstring(f.read(self.head_type.itemsize),
		#	dtype=self.head_type)
		#unused = f.read(256 - self._head.itemsize)

		#self._nhalos = self._head['num_halos'][0]

		if super_verbose == True:
			print "AHFCatalogue: loading halos...",
			sys.stdout.flush()

		self._load_ahf_halos(self._AHFFilename, snap)
		f.close()

		self._setup_children()

#		if make_grp is None:
#			make_grp = pp_cfg.rockstar_autogrp

		if make_grp:
			print "make grp not implemented yet"
#			self.make_grp()

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


	def _load_ahf_halos(self, f, snap):
		haloprops = np.loadtxt(f, dtype=self.halo_type, comments='#')
		self._nhalos = len(haloprops)
		self._haloprops = np.array(haloprops)
		# sort by number of particles to make compatible with AHF
		self._num_p_rank = np.flipud(self._haloprops[:]['num_p'].argsort(axis=0))

		for h in xrange(self._nhalos): # self._nhalos + 1?
			# AHF halos index start from 0... python start from 0.. lets be consitant :P
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

	def _load_ahf_particles(self, f, snap):
		NotImplementedError("Only halo loading implemented")

	@staticmethod
	def _can_load(snap, **kwargs):
		tipsy_dir = str("%soutput_%05d_tipsy/" % (snap.path(), snap.output_number()))
		if not os.path.isdir(tipsy_dir):
			return False

		if len(glob.glob('%s*.AHF_particles'%tipsy_dir)) == 0:
			return False
		
		fname = glob.glob('%s*.AHF_halos'%tipsy_dir)[0]
		if os.path.exists(fname):
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
