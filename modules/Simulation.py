'''
Based on Simulation.py from the Hamu project https://github.com/samgeen/Hamu
In the future I'd like all I/O to be done by an ObjectId, but that may require
moving to a database like MongoDB (Easily done with PyMongo) but this would be
quite a drain on resources....

TODO: Add features similar to Hamu i.e Automatic generation of axis labels for plots etc

@author: dsullivan, bthompson
'''
from __future__ import division
from .. import config
pymses_loaded = config.pymses_enabled
pynbody_loaded = config.pynbody_enabled
yt_loaded = config.yt_enabled
import numpy as np
import json, os, glob, uuid


def load(name):
	'''
	Load the simulation with this name from the json dump
	'''
	json_dir = config.json_dir
	filename = '%s/%s.json'%(json_dir, name)
	if os.path.isfile(filename):
		with open(filename, 'rb') as fp:
			data = json.load(fp)
		return Simulation(str(data['_name']), str(data['_path']), data['_oid'])
	else:
		raise Exception("No simulation with the name: %s"%name)


def init(name,path=None):
	'''
	Create a simulation from the current directory
	'''

	if path == None and config_path==False:
		path = os.getcwd()
	return create(name, path)

def new(name,rename=None):
	'''
	Create a simulation from the Simulation Directory
	'''
	path = config.simulation_dir + '/' + name
	if rename != None:
		name = rename # if you wish to call your simulation in the database by a different name
	return create(name, path)


def create(name, path):
	'''
	Create a Simulation object and write it to disk
	'''
	if os.path.isdir(path):
		json_dir = config.json_dir

		simulation = Simulation(name, path)
		simulation.save()

		data_dir = '%s/%s'%(json_dir, simulation._name)
		if not os.path.isdir(data_dir): os.mkdir(data_dir)
		return simulation
	else:
		raise Exception("Path does not exist: %s"%path)

def list():
	'''
	List all available simulations
	'''
	json_dir = config.json_dir

	for f in glob.glob('%s/*.json'%json_dir):
		print f


class Simulation():
	def __init__(self, name, path, oid=None):
		# This should never change
		if oid is None: self._oid = str(uuid.uuid4())
		else:
			print oid 
			self._oid = str(oid)

# calculate box size

		

		self._name = name
		self._path = path
		self._boxsize = self.box_size() #100   #100 #in Mpc h^-1

		# load more data
		json_dir = config.json_dir
		filename = '%s/%s.json'%(json_dir, name)
		with open(filename, 'rb') as fp:
			data = json.load(fp)
		self._periodic = True
		if "_periodic" in data:
			if data["_periodic"] is not None:
				self._periodic = data["_periodic"]
		self._parent = None # assign a parent simulation
		if "_parent" in data:
			if data["_parent"] is not None:
				self._periodic = data["_parent"]

		self._parent_domain = [0.50,0.50,0.50] #coordinate in parent simulation which is this simulations centre in code units 
		if "_parent_domain" in data:
			if data["_parent_domain"] != None:
				self._parent_domain = data["_parent_domain"]

	# add new methods based on "what modules (pymses, pynbody) are working

		if pynbody_loaded:
			from ramses_pp.analysis.halo_analysis.ahf import AHF
			for f in dir(AHF):
				if f[0] != "_": # ignore hidden functions
					self.func = self.call(getattr(AHF,f))  # loads any arbitary function
			self.run_ahf = self.call(getattr(AHF,dir(AHF)[2]))
			self.run_ahf_merger = self.call(getattr(AHF,dir(AHF)[3]))
			self.run_ahf_tracker = self.call(getattr(AHF,dir(AHF)[4]))

					
	# 
	def func(self,a):
		print("Not Defined")

	def call(self,func,*args,**kwargs):
		return lambda *args, **kwargs : func(self, *args, **kwargs)

	# end custom method calling

	def set_name(self, name):
		self._name = name

	def set_path(self, path):
		self._path = path

	def box_size(self):
		cmtokpc = 3.24077929e-22
		kpctompc = 0.001
		last_snap = self.num_snapshots()
		info = ("%s/output_%05d/info_%05d.txt" % (self._path, last_snap, last_snap))
		f = open(info, 'r')
		nline = 1  # read the last info file
		while nline <= 18:
			line = f.readline()
			if(nline == 11): h0 = np.float32(line.split("=")[1])
			if(nline == 16): lunit = np.float32(line.split("=")[1])
			nline += 1
		h = h0 / 100
		boxsize = lunit * cmtokpc * kpctompc * h
		return boxsize
		




	def set_periodic(self, periodic):
		if isinstance(boxsize, bool):
			self._periodic = periodic
		else:
			print "Invalid input, True or False"
			return

	def is_periodic(self):
		return self._periodic
		

	def jdefault(self, o):
		if isinstance(o, set):
			return list(o)
		return o.__dict__

	def to_JSON(self):
		return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

	def save(self):
		json = self.to_JSON()
		json_dir = config.json_dir
		filename = '%s/%s.json'%(json_dir, self._name)
		with open(filename, 'w') as fi:
			fi.write(json)
		fi.close()

	def delete(self, rm_dir=False):
		'''
		Delete the JSON record for this simulation
		'''
		json_dir = config.json_dir
		filename = os.path.join('%s/'%json_dir, self._name)
		os.remove(filename)

		if rm_dir:
			data_dir = self.data_dir()
			os.removedirs(data_dir)
		return True

	def exists(self, fname):
		filename = '%s/%s'%(self.data_dir(), fname)
		return os.path.isfile(filename)

	def list(self, string=None):
		if string:
			return glob.glob('%s/%s'%(self.data_dir(), string))
		else:
			return glob.glob('%s/*'%self.data_dir())



	def assign_parent_domain(self,domain):
		if len(domain) == 3 and sum(domain) <= 3.0:
			self._parent_domain = domain
			self.save()
			return
		# else
		print "invalid domain"
		return


	def assign_parent(self, name, domain):
		'''
		for most zoom runs, you will have based a simulation off it's parent. Usually it is a unigrid simulation
		this function will store
		domain = 3 vector describing the position shift
		name = name of parent simulation
		'''
		# check if it exists first
		json_dir = config.json_dir
		filename = '%s/%s.json'%(json_dir, name)
		if os.path.isfile(filename):
			self._parent = name
			self.assign_parent_domain(domain)
			print name, domain
			self.save()
			return	
		else:
			raise Exception("No simulation with the name: %s"%name)
			return

	def num_snapshots(self):
		return len(glob.glob('%s/output_00*'%self._path))

	def save_image(self, prefix, ioutput, plt, **kwargs):
		'''
		Save an image to the simulation directory
		'''
		filename = '%s/%s_%05d.png'%(self.data_dir(), prefix, ioutput)
		if os.path.isfile(filename) and config.override is not True:
			raise Exception("File: %s already exists. Override is disabled"%filename)
		plt.savefig(filename, **kwargs)

	def write_array(self, prefix, ioutput, array, **kwargs):
		'''
		Save a numpy array to disk under this simulation
		'''
		filename = '%s/%s_%05d.out'%(self.data_dir(), prefix, ioutput)
		if os.path.isfile(filename) and config.override is not True:
			raise Exception("File: %s already exists. Override is disabled"%filename)
		np.savetxt(filename, array, **kwargs)
		return True

	def load_array(self, prefix, ioutput, **kwargs):
		'''
		Load a numpy array from the reserved data directory
		'''
		filename = '%s/%s_%05d.out'%(self.data_dir(), prefix, ioutput)
		return np.loadtxt(filename, **kwargs)

	def data_dir(self):
		'''
		Return the data directory reserved for this simulation
		'''
		json_dir = config.json_dir
		return os.path.join('%s/'%json_dir, '%s/'%self._name)

	def name(self):
		'''
		Return the name of this simulation
		'''
		return self._name

	def path(self):
		'''
		Return the path to this simulation
		'''
		return self._path

	def snapshot(self, ioutput, module=config.default_module, **kwargs):
		'''
		Return a snapshot from the simulation
		'''
		if (module == 'yt') and yt_loaded:
			from .yt import YT
			return YT.load(self._path, self, ioutput=ioutput, **kwargs)
		elif (module == 'pymses') and pymses_loaded:
			from .pymses import Pymses
			return Pymses.load(self._path, self, ioutput, **kwargs)
		elif (module == 'pynbody') and pynbody_loaded:
			from .pynbody import Pynbody
			return Pynbody.load(self._path, self, ioutput, **kwargs)
		else:
			print 'yt loaded: ', yt_loaded
			print 'pymses loaded: ', pymses_loaded
			print 'pynbody loaded: ', pynbody_loaded
			raise Exception("Unknown module: %s or not loaded"%module)


	def initial_conditions(self):
		'''
		Returns the cosmology at z=0
		'''
		last_snap = self.num_snapshots()
		info = ("%s/output_%05d/info_%05d.txt" % (self._path, last_snap, last_snap))
		f = open(info, 'r')
		nline = 1  # read the last info file
		while nline <= 18:
			line = f.readline()
			if(nline == 10): faexp = np.float32(line.split("=")[1])
			if(nline == 11): h0 = np.float32(line.split("=")[1])
			if(nline == 12): omega_m_0 = np.float32(line.split("=")[1])
			if(nline == 13): omega_l_0 = np.float32(line.split("=")[1])
			if(nline == 14): omega_k_0 = np.float32(line.split("=")[1])
			if(nline == 15): omega_b_0 = np.float32(line.split("=")[1])
			if(nline == 16): lunit = np.float32(line.split("=")[1])
			if(nline == 17): dunit = np.float32(line.split("=")[1])
			if(nline == 18): tunit = np.float32(line.split("=")[1])
			nline += 1
		f.close()
		h = h0 / 100

		first_snap = 1
		info = ("%s/output_%05d/info_%05d.txt" % (self._path, first_snap, first_snap))
		f = open(info, 'r')
		nline = 1  # read the last info file
		while nline <= 18:
			line = f.readline()
			if(nline == 10): iaexp = np.float32(line.split("=")[1])
			nline += 1
		f.close()
		# store variables into a dictionary

		iz = 1.0/iaexp -1.0
		fz = 1.0/faexp -1.0

		initial_cons = {
			'iaexp':iaexp,		# initial a
			'faexp':faexp,		# final a
			'iz':iz, 			# initial z
			'fz':fz, 			# final z
			'omega_m_0':omega_m_0,
			'omega_l_0':omega_l_0,
			'omega_k_0':omega_k_0,
			'omega_b_0':omega_b_0,
			'h':h,
			'H0':h0,
			'lunit':lunit,
			'dunit':dunit,
			'tunit':tunit,
		}
		
		return initial_cons

	def merger_tree(self, finder=config.default_finder):
		'''
		Load a generic mergertree. Default is set in config if not overwritten
		'''

		from ..analysis.halo_analysis import trees
		if finder == 'rockstar':
			return trees.RockstarMergerTree(self)
		else:
			raise Exception("Unimplemented finder: %s" % finder)


	def redshift(self, z):
		'''
		Locate the snapshot closest to the given redshift
		'''
		from .utils import array_utils
		redshifts = self.avail_redshifts()
		idx = array_utils.argmin(redshifts, z)
		if config.verbose: print 'ioutput %05d closest to redshift %f'%(idx+1, z)

		return idx+1

	def avail_redshift(self, z):
		'''
		Return the closest available redshift
		'''
		redshifts = self.avail_redshifts()
		idx = self.redshift(z)
		return redshifts[idx]


	def redshift_deprecated(self, z):
		'''
		Locate the snapshot closest to the given redshift
		'''
		#First, gather list of redshifts
		num_snapshots = self.num_snapshots()
		redshifts = np.zeros(num_snapshots)

		i = 0
		for ioutput in range(1, num_snapshots+1):
			info = ("%s/output_%05d/info_%05d.txt" % (self._path, ioutput, ioutput))
			f = open(info, 'r')
			nline = 1
			while nline <= 10:
				line = f.readline()
				if(nline == 10): aexp = np.float32(line.split("=")[1])
				nline += 1
			redshift = 1.0/aexp -1.0
			redshifts[i] = float(redshift)
			i += 1
			f.close()

		idx = np.argmin(np.abs(redshifts - z))
		if config.verbose: print 'ioutput %05d closest to redshift %f'%(idx+1, z)

		return idx+1

	def avail_redshifts(self, zmin=None, zmax=1000):
		'''
		Return a list of the available redshifts between some range
		'''
		redshifts = []
		outputs = self.ordered_outputs()
		for output in outputs:
			info = '%s/info_%s.txt'%(output, output[-5:])
			f = open(info, 'r')
			nline = 1
			while nline <= 10:
				line = f.readline()
				if(nline == 10): aexp = np.float32(line.split("=")[1])
				nline += 1
			redshift = 1.0/aexp -1.0
			if (redshift >= zmin) and (redshift < zmax):
				redshifts.append(float(redshift))
			f.close()

		return np.array(redshifts)

 	def ordered_outputs(self):
 		'''
 		Return an ordered list of outputs
 		'''

		from .utils import string_utils
 		outputs = glob.glob('%s/output_00*'%self.path())
 		outputs.sort(key=string_utils.natural_keys)
 		return outputs


	def pos_to_codepos(self,val):
		''' simple dirty code to convert a position in Mpc/h into code units'''
		if val > self.box_size():
			print "value is larger than the box size"
			print "the boxsize is " + str(self.box_size())
			return
		else:
			new_val = val / self.box_size()
			print str(val) + " Mpc/h in code units is " + str(new_val)
			return new_val



### halomaker stuff

	def load_halomaker(self, subvol=False, ncpu=1):
		'''
		Returns the Halomaker Application Object
		'''
		from ramses_pp.analysis.halo_analysis.HaloMaker import HaloMaker
		return HaloMaker(self, subvol=False, ncpu=1)



#### end halomaker stuff

	def info_snap(self,ioutput):
		'''
		return the basic infomation of an individual snapshot without loading it
		'''
		num = self.num_snapshots()
		if ioutput > num or ioutput < 1:
			print "Snapshot needs to be in range of 1 and " + str(self.num_snapshots())
			raise e
			return

		infopath = ("%s/output_%05d" % (self._path, ioutput))
		if not os.path.isdir(infopath):
			print "Snapshot does not exist"
			raise e
			return
		
		info = infopath + ("/info_%05d.txt" % (ioutput))
		f = open(info,'r')
		nline = 1
		while nline <= 18:
			line = f.readline()
			if(nline == 10): aexp = np.float32(line.split("=")[1])
			if(nline == 16): lunit = np.float32(line.split("=")[1])
			if(nline == 17): dunit = np.float32(line.split("=")[1])
			if(nline == 18): tunit = np.float32(line.split("=")[1])
			nline += 1
		z = (1.0/aexp) - 1.0

		f.close()
		infodata = {
			'aexp': aexp,
			'lunit': lunit,
			'dunit': dunit,
			'z':z,
			}
		return infodata
	
	def info(self):
		'''
		List all snapshots with some basic info
		'''
		num = self.num_snapshots()

		print "## Output        aexp           z        unit_l        unit_d        unit_t"
		for ioutput in range(1, num+1):
			info = ("%s/output_%05d/info_%05d.txt" % (self._path, ioutput, ioutput))
			f = open(info, 'r')
			nline = 1
			while nline <= 18:
				line = f.readline()
				if(nline == 10): aexp = np.float32(line.split("=")[1])
				if(nline == 16): lunit = np.float32(line.split("=")[1])
				if(nline == 17): dunit = np.float32(line.split("=")[1])
				if(nline == 18): tunit = np.float32(line.split("=")[1])
				nline += 1
			z = (1.0/aexp) - 1.0
			print "    %5d  %10.5f  %10.5f  %12.5e  %12.5e  %12.5e" % (ioutput, aexp, z, lunit, dunit, tunit)
		f.close()
