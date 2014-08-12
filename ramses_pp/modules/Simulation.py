'''
Based on Simulation.py from the Hamu project https://github.com/samgeen/Hamu
In the future I'd like all I/O to be done by an ObjectId, but that may require
moving to a database like MongoDB (Easily done with PyMongo) but this would be
quite a drain on resources....

TODO: Add features similar to Hamu i.e Automatic generation of axis labels for plots etc

@author: dsullivan, bthompson
'''
from __future__ import division
from ramses_pp import config


pymses_loaded = config.pymses_enabled
pynbody_loaded = config.pynbody_enabled
yt_loaded = config.yt_enabled

if config.quick_import:
	if config.pymses_enabled:
		from .pymses import Pymses
	if config.pynbody_enabled:
		from .pynbody import Pynbody
	if config.yt_enabled:
		from .yt import YT
else:

	config.list_modules()

	if pymses_loaded:
		try:
			from .pymses import Pymses
		except ImportError as e:
			print 'Unable to import pymses'
			pymses_loaded = False
			print e
	if pynbody_loaded:
		try:
			from .pynbody import Pynbody
		except ImportError as e:
			print 'Unable to import pynbody'
			pynbody_loaded = False
			print e
	if yt_loaded:
		try:
			from .yt import YT
		except ImportError as e:
			print 'Unable to import yt'
			yt_loaded = False
			print e
		if (pymses_loaded or pynbody_loaded or yt_loaded) is False:
			raise RuntimeError("Could not import any modules!")
			
import numpy as np
import json, os, glob, uuid

from ramses_pp.applications import Halomaker

def load(name):
	'''
	Load the simulation with this name from the json dump
	'''
	json_dir = config.json_dir
	filename = '%s/%s.json'%(json_dir, name)
	if os.path.isfile(filename):
		with open(filename, 'rb') as fp:
			data = json.load(fp)
		return Simulation(str(data['_name']), str(data['_path']), str(data['_boxsize']), data['_halomaker_info'], data['_periodic'], data['_oid'])
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
	def __init__(self, name, path, boxsize=100, halomaker=None, periodic=True, oid=None):
		# This should never change
		if oid is None: self._oid = str(uuid.uuid4())
		else:
			print oid 
			self._oid = str(oid)

# calculate box size

		

		self._name = name
		self._path = path
		self._boxsize = self.box_size() #100   #100 #in Mpc h^-1
		self._halomaker_info = {  # store input parameters for HaloMaker
				'method': 'MSM',  
				'b': 0.2,
				'cdm' : ".false.",
				'npart' : '20',
				'adaptahop' : {
					'nvoisins' : 32,
					'nhop' : 16,
					'rhot' : 80,
					'fudge' : 4.0,
					'fudgepsilon' : 0.0,
					'alphap' : 3.0,
					'megaverbose' : ".false.",	
					} ,
				'verbose' : ".true.",
			} ,
		self._periodic = True,

	# add new methods based on "what modules (pymses, pynbody) are working

		if pynbody_loaded:
			from ..applications.ahf import AHF
			for f in dir(AHF):
				if f[0] != "_": # ignore hidden functions
					self.func = self.call(getattr(AHF,f))  # loads any arbitary function
			self.run_ahf = self.call(getattr(AHF,dir(AHF)[2]))
			self.run_ahf_merger = self.call(getattr(AHF,dir(AHF)[3]))
			self.run_ahf_tracker = self.call(getattr(AHF,dir(AHF)[4]))

					

	def func(self,a):
		print("Not Defined")

	def call(self,func,*args,**kwargs):
		return lambda *args, **kwargs : func(self, *args, **kwargs)

	def set_name(self, name):
		self._name = name

	def set_path(self, path):
		self._path = path

#	def set_boxsize(self, boxsize):
#		if isinstance(boxsize, int):
#			self._boxsize = int(boxsize)  # may mess up with your simulation if this is incorrect
#		else:
#			print "Invalid boxsize, boxsize needs to be an integer"
#			return
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

	def num_snapshots(self):
		return len(glob.glob('%s/output_00*'%self._path))

	def save_image(self, prefix, ioutput, plt, **kwargs):
		'''
		Save an image to the simulation directory
		'''
		filename = '%s/%s_%05d.png'%(self.data_dir(), prefix, ioutput)
		if os.path.isfile(filename) and config.override is not True:
			raise Exception("File: %s already exists. Override is disabled"%filename)
		plt.savefig(filename)

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

	def snapshot(self, ioutput, module=config.default_module):
		'''
		Return a snapshot from the simulation
		'''
		if (module == 'yt') and yt_loaded:
			from .yt import YT
			return YT.load(self._path, ioutput)
		elif (module == 'pymses') and pymses_loaded:
			from .pymses import Pymses
			return Pymses.load(self._path, ioutput)
		elif (module == 'pynbody') and pynbody_loaded:
			from .pynbody import Pynbody
			return Pynbody.load(self._path, ioutput)
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

		h = h0 / 100

		first_snap = 1
		info = ("%s/output_%05d/info_%05d.txt" % (self._path, first_snap, first_snap))
		f = open(info, 'r')
		nline = 1  # read the last info file
		while nline <= 18:
			line = f.readline()
			if(nline == 10): iaexp = np.float32(line.split("=")[1])
			nline += 1

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
		Locate the snapshot closest to the given redshift, deprecated since it assumes you have all the snapshots in one location, replaced by redshift and avail_redshift
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

		return np.array(redshifts)

 	def ordered_outputs(self):
 		'''
 		Return an ordered list of outputs
 		'''

		from .utils import string_utils
 		outputs = glob.glob('%s/output_00*'%self.path())
 		outputs.sort(key=string_utils.natural_keys)
 		return outputs



### halomaker stuff

	def halomaker_info(self):
		return self._halomaker_info

	def halomaker(self, subvol=False, ncpu=1):
		return Halomaker.Halomaker(self, subvol=False, ncpu=1)
		


#### end halomaker stuff

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
