'''
Based on Simulation.py from the Hamu project https://github.com/samgeen/Hamu
In the future I'd like all I/O to be done by an ObjectId, but that may require
moving to a database like MongoDB (Easily done with PyMongo) but this would be
quite a drain on resources....

TODO: Add features similar to Hamu i.e Automatic generation of axis labels for plots etc

@author: dsullivan, bthompson
'''
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

def init(name):
	'''
	Create a simulation from the current directory
	'''
	path = os.getcwd()
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

		self._name = name
		self._path = path

	def set_name(self, name):
		self._name = name

	def set_path(self, path):
		self._path = path

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
			return YT.load(self._path, ioutput, **kwargs)
		elif (module == 'pymses') and pymses_loaded:
			from .pymses import Pymses
			return Pymses.load(self._path, ioutput, **kwargs)
		elif (module == 'pynbody') and pynbody_loaded:
			from .pynbody import Pynbody
			return Pynbody.load(self._path, ioutput, **kwargs)
		else:
			print 'yt loaded: ', yt_loaded
			print 'pymses loaded: ', pymses_loaded
			print 'pynbody loaded: ', pynbody_loaded
			raise Exception("Unknown module: %s or not loaded"%module)

	def merger_tree(self, finder='rockstar'):
		'''
		Load a generic merger tree - default to rockstar if not overridden
		'''
		from ..analysis.halo_analysis import trees
		if finder == 'rockstar':
			return trees.RockstarMergerTree(self)
		else:
			raise Exception("Unimplemented finder: %s"%finder)

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

		return np.array(redshifts)

	def ordered_outputs(self):
		'''
		Return an ordered list of outputs
		'''
		from .utils import string_utils
		outputs = glob.glob('%s/output_00*'%self.path())
		outputs.sort(key=string_utils.natural_keys)
		return outputs

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
