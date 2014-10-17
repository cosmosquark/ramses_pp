'''
heavily based on halos.py from pynbody and catalogUtil.py within VIDE, so credit to those guys
Rewritten to allow loading of Rockstar catalogues when using any module

Currently depends on yt-3.0 for unit handling. I'll try and code this
dependency out sometime in the future if needed.

@author bthompson
=======
'''

super_verbose=True


import yt
import numpy as np
import sys, os.path, glob
import weakref, copy
import re, uuid
from ramses_pp import config as pp_cfg
import abc
import pickle
from netCDF4 import Dataset

import videtools
import shelve

NetCDFFile = Dataset
ncFloat = 'f8'

class DummyVoid(object):
	def __init__(self):
		self.properties = {}

class Void():

	def __init__(self, void_id, void_catalogue, *args):
		self.properties = {}
		self._void_catalogue = void_catalogue
		self._void_id = void_id
		self._descriptor = "void_" + str(void_id)

	def __getitem__(self, item):
		snap = self._void_catalogue._base()
		ds = snap.raw_snapshot()
		unit = self._void_catalogue.units[item]
		return ds.arr(self.properties[item], unit)

	def field_list(self):
		return self._halo_catalogue.void_type

	


class VoidCatalogue():
	'''
	Generic Void Catalogue
	'''


	def __init__(self,finder,path,ioutput):
		self._voids = {}
		self._finder = finder
		self._simpath = path
		self._ioutput = ioutput

	def calc_item(self, i):
		if i in self._voids:
			return self._voids[i]
		else:
			v = self._get_void(i)
			self._voids[i] = v
			return v

	def __len__(self):
		return len(self._voids)

	def __iter__(self):
		return self._void_generator()

	def __getitem__(self,item):
		if isinstance(item, slice):
			indices = item.indices(len(self._voids))
			[self.calc_item(i) for i in range(*indices)]
			return self._voids[item]
		else:
			return self.calc_item(item)

	def _void_generator(self):
		i = 0
		while True:
			try:
				yield self[i]
				i += 1
				if i > len(self._voids)-1:
					break
			except RuntimeError:
				break

	def __contains__(self,voidid):
		return self.contains(voidid)

	@staticmethod
	def can_load(self):
		return False

	@staticmethod
	def can_run(self):
		return False

	def sort(self, field, reverse=True):
		'''
		Sort halos by a given field
		'''
		return sorted(self, key=lambda x: x[field], reverse=reverse)


class VIDECatalogue(VoidCatalogue):
	'''
	Class to handle void catalogues, adapted from the Rockstar catalogue reader (by Peter Behroozi) and VIDE (by Paul Sutter et al)
	'''

	void_type = np.dtype([('id',np.int64),('FileVoidID',np.int64),('CoreParticle',np.int64),('CoreDens','f'),
				 ('ZoneVol','f'), ('Zone#Part',np.int64), ('Void#Zones',np.int64), ('VoidVol','f'),
				 ('Void#Part',np.int64),('VoidDensContrast','f'),('VoidProb','f'), ('VoidR','f'),
				 ('Reff','f'),('Volume','f'),('Macrocenter','f',3),('RA','f'),('Dec','f'),
				 ('Redshift','f'),('ParentID',np.int64),
				 ('TreeLevel',np.int64),('NumChildren',np.int64),('CentralDen','f'),
				 ('Ellipticity','f'),('EigenVals','f',3), 
				 ('EigenVecs','f',(3,3)) ])


        #cm = comoving
        #nocm = physical
        units = {'id':'dimensionless',
                        'FileVoidID':'dimensionless',
                        'CoreParticle':'dimensionless',
                        'CoreDens':'dimensionless',
                        'ZoneVol':'dimensionless',
                        'Zone#Part':'dimensionless',
                        'Void#Zones':'dimensionless',
                        'VoidVol':'dimensionless',
                        'Void#Part':'dimensionless',
                        'VoidDensContrast':'dimensonless',
                        'VoidProb':'dimensionless',
			'VoidR' : 'dimensionless',
                        'Reff':'Mpc/h',
			'Volume': 'Mpc ** 3 / h**3',
			'Macrocenter' : 'Mpc/h',
                        'RA':'dimensionless',
                        'Dec':'dimensionless',
			'Redshift':'dimensionless',
                        'ParentID':'dimensionless',
                        'TreeLevel':'dimensionless',
                        'NumChildren':'dimensionless',
                        'CentralDen':'dimensionless',
                        'Ellipticity':'dimensionless',
                        'EigenVals':'dimensionless',
                        'EigenVecs':'dimensionless',
			}

	def __init__(self, snap, particles=True, subsample=1.0, subvol="00", untrimmed=False, dataportion="central"):
#		if not self._can_load(snap):
#			raise Exception("Cannot locate/load VIDE catalogue")

		self._base = weakref.ref(snap)
		self._redshift = '%0.2f' % snap.cosmology()['z']
		VoidCatalogue.__init__(self,"VIDE",snap.path(),snap.output_number())
		if particles == True:
			self._voidtype = 'particles'
			self._prefix = 'sim_'
		else:
			self._voidtype = 'halos'
			self._prefix = 'sim_halos_'

		self._sim_name = os.path.split(snap.path())[1]
		self._short_name = str(self._prefix) + self._sim_name + "ss" + str(subsample)
		self._full_name = self._short_name + "_z" + str(self._redshift) + "_d" + subvol
		self._catalogue_path = snap.path() + "/vide/voids/" + self._short_name + "/sample_" + self._full_name + "/"
		print "Loading void info from " + str(self._catalogue_path) + "...."
#		sample =  shelve.open(self._catalogue_path+"sample_info.dat")
		f = open(self._catalogue_path+"sample_info.txt", 'r')
		nline = 1
		while nline <= 16:
			line = f.readline()
			if(nline == 3): self._datatype = str(line.split(":")[1])  # Data type
			if(nline == 4): self._zrange = str(line.split(":")[1])  # Redshift range
			if(nline == 5): self._islightcone = np.int32(line.split(":")[1])  # Particles placed on lightcone
			if(nline == 6): self._ispeculiarvel = np.int32(line.split(":")[1])  # Peculiar velocities included
			if(nline == 7): self._addsubsampfract = np.int32(line.split(":")[1])  # Additional subsampling fraction
			if(nline == 10): self._numsubvols = np.int32(line.split(":")[1])  # Number of simulation subvolumes
			if(nline == 11): self._subvolind = str(line.split(":")[1])  # My subvolume index
			if(nline == 12): self._estvol =	np.float32(line.split(":")[1]) # Estimated volume (cubic Mpc/h)
			if(nline == 13): self._realtracers = np.int32(line.split(":")[1]) # Number of real (non-boundary) tracers
			if(nline == 14): self._notracers = np.int32(line.split(":")[1]) # Total number of tracers
			if(nline == 15): self._esttracersep = np.float32(line.split(":")[1]) # in mpc/h
			if(nline == 16): self._minvoidsize = np.float32(line.split(":")[1]) # in mpc/h
			nline = nline + 1

#    			sample = pickle.load(input)
#		self.raw_sampleInfo = sample
		
		infoFile = self._catalogue_path + "zobov_slice_" + self._full_name + ".par"
		print "accessing " + str(infoFile) + "..."		

		File = NetCDFFile(infoFile, 'r')
		ranges = np.zeros((3,2))
		boxLen = np.zeros(3)
		
		ranges[0][0] = getattr(File, 'range_x_min')
		ranges[0][1] = getattr(File, 'range_x_max')
		ranges[1][0] = getattr(File, 'range_y_min')
		ranges[1][1] = getattr(File, 'range_y_max')
		ranges[2][0] = getattr(File, 'range_z_min')
		ranges[2][1] = getattr(File, 'range_z_max')
		self.isObservation = getattr(File, 'is_observation')
		self.maskIndex = getattr(File, 'mask_index')

		File.close()

		## store the boxlen infomation

 		boxLen[0] = ranges[0][1] - ranges[0][0]
  		boxLen[1] = ranges[1][1] - ranges[1][0]
  		boxLen[2] = ranges[2][1] - ranges[2][0]
  		self.ranges = ranges
		self.boxlen = boxLen

		volNorm = videtools.getVolNorm(self)
		self.volNorm = volNorm[0]


		if untrimmed:
			prefix = "untrimmed_"
		else:
			prefix = ""

		self._nvoids = 0 # will update this soon
		self._load_vide_voids(self._catalogue_path, prefix, dataportion, self._full_name)
		


	def _get_void(self, i):
		if self.base is None:
			raise RuntimeError("Parent snapshot has been deleted")
		return self._voids[i]


	

	def _load_vide_voids(self, catalogue_path, prefix, dataportion, full_name):
		fileName = catalogue_path + "/" + prefix + "voidDesc_" + dataportion + "_" + full_name + ".out"
		desc = np.loadtxt(fileName, comments="#", skiprows=2)

		fileName = catalogue_path + "/" + prefix + "macrocenters_" + dataportion + "_" + full_name + ".out"
		macrocenters = np.loadtxt(fileName)

		fileName = catalogue_path + "/" + prefix + "sky_positions_" + dataportion + "_" + full_name + ".out"
		sky = np.loadtxt(fileName, comments="#")

		fileName = catalogue_path + "/" + prefix + "shapes_" + dataportion + "_" + full_name + ".out"
		shapes = np.loadtxt(fileName, comments="#")

		fileName = catalogue_path + "/" + prefix + "centers_" + dataportion + "_" + full_name + ".out"
		derived = np.loadtxt(fileName, comments="#", usecols = (4,6,10,11,12,13))  # why are we getting volume and radius again (col 4 and 6 are the radius and volume respectively.. but only to 2 decimal places!!!!
		self._nvoids = len(desc)

		# lets merge this infomation together
		self._voidprops = np.array(np.zeros(int(self._nvoids),dtype=self.void_type))

		for vn in range(0,self._nvoids):

			self._voidprops[vn][0] = desc[vn][0]
			self._voidprops[vn][1] = desc[vn][1]
			self._voidprops[vn][2] = desc[vn][2]
			self._voidprops[vn][3] = desc[vn][3]
			self._voidprops[vn][4] = desc[vn][4]
			self._voidprops[vn][5] = desc[vn][5]
			self._voidprops[vn][6] = desc[vn][6]
			self._voidprops[vn][7] = desc[vn][7]
			self._voidprops[vn][8] = desc[vn][8]
			self._voidprops[vn][9] = desc[vn][9]
			self._voidprops[vn][10] = desc[vn][10]

		# will leave self._voidprops[11,:] for a sec (i.e voidR)
			self._voidprops[vn][14][:] = macrocenters[vn][1:]
			self._voidprops[vn][15] = sky[vn][0] # RA
			self._voidprops[vn][16] = sky[vn][1] # dec
			self._voidprops[vn][17] = sky[vn][2] # z
			
			self._voidprops[vn][12] = derived[vn][0]
			self._voidprops[vn][13] = derived[vn][1]
			self._voidprops[vn][18] = derived[vn][2]
			self._voidprops[vn][19] = derived[vn][3]
			self._voidprops[vn][20] = derived[vn][4]
			self._voidprops[vn][21] = derived[vn][5]
			self._voidprops[vn][22] = shapes[vn][1]
			self._voidprops[vn][23][0] = shapes[vn][2]
			self._voidprops[vn][23][1] = shapes[vn][3]
			self._voidprops[vn][23][2] = shapes[vn][4]
	
			self._voidprops[vn][24][0][0] = shapes[vn][5]
			self._voidprops[vn][24][0][1] = shapes[vn][6]
			self._voidprops[vn][24][0][2] = shapes[vn][7]
	
			self._voidprops[vn][24][1][0] = shapes[vn][8]
			self._voidprops[vn][24][1][1] = shapes[vn][9]
			self._voidprops[vn][24][1][2] = shapes[vn][10]
	
			self._voidprops[vn][24][2][0] = shapes[vn][11]
			self._voidprops[vn][24][2][1] = shapes[vn][12]
			self._voidprops[vn][24][2][2] = shapes[vn][13]

		
			self._voidprops[vn][11] = np.power((desc[vn][7] / self.volNorm * 3. / 4. / np.pi),1./3.)
			self._voids[vn] = Void(vn,self)
			self._voids[vn].properties = self._voidprops[vn]

#		self.desc = desc
#		self.macrocenters = macrocenters
#		self.sky = sky
#		self.shape = shape

		

		 
		

		
