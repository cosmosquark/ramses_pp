'''
Based on Snapshot.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import abc
from .. import config
import numpy as np
import os
import glob

#Abstract snapshot class for loading a ramses output

class Snapshot():
	'''
	Abstract class for loading a single ramses snapshot. Pymses/Pynbody specific loaders are based on this class.
	'''
	__metaclass__ = abc.ABCMeta

	def __init__(self, path, ioutput, sim_type, **kwargs):
		'''
		Constructor
		'''
		self._path = path
		self._type = sim_type
		self._ioutput = ioutput
		self._attributes = {}
		self._info = None
		for key in kwargs.keys():
			self._attributes[key] = kwargs[key]

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


	@abc.abstractmethod
	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return

	def cube_positions(self, n):
		'''
		Return a list of cube positions for a cartesian
		domain decomp
		n - number of cubes per dimension
		'''
		ndim = self.ndim()
		boxlen = self.boxlen() # code units
		cubes = []
		dx = float(boxlen)/float(n)
		for i in range(0, n):
			#x coordinate of the cube centre
			cen_i = dx * (i + 0.5)
			if ndim > 1:
				for j in range(0, n):
					#y coordinate of the cube centre
					cen_j = dx * (j + 0.5)
					if ndim > 2:
						for k in range(0, n):
							# z coordinate of the cube centre
							cen_k = dx * (k + 0.5)
							cubes.append([cen_i, cen_j, cen_k])
					else:
						cubes.append([cen_i, cen_j])
			else:
				cubes.append([cen_i])
		return np.array(cubes), dx

	def num_snapshots(self):
		return len(glob.glob('%s/output_00*'%self._path))

	def has_key(self, key):
		'''
		Shortcut to self._attributes.has_key
		'''
		return self._attributes.has_key(key)

	def put(self, key, value):
		'''
		Add an item to the attributes dict
		'''
		self._attributes[key] = value

	def get(self, key):
		'''
		Return the key-value pair from the attributes dict
		'''
		return self._attributes[key]

	def info_path(self):
		'''
		Return the path to the info_000*.txt file
		'''
		return "%s/output_%05d/info_%05d.txt" % (self.path(), self.output_number(), self.output_number())

	def info_dict(self):
		'''
		Load info object using pymses
		'''
		from ramses_pp.modules.pymses import Pymses
		self._info = Pymses.info_dict(self)
		return self._info

	def cpu_list(self, bounding_box):
		'''
		Return the list of CPUs which cover the bounding box
		- bounding box: (2, ndim) ndarray containing min/max bounding box
		'''
		from ramses_pp.modules.pymses import Pymses
		dom_decomp = Pymses.hilbert_dom_decomp(self)
		return dom_decomp.map_box(bounding_box)

	def header_path(self):
		'''
		Return the path to the header_000*.txt file
		'''
		return "%s/output_%05d/header_%05d.txt" % (self.path(), self.output_number(), self.output_number())

	@abc.abstractmethod
	def ncpu(self):
		'''
		Return the number of CPUs used
		'''
		return

	@abc.abstractmethod
	def current_redshift(self):
		'''
		Return the current redshift
		'''
		return

	@abc.abstractmethod
	def cosmology(self):
		'''
		Return an object with cosmological parameters
		'''
		return

	@abc.abstractmethod
	def info(self):
		'''
		Return info object
		'''
		return

	def box_size(self):
		'''
		Return comoving boxsize in Mpc/h
		'''
		from .utils import constants as C
		with open(self.info_path(), 'r') as f:
			nline = 1  # read the last info file
			while nline <= 18:
				line = f.readline()
				if(nline == 11): h0 = np.float32(line.split("=")[1])
				if(nline == 16): lunit = np.float32(line.split("=")[1])
				nline += 1
		#f.close()
		h = h0 / 100
		boxsize = lunit * C.cmtokpc * C.kpctompc * h
		return boxsize

	def particles(self):
		'''
		Return total particle counts
		'''
		with open(self.header_path(), 'r') as f:
			nline = 1
			while nline <= 8:
				line = f.readline()
				if (nline == 2): npart = int(line)
				if (nline == 4): ndark = int(line)
				if (nline == 6): nstar = int(line)
				if (nline == 8): nsink = int(line)
				nline += 1

		#f.close()
		particles = {
			"total":npart,
			"dark":ndark,
			"star":nstar,
			"sink":nsink,
			}
		return particles
				

	def grid(self):
		'''
		Return level_min and level_max
		'''
		with open(self.info_path(), 'r') as f:
			nline = 1  # read the last info file
			while nline <= 18:
				line = f.readline()
				if(nline == 3): min = np.float32(line.split("=")[1])
				if(nline == 4): max = np.float32(line.split("=")[1])
				nline += 1

			grid = {
				"min":min,
				"max":max,
			}
#		f.close()
		return grid

	def is_zoom(self):
		'''
		Determines if this is a zoom simulation or not:
		A good proxy for a zoom run is the DM particle count, 
		since the result of a cube root follow by a log of base 2 will *NOT* be an integer
		'''

		particles = self.particles()
		dark = particles["dark"]
		# cube root
		dark_cube = np.float32(np.power(dark,(1.0/3.0)))
		dark_log = np.log2(dark_cube)
		if float(dark_log).is_integer(): # if an integer, then not a zoom
			return False
		else:
			return True
		

	def a_dot(self):
		cosmology = self.cosmology()
		omega_m_0 = cosmology['omega_m_0']
		omega_l_0 = cosmology['omega_l_0']
		omega_k_0 = cosmology['omega_k_0']
		aexp = cosmology['aexp']
		h0 = cosmology['h']*100
		a_dot = h0 * np.sqrt(omega_m_0 * (aexp ** -3) + omega_k_0 * (aexp ** -2) + omega_m_0)
		return a_dot

	def a_dot_recip(self):
		return 1. / self.a_dot()

	def hzoverh0(self):
		""" returns: H(a) / H0  = [omegam/a**3 + (1-omegam)]**0.5 """
		cosmology = self.cosmology()
		omega_m_0 = cosmology['omega_m_0']
		aexp = cosmology['aexp']
		return numpy.sqrt(omega_m_0 * numpy.power(aexp, -3) + (1. - omega_m_0))

	def _lingrowthintegrand(self):
		""" (e.g. eq. 8 in lukic et al. 2008)   returns: da / [a*H(a)/H0]**3 """
		cosmology = self.cosmology()
		return numpy.power((cosmology['aexp'] * self.hzoverh0()), -3)


	def _lingrowthfac(self, return_norm=False):
		"""
		returns: linear growth factor, b(a) normalized to 1 at z=0, good for flat lambda only
		a = 1/1+z
		b(a) = Delta(a) / Delta(a=1)   [ so that b(z=0) = 1 ]
		(and b(a) [Einstein de Sitter, omegam=1] = a)
		Delta(a) = 5 omegam / 2 H(a) / H(0) * integral[0:a] [da / [a H(a) H0]**3]
		equation  from  peebles 1980 (or e.g. eq. 8 in lukic et al. 2008) 

		"""

		cosmology = self.cosmology()
		omega_m_0 = cosmology['omega_m_0']
		omega_l_0 = cosmology['omega_l_0']
		aexp = cosmology['aexp']

		import scipy.integrate
		if (abs(omega_m_0 + omega_l_0 - 1.) > 1.e-4):
			raise RuntimeError, "Linear growth factors can only be calculated for flat cosmologies"
		
		lingrowth = scipy.integrate.quad(self._lingrowthintegrand(), 0., a, (omega_m_0))[0]
		lingrowth *= 5. / 2. * omega_m_0 * self.hzoverh0(a, omega_m_0)

		# then calc. for z=0 (for normalization)

		a0 = 1.0
		lingrowtha0 = scipy.integrate.quad( _lingrowthintegrand, 0., a0, (omega_m_0))[0]
		lingrowtha0 *= 5. / 2. * omega_m_0 * self.hzoverh0()
		
		lingrowthfactor = lingrowth / lingrowtha0
		if return_norm:
			return lingrowthfactor, lingrowtha0
		else:
			return lingrowthfactor

	def linear_growth_factor(self):
		"""Calculate the linear growth factor b(a), normalized to 1
		at z=0, for the cosmology of snapshot f.

		The output is dimensionless. If a redshift z is
		specified, it is used in place of the redshift in
		output f."""

		return self._lingrowthfac(return_norm=False)


	def tform(self, tform, rt=False):
		'''
		Return the time since the big bang when this particle formed IN Gyr.
		For lookback time, do simu.current_time - tform
		'''

		if hasattr(self, '_friedman') is False:
			self.integrate_friedman(store=True)

		cosmology = self.cosmology()
		h0 = cosmology['h']*100

		#These are all just pointers, so still memory efficient
		friedman = self._friedman
		axp_out = friedman['axp_out']
		#hexp_out = friedman['hexp_out']
		tau_out = friedman['tau_out']
		t_out = friedman['t_out']
		#age_tot = friedman['age_tot']
		#age_simu = friedman['age_simu']
		time_simu = friedman['time_simu']

		#unit_t = self.info()['unit_t']

		if rt:
			return (time_simu - tform)/(h0*1e5/3.08e24)/(365.*24.*3600.*1e9)

		ntable = len(axp_out)-1
		if config.verbose: print 'ntable = %d'%ntable

		output = np.zeros(len(tform))
		for j in range(0, len(tform)-1):

			i = 1
			while ( (tau_out[i] > tform[j]) and (i < ntable) ):
				i+=1

			#Interpolate time
			time = t_out[i] * (tform[j] - tau_out[i-1]) / (tau_out[i] - tau_out[i-1]) + \
				t_out[i-1] * (tform[j] - tau_out[i]) / (tau_out[i-1] - tau_out[i])

			time = max( (time_simu - time)/(h0*1e5/3.08e24)/(365*24*3600*1e9), 0 )
			output[j] = time

			#output[j] = (time_simu - time)*unit_t/(365 * 24 * 3600 * 1e9)

		return output

	def integrate_friedman(self, aexp=None, store=False):
		from ..fortran import friedman as fm

		cosmology = self.cosmology()
		omega_m_0 = cosmology['omega_m_0']
		omega_l_0 = cosmology['omega_l_0']
		omega_k_0 = cosmology['omega_k_0']
		if aexp == None: aexp = cosmology['aexp']
		h0 = cosmology['h']*100

		alpha=1e-6
		axpmin=1e-3
		ntable=1000

		axp_out, hexp_out, tau_out, t_out, age_tot = fm.friedman(omega_m_0, omega_l_0, omega_k_0, alpha, axpmin, ntable)

		#Find neighbouring expansion factor
		i = 1
		while ( (axp_out[i] > aexp) and (i < ntable) ):
			i+=1

		#Interpolate time
		time_simu = t_out[i] * (aexp - axp_out[i-1])/(axp_out[i]-axp_out[i-1]) + \
			t_out[i-1]*(aexp - axp_out[i])/(axp_out[i-1] - axp_out[i])

		age_simu = (time_simu+age_tot)/(h0*1e5/3.08e24)/(365*24*3600*1e9)

		friedman = {
				'axp_out':axp_out,
				'hexp_out':hexp_out,
				'tau_out':tau_out,
				't_out':t_out,
				'age_tot':age_tot,
				'age_simu':age_simu,
				'time_simu':time_simu
			}

		if store:
			self._friedman = friedman

		return friedman

	def halos(self, finder=config.default_finder,halo=None):
		'''
		Load a generic halo catalogue
		'''
		from ..analysis.halo_analysis import halos
		if finder=='rockstar':
			return halos.RockstarCatalogue(self)
		elif finder=="AHF":
			# if halo = an id, then track back the halo over time
			return halos.AHFCatalogue(self,halo=halo)
		elif finder=="halomaker_simple":
			return halos.HaloMakerSimpleCatalogue(self)
		else:
			raise Exception("Unimplemented finder: %s"%finder)


	def swap_modules(self, simulation, module):
		'''
		Swap to another swap_module
		'''
		return simulation.snapshot(self._ioutput, module=module)
		

	def voids(self, finder=config.void_finder,particles=True, subsample=1.0, subvol="00", untrimmed=False, dataportion="central", parent=True):
		'''
		Load a generic void catalogue
		'''

		from ..analysis.void_analysis import voids

		finder="VIDE"
		print finder
		if finder=="VIDE":
			# get parent simulation
			# fuck it.. just going to hack it in for now
			import Simulation
#			simparent = Simulation.load("void_sizes_n07_100")
#			snapparent = simparent.snapshot(simparent.num_snapshots())	
			voidcat =  voids.VIDECatalogue(self,particles=True, subsample=1.0, subvol="00", untrimmed=False, dataportion="central")
#			for v in range(0,len(voidcat)):
#				voidcat[v]["Macrocenter"] = self.shift_parent(voidcat[v]["Macrocenter"])
			return voidcat

		else:
			raise Exception("Unimplemented finder: %s"%finder)


	def vide_input(self,lightcone=False, pecval = False, threads=4, divisions=1, slices=1, subvolumes = 1, particles=True, subsamples = [1.0], halo_min_masses = ["none"], finder=None):
		'''
		Because VIDE is OpenMP and not MPI parallalised, I would not recommend running VIDE on a cluster.
		to get around this, run vide on a computer with a high amount of RAM.
		this code will generate inputs to VIDE for you which you can
	
		this code just generates VIDE inputs for RAMSES for a single snapshot... this also better doccuments the RAMSES modules within VIDE

		useful references
		http://cosmicvoids.net
		https://bitbucket.org/cosmicvoids/vide_public

		config.vide_catalogue_root

		input file is written to the simulation directory

		other things
		- lightcone = place particles on the lightcone (z-axis in sims)?
		- pecval = add peculiar velocities?
		- threads = optimization: maximum number of parallel threads to use
		- divisions = optimization: number of subdivisions of the box
		- slices = how many independent slices along the z-axis?
		- subvolumes = how many subdivisions along the x- and y- axis? ( = 2 will make 4 subvolumes for each slice, = 3 will make 9, etc.)
		- subsamples = list of desired subsamples - these are in unts of h Mpc^-3!. subsampling basically removes a random fraction of particles from VIDE (e.g if subSamples = 0.5, then subsampling will remove half of the particles.. this is also a list so you can run multiple subsamples (e.g [1.0,0.8,0.6])
		- particles ... if true.. run void analysis on particles... else do it on halos
		- halo_masses = a list of halo masses to consider (e.g ["none", 1.2e13]), if "none" then run the analysis for all halos without a mass cut


		NB: this routine is designed only for runs on an individual snapshot... if you want to make runs for all snapshots, there will eventually be a routine in the main simulation object

		I would leave subvolumes and slices as 1, but this has more use when your considering lightcones.. of which I need to do more with (~ben)

		Halos side is not tested yet.. the first line in the halo file should have the finder commented into it
		
		also, since VIDE has no concept of particle mass?.. I highly recommend not running VIDE on zoom runs
		'''

		if self.is_zoom():
			print "VIDE support for zoom sims is questionable"
			return None

		sim_path = self.path()
		ioutput = self.output_number()

		if not os.path.isdir('%s/vide'%(sim_path)):
			os.mkdir('%s/vide'%(sim_path))

		if not os.path.isdir('%s/vide/halos'%(sim_path)):
			os.mkdir('%s/vide/halos'%(sim_path))

		if particles:
			filename = os.path.join('%s/vide/vide_input_particles_%05d.py'%(sim_path,ioutput))
		else:

			# does the halo file exist? stored in simulation/vide_halos/output_NNNNN.txt
			if not os.path.isdir(os.path.join('%s/vide/halos'%(sim_path))):
				print "vide_halos directory does not exist"
				return None
			if finder == None:
				if not os.path.isfile('%s/vide/halos/output_%05d.txt' % (sim_path, ioutput)):
					print "vide halos input file does not exist"
					return None
				vide_halos_location = str('%s/vide/halos'%(sim_path)) + ("/output_%05d.txt" % ioutput)
				f.open(vide_halos,'r')
				line = f.readline()
				finder = str(line.split("=")[1])
				f.close()
				copy_file = os.path.join('%s/vide/vide_input_halos_%05d.py'%(sim_path,ioutput))
				copy_destination = os.path.join('%s/vide/vide_input_halos_%s_%05d.py'%(sim_path,finder,ioutput))
				os.system(("cp %s %s" % (copy_file, copy_destination)))
			else: # in most cases, we will have a finder
				if not os.path.isfile('%s/vide/halos/output_%s_%05d.txt' % (sim_path, finder, ioutput)):
					print "vide halos input file for " + str(finder) + " does not exist"
					return None

			vide_halos = "vide/halos/output_%s_NNNNN.txt" % str(finder)
			filename = os.path.join('%s/vide/vide_input_halos_%s_%05d.py'%(sim_path,str(finder),ioutput))


		sim_name = os.path.split(sim_path)[1] 
		vide_catalogDir = os.path.join('%s/%s/'%(config.vide_catalogue_root ,sim_name))

		# calculate redshift
		redshift = self.current_redshift()
		if redshift < 0.0:
			redshift = 0.0

		cosmology = self.cosmology()

		dark_part = self.particles()['dark']
		part_1d = str(int(np.rint(np.power(dark_part,(1.0/3.0)))))  # that rint is crucial otherwise int(256.0) = 255!


		f = open(filename,'w')
		f.write('# remember, this is free software.. GNU stuff.. no warranty... GNU licensing.. etc same as any other free software.. free as in free beer \n')
		f.write('# if you wish to know what any of these lines do, then look at the sample ramses file or sample files online https://bitbucket.org/cosmicvoids/vide_public \n')
	
		# configuration

		f.write('#CONFIGURATION \n')
		f.write('continueRun = False \n')
		f.write('startCatalogStage = 1 \n')
		f.write('endCatalogStage = 3 \n')
		f.write('catalogDir = \"%s\" \n' % vide_catalogDir)
		f.write('voidOutputDir = \"%s/vide/voids/\" \n' % vide_catalogDir)
		f.write('logDir = \"%s/vide/logs/\" \n' % vide_catalogDir)
		f.write('figDir = \"%s/vide/figs/\" \n' % vide_catalogDir)
		f.write('scriptDir = \"%s/vide/scripts/\" \n' % vide_catalogDir)
		f.write('dataType = \"simulation\" \n')
		f.write('dataFormat = \"ramses\" \n')
		f.write('dataUnit = 1\n') ## assuming you are using your data units between 0.0 and 1.0.. care
		if lightcone:
			f.write('useLightCone = True \n')
		else:
			f.write('useLightCone = False \n')
		if pecval:
			f.write('doPecVel = True \n')
		else:
			f.write('doPecVel = False \n')
		f.write('numZobovThreads = %01d \n' % threads)
		f.write('numZobovDivisions = %01d \n' % divisions)
		if particles:
			f.write('prefix = \"sim_%s\" \n' % sim_name)
		else:
			f.write('prefix = \"sim_halos_%s_%s\" \n' % (finder, sim_name))
		f.write('numSlices = %01d \n' % slices)
		f.write('numSubvolumes = %01d \n' % subvolumes)
		
		# particles, things below this should not need to be edited	

		if particles:
			f.write('#PARTICLES \n')
			f.write('particleFileBase = \"output_NNNNN\" \n')
			f.write('particleFileDummy = \"NNNNN\" \n')
			f.write('fileNums = [\"%05d\"] \n' % ioutput)	
			f.write('redshifts = [\"%s\"] \n' % str(redshift))
			f.write('subSamples = %s \n' % str(subsamples))
			f.write('doSubSampling = False \n')
			f.write('doSubSamplingInPrep = False \n')
			f.write('subSampleMode = \'relative\'\n')
			f.write('shiftSimZ = False \n')

		# find voids in halos
		else:
			f.write('#HALOS \n')
			f.write('#FINDER = \n' % finder)
			f.write('haloFileBase = %s \n' % vide_halos )
			f.write('haloFileDummy = \'NNNNN\' \n')
			f.write('minHaloMasses = %s \n' % halo_min_masses)
			f.write('haloFileMCol  = 6 \n')
			f.write('haloFileXCol  = 0 \n')
			f.write('haloFileYCol  = 1 \n')
			f.write('haloFileZCol  = 2 \n')
			f.write('haloFileVXCol = 3 \n')
			f.write('haloFileVYCol = 4 \n')
			f.write('haloFileVZCol = 5 \n')
			f.write('haloFileColSep = ',' \n')
			f.write('haloFileNumComLines = 1 \n') # we have left the finder name within the file

		# simulation infomation

		f.write('numPart = %s*%s*%s \n' % (part_1d,part_1d,part_1d) )
		f.write('lbox = %f \n' % self.box_size())
		f.write('omegaM = %s \n' % cosmology['omega_m_0'])
		f.write('hubble = %s \n' % cosmology['h'])
		f.close()
		print "VIDE input file saved as ", filename
		return

