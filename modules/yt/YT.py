'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

from .. import Snapshot
import sys, os
from ramses_pp import config as cfg 

#from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
import yt
#from yt import derived_field
#from yt.units.yt_array import YTQuantity
from yt.mods import *
from yt.analysis_modules.halo_finding.api import *
from yt.utilities.physical_constants import \
	boltzmann_constant_cgs, \
	mass_hydrogen_cgs, \
	mass_sun_cgs, \
	mh, \
	G

from yt.data_objects.particle_filters import add_particle_filter

verbose = True

#def rho_crit_now(data):
#	#G = 6.6743E-8 # cgs
#	H0 = YTQuantity((data.pf['H0'] * 1000) / 3.08567758E22, 's**-1') # s^-1
#	rho_crit_0 = (3 * H0**2) / (8 * np.pi * G.in_cgs())
#	return rho_crit_0


#Add some useful fields
def _Temperature(field, data):
	rv = data["pressure"]/data["density"]
	rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
	return rv

def _number_density(field, data):
	rv = data["density"]/mass_hydrogen_cgs
	return rv

def _over_density(field, data):
	omega_baryon_now = data.ds.parameters['omega_b']
	#return data['Density'] / (omega_baryon_now * rho_crit_now * (data.pf.hubble_constant**2) * ((1+data.pf.current_redshift)**3))
	return data['density'] / (omega_baryon_now * rho_crit_now(data) * ((1+data.ds.current_redshift)**3))

def load(path, ioutput, **kwargs):
	return YTSnapshot(path, ioutput, **kwargs)

def star_filter(pfilter,data):
	filter = np.logical_and(data["particle_age"] != 0, data["particle_age"] != None)
	return filter

def dark_filter(pfilter,data):
	filter = np.logical_and(data["particle_age"] == 0, data["particle_age"] != None)
	return filter

def young_star_filter(pfilter,data):
	filter = np.logical_and(data["particle_age"] != 0, data["particle_birth_epoch"].in_units("Gyr") <= data["particle_birth_epoch"].in_units("Gyr").min() + data.ds.arr("200","Myr"))
	return filter


class YTSnapshot(Snapshot.Snapshot):
	def __init__(self, path, ioutput, **kwargs):
		Snapshot.Snapshot.__init__(self, path, ioutput, "yt", **kwargs)
		'''
		Load the snapshot
		'''

		#self._snapshot = yt.load(os.path.join('%s/output_%05d/info_%05d.txt'%(path, ioutput, ioutput)))
		if "patch" in kwargs:
                        patch = kwargs.get("patch","default")
			print patch
		
		try:
			stars = kwargs.get("stars",False)
		except KeyError:
			stars = False

		try:
			dark = kwargs.get("dark",False)
		except KeyError:
			dark = False

		if "patch" in kwargs:
			self._snapshot = yt.mods.load(os.path.join('%s/output_%05d/info_%05d.txt'%(path, ioutput, ioutput)), patch=patch)
		else:
			self._snapshot = yt.mods.load(os.path.join('%s/output_%05d/info_%05d.txt'%(path, ioutput, ioutput)))

## snapshot filter methods ... for example, filtering out DM particles (creation time != 0)
		if stars == True:
			add_particle_filter("stars", function=star_filter, filtered_type="all", requires=["particle_age"])
			self._snapshot.add_particle_filter("stars")

		if dark == True:
			add_particle_filter("dark", function=dark_filter, filtered_type="all", requires=["particle_age"])
			self._snapshot.add_particle_filter("dark")
	#		add_particle_filter("stars", function=Stars, filtered_type='all', requires=["particle_type"])
			

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
		#There was an inconsitency here - but changing could break other code
		#return os.path.join(os.path.dirname(self._path), 'output_%05d/'%self._ioutput)
		return self._path


	def ncpu(self):
		'''
		Return the number of CPUs used
		'''
		ds = self.raw_snapshot()
		return ds.parameters['ncpu']

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

	def current_redshift(self):
		'''
		Return the current redshift
		'''
		ds = self.raw_snapshot()
		return ds.current_redshift
		

	def cosmology(self):
		'''
		Return an object with cosmological parameters
		'''
		ds = self.raw_snapshot()
		info = ds.parameters
		omega_m_0 = info['omega_m']
		omega_l_0 = info['omega_l']
		omega_k_0 = info['omega_k']
		omega_b_0 = info['omega_b']
		aexp = info['aexp']
		h  = info['H0']/100
		if aexp > 1.0:
			aexp = 1.0
		z = 1.0/aexp - 1.0

		cosmology = {
			'omega_m_0':omega_m_0,
			'omega_l_0':omega_l_0,
			'omega_k_0':omega_k_0,
			'omega_b_0':omega_b_0,
			'h':h,
			'aexp':aexp,
			'z':z,
		}
		return cosmology

	def info(self):
		'''
		Return info object
		'''
		return self.raw_snapshot().parameters

	def hop_path(self):
		hop_dir = os.path.join('%s/'%self._path, 'hop_dir')
		return hop_dir

	def rockstar_path(self):
		rockstar_dir = os.path.join('%s/'%self._path, cfg.rockstar_base)
		return rockstar_dir


	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return self._snapshot

	def raw_snapshot_all(self):
		'''
		Return the raw snapshot object
		'''
		return self._snapshot.all_data()

	def halos(self, finder=cfg.default_finder, **kwargs):
		from ramses_pp.analysis.halo_analysis import halos
		if finder=='rockstar':
			return halos.RockstarCatalogue(self)
		elif finder=='AHF':
			return halos.AHFCatalogue(self, **kwargs)
		elif finder == 'hop':
			return self.hop_halos()

		elif finder == "halomaker_simple":
			return halos.HaloMakerSimpleCatalogue(self)

	def hop_halos():
		ds = self._snapshot

		#Check if HOP file exists (note we will adopt a naming convention here)
		prefix = 'out_%05d_hop'%self._ioutput
		hop_dir = self.hop_path()

		#if not os.path.isdir(hop_dir): os.mkdir(hop_dir)

		if (verbose): print 'Loading halos from directory: %s'%hop_dir

		#Return the halos
		if self.halo_cat_exists() is False:
			dump_fname = '%s/%s'%(hop_dir, prefix)
			print 'No hop catalogue found, running hop...'
			print 'Will dump halos in: %s'%dump_fname
			halos = HaloFinder(ds, threshold=200)
			#Dump the halos for next time
			halos.dump(dump_fname)
 		else:
			halo_file = '%s/%s'%(hop_dir, prefix)
			if (verbose): print 'Loading halo file: %s'%halo_file
			halos = LoadHaloes(ds, halo_file)
			if (verbose): print 'Loaded %d halos'%len(halos)

		return halos

	def halo_cat_exists(self):
		#Check if HOP file exists (note we will adopt a naming convention here)
		hop_dir = self.hop_path()

		#if not os.path.isdir(hop_dir): os.mkdir(hop_dir)
		if (verbose): print 'Loading halos from directory: %s'%hop_dir

		#Check if out_%05d_hop.h5, .out and .txt exist
		prefix = 'out_%05d_hop'%self._ioutput
		extensions = ['h5', 'txt', 'out']
		exists = False

		for ext in extensions:
			if os.path.isfile('%s/%s.%s'%(hop_dir, prefix, ext)):
				exists = True
				break

		return exists

	def shift_parent(self, pos):
		domain=self._simulation._parent_domain
		if domain:
			x_min = self.raw_snapshot().arr([0.0,0.0,0.0],"code_length")
			x_max = self.raw_snapshot().arr([1.0,1.0,1.0],"code_length")
			x_cent = self.raw_snapshot().arr([0.5,0.5,0.5],"code_length")
#			shift = x_cent - pos
			domain = self.raw_snapshot().arr(domain,"code_length")
			units = pos.units # store the current units
			pos_new = pos.convert_to_units("code_length")
			x_max.convert_to_units("code_length")
##			print domain
			shift = x_cent - domain
##			print shift
			pos_new = pos_new + shift
			print pos_new, x_max
			# check for periodic boundries etc
			for i in range(0,2):
				thing = pos_new[i].value
				thing %= 1.0
				print thing
				pos_new[i] = thing
##			print pos_new
			pos_new = pos_new.convert_to_units(units) # convert back to origional units
			return pos_new

		else:
			print "No parent domain/center set, returning"
			return None


	def field_list(self):
		ds = self.raw_snapshot_all()
		for field in ds.derived_field_list:
			print field

	def trim_snapshot(self,size=1,units='Mpc',shape="sphere", region=[0.5,0.5,0.5]):
		''' by default, trim your simulation down to a 1 MPc sphere situated at the center '''
		raw_snap = self._snapshot
		if shape=="sphere":
			sp = raw_snap.sphere(region,(size,units))
			return sp
		else:
			print "invalid shape"
			return None

	def _a_dot(self):
		cos = self.cosmology()
		a_dot = (cos["h"]*100) * cos["aexp"] * np.sqrt(cos["omega_m_0"] * (cos["aexp"] ** -3) + cos["omega_k_0"]* (cos["aexp"] ** -2) + cos["omega_l_0"])
		return a_dot
