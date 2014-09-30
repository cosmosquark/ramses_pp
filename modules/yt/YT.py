'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

from .. import Snapshot
import sys, os
#from ramses_pp import config

#from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
import yt
from yt.mods import *
from yt.analysis_modules.halo_finding.api import *
from yt.utilities.physical_constants import \
	boltzmann_constant_cgs, \
	mass_hydrogen_cgs, \
	mass_sun_cgs, \
	mh

from yt.data_objects.particle_filters import add_particle_filter

verbose = True
#simulation_dir = config.simulation_dir

def rho_crit_now(data, units='cgs'):
	if units == 'SI':
		G = 6.6743E-11
	else:
		G = 6.6743E-8 # cgs
	H0 = (data.ds['H0'] * 1000) / 3.08567758E22 # s^-1
	rho_crit_0 = (3 * H0**2) / (8 * np.pi * G)
	return rho_crit_0

#Add some useful fields
def _Temperature(field, data):
	rv = data["Pressure"]/data["Density"]
	rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
	return rv

def _NumDens(field, data):
	rv = data["Density"]/mass_hydrogen_cgs
	return rv

def _OverDensity(field, data):
	omega_baryon_now = data.ds.parameters['omega_b']
	#return data['Density'] / (omega_baryon_now * rho_crit_now * (data.pf.hubble_constant**2) * ((1+data.pf.current_redshift)**3))
	return data['Density'] / (omega_baryon_now * rho_crit_now(data) * ((1+data.ds.current_redshift)**3))

def load(folder, ioutput, **kwargs):
	add_field(("gas", "Temperature"), function=_Temperature, units=r"\rm{K}")
	add_field(("gas", "Number Density"), function=_NumDens, units=r"\rm{cm}^{-3}")
	add_field(("gas", "Baryon Overdensity"), function=_OverDensity,
          units=r"")
	return YTSnapshot(folder, ioutput, **kwargs)

def star_filter(pfilter,data):
	filter = np.logical_and(data.particles.source["particle_age"] != 0, data.particles.source["particle_age"] != None)
	return filter

def dark_filter(pfilter,data):
	filter = np.logical_and(data.particles.source["particle_age"] == 0, data.particles.source["particle_age"] != None)
	return filter

class YTSnapshot(Snapshot.Snapshot):
	def __init__(self, folder, ioutput, **kwargs):
		Snapshot.Snapshot.__init__(self, "yt")
		'''
		Load the snapshot
		'''
		print folder
		self._path = folder #simulation_dir + "/" + folder
		self._ioutput = ioutput
		self._snappath = os.path.join('%s/output_%05d/'%(self._path, ioutput))
		
		try:
			stars = kwargs.get("stars",False)
		except KeyError:
			stars = False

		try:
			dark = kwargs.get("dark",False)
		except KeyError:
			dark = False

		self._snapshot = yt.mods.load(os.path.join('%s/output_%05d/info_%05d.txt'%(folder, ioutput, ioutput)))

## snapshot filter methods ... for example, filtering out DM particles (creation time != 0)
		if stars == True:
			print "stars"
			add_particle_filter("stars", function=star_filter, filtered_type="all", requires=["particle_age"])
			self._snapshot.add_particle_filter("stars")

		if dark == True:
			print "dark"
			add_particle_filter("dark", function=dark_filter, filtered_type="all", requires=["particle_age"])
			self._snapshot.add_particle_filter("dark")
	#		add_particle_filter("stars", function=Stars, filtered_type='all', requires=["particle_type"])
			

		#Implement abstract methods from Snapshot.py

	def output_number(self):
		'''
		Return the output number for this snapshot
		'''
		return self._ioutput

	def folder(self):
		return self._folder

	def path(self):
		'''
		Return the path to this simulation directory
		'''
		return self._path

	def snappath(self):
		'''
		Return the path to this snapshot
		'''
		return self._snappath

	def ncpu(self):
		'''
		Return the number of CPUs used
		'''
		ds = self.raw_snapshot()
		return ds.parameters['ncpu']

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
		return self.raw_snapshot().parameters

	def hop_path(self):
		hop_dir = os.path.join('%s/'%self._path, 'hop_dir')
		return hop_dir

	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return self._snapshot

	def raw_snapshot_all(self):
		'''
		Return all the raw snapshot object data
		'''
		raw_snap = self.raw_snapshot()
		raw_snap_all = raw_snap.all_data()
		raw_snap_all.ds    # patch to get the pf functionality working nicely... bit of a hack
		return raw_snap_all

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



	#Return the HOP halo catalogue. Can override run_hop to force re-running
	def halos(self, finder="AHF", run_finder=False):
		ds = self._snapshot
		from ...analysis.halo_analysis import halos
		#Check if HOP file exists (note we will adopt a naming convention here)
		if finder == "hop":
			hop_dir = self.hop_path()
	
			#if not os.path.isdir(hop_dir): os.mkdir(hop_dir)
	
			if (verbose): print 'Loading halos from directory: %s'%hop_dir
	
			#Check if out_%05d_hop.h5, .out and .txt exist
			prefix = 'out_%05d_hop'%self._ioutput
			extensions = ['h5', 'txt', 'out']
	
			if run_finder == False:
				for ext in extensions:
					if not os.path.isfile('%s/%s.%s'%(hop_dir, prefix, ext)):
						run_finder = True
						break

			#Return the halos
			if run_finder:
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


		elif finder == "AHF":
			return halos.AHFCatalogue(self)

		elif finder == "rockstar":
			return halos.RockstarCatalogue(self)

		else:
			raise Exception("Unimplemented finder: %s" %finder)

# YT with RAMSES has no way of determining whether particles are stars or dark matter.. for that, we need to make ammends.. for example, stars have a birth time, dark matter particles do not.


		
# standard cosmology

	def _a_dot(self):
		cos = self.cosmology()
		a_dot = (cos["h"]*100) * cos["aexp"] * np.sqrt(cos["omega_m_0"] * (cos["aexp"] ** -3) + cos["omega_k_0"]* (cos["aexp"] ** -2) + cos["omega_l_0"])
		return a_dot
