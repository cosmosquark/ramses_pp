'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''
#import pymses
from pymses.filters import RegionFilter
from pymses.utils.regions import Sphere
from pymses.utils import constants as C
from pymses.filters import CellsToPoints
from .. import Profile
import numpy as np
from ramses_pp import config

def load(snapshot):
	return PymsesProfile(snapshot)

class PymsesProfile(Profile.Profile):
	def __init__(self, snapshot):
		Profile.Profile.__init__(self, snapshot, 'pymses')
		'''
		TODO Verift snapshot is a Pymses snapshot
		'''
		self._snapshot = snapshot

	def mstar_halomass_tree(self, tree, ds):

		snapshot = self._snapshot
		ro = snapshot.raw_snapshot()
		parts = ro.particle_source(["mass", "epoch"])
		amr = ro.amr_source(['rho'])

		#Figure out the smallest possible cell size in code units
		boxlen = ro.info['boxlen']
		lmax = ro.info['levelmax']

		min_dx = boxlen/float(2**(lmax))

		offset = tree.get_offset()
		snap_num = self._snapshot.output_number() - offset

		halos = tree.sub_catalogue('snap_num', snap_num)

		#idx = np.where(tree._tree[:]['snap_num'] == self._snapshot.output_number()-30)[0]
		#halos = tree[idx] # Is this right?

		print 'Loaded halos have snap_num:', halos[1]['snap_num']

		mstar = []
		mhalo = []

		i = 1
		for halo in halos:
			#Check if distinct halo
			if not halo['pid'].value == -1:
				print 'Halo %d skipped: PID = %d'%(halo['id'].value, halo['pid'].value)
				continue

			#Unit conversion
			cen = halo['pos']
			cen = ds.arr(cen.value, cen.unit)
			r = halo['Rvir']
			r = ds.arr(r.value, r.unit)

			#Create the sphere
			cen = cen.in_units('code_length').value
			r = r.in_units('code_length').value
			r_amr = float(r)

			#If the halo radius is smaller than the min cell size,
			#we need to try and approximate the masses enclosed.
			#Not sure how correct this is...
			factor = 1
			if r < min_dx:
				r_amr += min_dx
				vol_halo = (4./3.)*np.pi*r**3
				vol_cell = (4./3.)*np.pi*r_amr**3
				factor = vol_halo/vol_cell

			#Define the regions
			region_part = Sphere(cen, r)
			region_amr = Sphere(cen, r_amr)

			#Filter the particles
			filt_parts = RegionFilter(region_part, parts)			
			part_source = filt_parts.flatten()

			dm = np.where(part_source['epoch'] == 0)[0]
			stars = np.where(part_source['epoch'] != 0)[0]

			part_mass = part_source['mass']*ro.info['unit_mass'].express(C.Msun)

			if len(part_mass[dm]) < 50:
				print 'Discarding halo with less than 200 particles. Mvir=%e'%halo['Mvir'].value
				continue

			#Filter the AMR data
			filt_amr = RegionFilter(region_amr, amr)
			cell_source = CellsToPoints(filt_amr)
			cells = cell_source.flatten()
			rho = cells['rho']*ro.info['unit_density'].express(C.Msun/C.kpc**3)
			vol = (cells.get_sizes()*ro.info['unit_length'].express(C.kpc))**3
			cell_mass = rho*vol

			#Compute the total mass enclosed
			gas_mass = np.sum(cell_mass)*factor # Multiply by factor incase we had to inflate the sphere
			stellar_mass = np.sum(part_mass[stars])
			particle_mass = np.sum(part_mass)
			total_mass = gas_mass + particle_mass

			mstar.append(stellar_mass)
			mhalo.append(total_mass)

			i+=1
			print i
			#print i
			if(config.verbose and (i%100)==0): print 'Processed %d halos...'%i
		return np.array(mstar), np.array(mhalo)

	def fgas_halomass_tree(self, tree, ds):

		snapshot = self._snapshot
		ro = snapshot.raw_snapshot()
		amr = ro.amr_source(['rho'])
		parts = ro.particle_source(["mass", "epoch"])

		#Figure out the smallest possible cell size in code units
		boxlen = ro.info['boxlen']
		lmax = ro.info['levelmax']

		min_dx = boxlen/float(2**(lmax))

		offset = tree.get_offset()
		snap_num = self._snapshot.output_number() - offset

		halos = tree.sub_catalogue('snap_num', snap_num)

		#idx = np.where(tree._tree[:]['snap_num'] == self._snapshot.output_number()-30)[0]
		#halos = tree[idx] # Is this right?

		print 'Loaded %d halos'%(len(halos))
		print 'Loaded halos have snap_num:', halos[0]['snap_num']

		fgas = []
		mhalo = []

		i = 1
		for halo in halos:
			#Check if distinct halo
			if not halo['pid'].value == -1:
				#print 'Halo %d skipped: PID = %d'%(halo['id'].value, halo['pid'].value)
				continue

			#Unit conversion
			cen = halo['pos']
			cen = ds.arr(cen.value, cen.unit)
			r = halo['Rvir']
			r = ds.arr(r.value, r.unit)

			#Create the sphere
			cen = cen.in_units('code_length').value
			r = r.in_units('code_length').value
			r_amr = float(r)

			#If the halo radius is smaller than the min cell size,
			#we need to try and approximate the masses enclosed.
			#Not sure how correct this is...
			factor = 1
			if r < min_dx:
				r_amr += min_dx
				vol_halo = (4./3.)*np.pi*r**3
				vol_cell = (4./3.)*np.pi*r_amr**3
				factor = vol_halo/vol_cell

			#Define the regions
			region_part = Sphere(cen, r)
			region_amr = Sphere(cen, r_amr)

			#Filter the particles
			filt_parts = RegionFilter(region_part, parts)			
			part_source = filt_parts.flatten()

			dm = np.where(part_source['epoch'] == 0)

			part_mass = part_source['mass']*ro.info['unit_mass'].express(C.Msun)

			if len(part_mass[dm]) < 200:
				#print 'Discarding halo with less than 200 particles. Mvir=%e'%halo['Mvir'].value
				continue

			#Filter the AMR data
			filt_amr = RegionFilter(region_amr, amr)
			cell_source = CellsToPoints(filt_amr)
			cells = cell_source.flatten()
			rho = cells['rho']*ro.info['unit_density'].express(C.Msun/C.kpc**3)
			vol = (cells.get_sizes()*ro.info['unit_length'].express(C.kpc))**3
			cell_mass = rho*vol

			#Compute the total mass enclosed
			gas_mass = np.sum(cell_mass)*factor # Multiply by factor encase we had to inflate the sphere
			particle_mass = np.sum(part_mass)
			total_mass = gas_mass + particle_mass

			fgas.append(gas_mass/total_mass)
			mhalo.append(total_mass)

			i+=1
			if(config.verbose and (i%100)==0): print 'Processed %d halos...'%i
		return np.array(fgas), np.array(mhalo)

	def fgas_halomass(self, halos=None):

		snapshot = self._snapshot
		ro = snapshot.raw_snapshot()
		amr = ro.amr_source(['rho'])
		parts = ro.particle_source(["mass"])

		#Figure out the smallest possible cell size in code units
		boxlen = ro.info['boxlen']
		lmax = ro.info['levelmax']

		min_dx = boxlen/float(2**(lmax))

		if config.verbose: print 'min_dx=', min_dx

		if (halos == None):
			halos = snapshot.halos()

		if config.verbose: print 'Processing %d halos'%len(halos)

		fgas = []
		mhalo = []

		#z = self._snapshot.current_redshift()
		#a = 1./(1.+z)

		i = 0
		for halo in halos:
			if halo['num_p'] < 150:
				break

			#if tree:
				#Find this halo and count number of progenitors.
				#If too high, skip (frequent mergers -> tidal stripping)

				#idx_tree_halo = (tree._tree[:]['orig_halo_id'] == halo['id'].value) &\
				#		 (tree._tree[:]['snap_num'] == snapshot.rockstar_output_number())
				#print idx_tree_halo
				#tree_halo = tree._tree[idx_tree_halo]
				#assert(len(tree_halo)==1)
				#all_halos = []
				#tree.find_progs(np.array([tree_halo]), all_halos)
				#num_progs = len(all_halos)
				#print 'Halo %d has %d progenitors'%(halo['id'], num_progs)
				#if num_progs > 150: continue

			#Create the sphere
			cen = halo['pos'].in_units('code_length').value
			r = halo['Rvir'].in_units('code_length').value
			r_amr = float(r)

			#If the halo radius is smaller than the min cell size,
			#we need to try and approximate the masses enclosed.
			#Not sure how correct this is...
			factor = 1
			if r < min_dx:
				r_amr += min_dx
				vol_halo = (4./3.)*np.pi*r**3
				vol_cell = (4./3.)*np.pi*r_amr**3
				factor = vol_halo/vol_cell

			#Define the regions
			region_part = Sphere(cen, r)
			region_amr = Sphere(cen, r_amr)

			#Filter the particles
			filt_parts = RegionFilter(region_part, parts)
			part_source = filt_parts.flatten()
			part_mass = part_source['mass']*ro.info['unit_mass'].express(C.Msun)

			#Filter the AMR data
			filt_amr = RegionFilter(region_amr, amr)
			cell_source = CellsToPoints(filt_amr)
			cells = cell_source.flatten()
			rho = cells['rho']*ro.info['unit_density'].express(C.Msun/C.kpc**3)
			vol = (cells.get_sizes()*ro.info['unit_length'].express(C.kpc))**3
			cell_mass = rho*vol

			#Compute the total mass enclosed
			gas_mass = np.sum(cell_mass)*factor # Multiply by factor incase we had to inflate the sphere
			particle_mass = np.sum(part_mass)
			total_mass = gas_mass + particle_mass

			fgas.append(gas_mass/total_mass)
			mhalo.append(total_mass)

			i+=1
			#print i
			if(config.verbose and (i%100)==0): print 'Processes %d halos...'%i
		return np.array(fgas), np.array(mhalo)

