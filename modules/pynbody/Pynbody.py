
'''
Based on Pymses.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import pynbody
from pynbody.units import Unit
from .. import Snapshot
import sys, os
import numpy as np
from ... import config


# multithreading
pynbody.ramses.multiprocess_num = 16
pynbody.config['number_of_threads'] = 16

def load(folder, ioutput=None, **kwargs):
	return PynbodySnapshot(path, ioutput, **kwargs)


class PynbodySnapshot(Snapshot.Snapshot):
	def __init__(self, path, ioutput, gas=True, **kwargs):
		# gas is default to FALSE since pynbody does weird things with gas conversion to tipsy (pymses/yt is better for the gas data)
		Snapshot.Snapshot.__init__(self, path, ioutput, "Pynbody", **kwargs)
		'''
		Load the snapshot using pymses.RamsesOutput
		TODO Verify snapshot exists
		Specify the halo arg to load only the required number of CPUs
		'''
		if kwargs.has_key('halo'):
			halo = kwargs['halo']
			cpus = self.cpu_list(halo.get_bounding_box())
			self._snapshot = pynbody.load('%s/output_%05d'%(path, ioutput), cpus=cpus **kwargs)
		else:
			self._snapshot = pynbody.load('%s/output_%05d'%(path, ioutput), **kwargs)

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
		return self._path
	

	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return self._snapshot

	def ncpu(self):
		'''
		Return the number of CPUs used
		'''
		s = self.raw_snapshot()
		return s.ncpu

	def current_redshift(self):
		'''
		Return the current redshift
		'''
		s = self.raw_snapshot()
		aexp = s.properties['a']
		return (1/aexp) - 1

	def cosmology(self):
		'''
		Return an object with cosmological parameters
		'''
		s = self.raw_snapshot()

		#We need the ramses snapshot for this
		if (type(s) == pynbody.tipsy.TipsySnap):
			ramses_path = os.path.abspath(os.path.join(self.path(), '../../../'))
			snap = load(ramses_path, self._simulation, self.output_number())
			s = snap.raw_snapshot()

		info = s.properties
		aexp = info['a']
		omega_m_0 = info['omegaM0']
		omega_l_0 = info['omegaL0']
		h = info['h']
		aexp = info['a']

# pynbody does not natively load omega_k and omega_b... processing these now from the last snapshot (since these values do not change)
		n = self.output_number()
		info_file = ("%s/output_%05d/info_%05d.txt" % (self._path, n, n))
		f = open(info_file, 'r')
		nline = 1
		while nline <= 18:
			line = f.readline()
			if(nline == 14): omega_k_0 = np.float32(line.split("=")[1])
			if(nline == 15): omega_b_0 = np.float32(line.split("=")[1])
			nline += 1
	
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
		#return self.raw_snapshot().properties
		return self.raw_snapshot()._info

	def analytical_mass_function(self, kern='ST'):
		'''
		Return the Sheth-Tormen mass function
		Units: Mpc^-3 h^3 and Msun/h
		'''
		
		s = self._snapshot
		M, sigma, N = pynbody.analysis.halo_mass_function(s, kern=kern)
		return M, sigma, N


	def tipsy_dir(self):
		path = self.path()
		ioutput = self.output_number()
		output_str = 'output_%05d'%ioutput

		tipsy_dir = '%s/%s/%s_tipsy/'%(path, output_str, output_str)
		return tipsy_dir

	def tipsy_fname(self):
		ioutput = self.output_number()
		output_str = 'output_%05d'%ioutput

		ftipsy = '%s/%s_fullbox.tipsy'%(self.tipsy_dir(), output_str)
		return ftipsy

	def tipsy_units(self):
		import numpy as np
		from ..utils import constants
		s = self.raw_snapshot()

		#We need the ramses snapshot for this
		if (type(s) == pynbody.tipsy.TipsySnap):
			ramses_path = os.path.abspath(os.path.join(self.path(), '../../../'))
			snap = load(ramses_path, self._simulation, self.output_number())
			s = snap.raw_snapshot()

		# figure out the units			 
		cmtokpc = constants.cmtokpc
		G_u = constants.G_u

		lenunit  = s._info['unit_l']/s.properties['a']*cmtokpc
		massunit = pynbody.analysis.cosmology.rho_crit(s, z=0, unit='Msol kpc^-3')*lenunit**3
		timeunit = np.sqrt(1/G_u * lenunit**3/massunit)
		return lenunit, massunit, timeunit

	def tipsy(self, convert=True, gas=True, rewrite=False):
		# note RAMSES-CH does not handle very well with the pynbody conversion... tbh gas converted to particles is bad news anyway
		#Grab the tipsy output for this snapshot, if it exists
		print "if you are using RAMSES-CH, this will probably break, unless you set has_gas = False in RamsesSnap__init__ within pynbody"
		ftipsy = self.tipsy_fname()

		if os.path.exists(ftipsy) and rewrite == False:
			t = load(ftipsy)
			setattr(t, '_ioutput', self.output_number())
			return t
		else:
			if convert:
				"""
				Credit to Rok Roskar (roskar@physik.uzh.ch)
				"""
				if rewrite and os.path.exists(self.tipsy_dir()):
					print "deleting tipsy conversion"
					rmdir = self.tipsy_dir()
					os.system('rm -rf %s'%(rmdir))

				s = self.raw_snapshot()

				newdir = self.tipsy_dir()
				if not os.path.exists(newdir):
					os.mkdir(newdir)
				newfile = self.tipsy_fname()

				lenunit, massunit, timeunit = self.tipsy_units()

				l_unit = Unit('%f kpc'%lenunit)
				m_unit = Unit('%f Msol'%massunit)
				#t_unit = Unit('%f Gyr'%timeunit)

				#l_unit = Unit('Mpc')
				#m_unit = Unit('1e10 Msol')
				t_unit = Unit('%f Gyr'%timeunit)
				v_unit = l_unit/t_unit

				#write_param_file(newfile, s, lenunit, massunit)
				f = open('%s.param'%newfile, 'w')
				f.write('dHubble0 = %f\n'%(s.properties['h']*100))
				f.write('dOmega0 = %f\n'%s.properties['omegaM0'])
				f.write('dLambda = %f\n'%s.properties['omegaL0'])
				f.write('dMsolUnit = %f\n'%massunit)
				f.write('dKpcUnit = %f\n'%lenunit)
				f.close()
			 
				#newfile = "%s_tipsy/%s_fullbox.tipsy"%(s.filename[:-1],outname)
				print 'Writing file %s'%newfile
				#s['mass'].convert_units('%f Msol'%massunit)
				s['mass'].convert_units('%f Msol'%massunit)
				if gas == True and len(s.g) > 0:
					s.g['rho'].convert_units(m_unit/l_unit**3) # get gas variables
					s.g['temp']
					s.g['metals'] = s.g['metal']
					s['eps'] = s.g['smooth'].min()   # smooth the gas
					s['eps'].units = s['pos'].units
					del(s.g['metal'])
					del(s['smooth'])

				s['pos'].convert_units(l_unit)
				s['vel'].convert_units(v_unit)
				
				s.write(filename='%s'%newfile, fmt=pynbody.tipsy.TipsySnap, binary_aux_arrays = True)
				t = load(newfile)
				setattr(t, '_ioutput', self.output_number())
				return t
			else:
				raise Exception("Tipsy file not found: %s"%ftipsy)

### AHf halos using the tipsy outputs

# see here for more doccumentation http://pynbody.github.io/pynbody/tutorials/halos.html


	def tipsy_halos(self, LgridDomain=128, LgridMax=1073741824, VescTune=1.5, Dvir=200, nmin_per_halo = 50,
		 MaxGatherRad=3.0, num_threads=16, run_ahf=False, rewrite_tipsy=False, gas=False, configloc = True):
		import glob
		s = self.raw_snapshot()
		isRamses = (type(s) == pynbody.ramses.RamsesSnap)

		fname = self.tipsy_fname() if isRamses else self.path()

		if run_ahf:
			#Remove the AHF files
			ahf_files = glob.glob('%s.*.AHF_*'%fname)
			for f in ahf_files:
				os.remove(f)

		if self.halo_cat_exists() is False:
			#We need to run AHF
			from ..utils import cosmo
			z = self.current_redshift()
			omega_m_z = cosmo.omega_z(s.properties['omegaM0'], z)

			#First, convert to tipsy
			if os.path.exists(fname) is False and isRamses: self.tipsy(gas=gas, rewrite=rewrite_tipsy)
			
			lenunit, massunit, timeunit = self.tipsy_units()

			l_unit = Unit('%f kpc'%lenunit)
			t_unit = Unit('%f Gyr'%timeunit)
			v_unit = l_unit/t_unit

			f = open('%s.AHF.input'%fname,'w')
			f.write('[AHF]\n')
			f.write('ic_filename = %s\n'%fname)
			f.write('ic_filetype = 90\n')
			f.write('outfile_prefix = %s\n'%fname)
			f.write('LgridDomain = %d\n'%LgridDomain)
			f.write('LgridMax = %d\n'%LgridMax)
			f.write('NperDomCell = 5\n')
			f.write('NperRefCell = 5\n')
			f.write('VescTune = %1.1f\n'%VescTune)
			f.write('NminPerHalo = %d\n'%nmin_per_halo)
			f.write('RhoVir = 0\n')
			f.write('Dvir = %f\n'%Dvir)
			f.write('MaxGatherRad = %1.1f\n'%MaxGatherRad)
			f.write('[TIPSY]\n')
			f.write('TIPSY_BOXSIZE = %e\n'%(s.properties['boxsize'].in_units('Mpc')*s.properties['h']/s.properties['a']))
			f.write('TIPSY_MUNIT   = %e\n'%(massunit*s.properties['h']*(1/omega_m_z)))
			f.write('TIPSY_OMEGA0  = %f\n'%s.properties['omegaM0'])
			f.write('TIPSY_LAMBDA0 = %f\n'%s.properties['omegaL0'])			

		 #   velunit = Unit('%f cm'%s._info['unit_l'])/Unit('%f s'%s._info['unit_t'])			

			f.write('TIPSY_VUNIT   = %e\n'%v_unit.ratio('km s^-1 a', **s.conversion_context()))		
		
			# the thermal energy in K -> km^2/s^2		 
			f.write('TIPSY_EUNIT   = %e\n'%((pynbody.units.k/pynbody.units.m_p).in_units('km^2 s^-2 K^-1')*5./3.))
			f.close()			

			os.environ['OMP_NUM_THREADS'] = str(num_threads)
			if not configloc:
				print "running AHF from your applications dir"
				os.system("~/apps/bin/AHF-v1.0-084 %s.AHF.input"%fname)
			else:
				print "running AHF from your compiled location"
				appstring = config.applications_dir + "/AHF-v1.0-084"
				filename = "%s.AHF.input" %fname
				print filename
				print appstring
				exe_string = appstring + " " + filename
				os.system(exe_string)

		#Return the halos. We need the tipsy snap now to load them
		if isRamses:
			#print 'Warning: Analysis of RAMSES output with tipsy halos can result in unexpected results...'
			s = self.tipsy().raw_snapshot()

		halos = s.halos()
		print 'Loaded %d halos'%len(halos)
		return halos

	def halos(self, finder=config.default_finder, filename=None):
		'''
		Load a generic halo catalogue - default to rockstar if not overridden
		Override the snapshot method, and force loading using yt for unit coherence
		'''
		simulation = self._attributes['simulation']
		yt_snap = self.swap_modules(simulation, 'yt')
		from ...analysis.halo_analysis import halos
		if finder=='rockstar':
			return halos.RockstarCatalogue(yt_snap)
		elif finder=="AHF":
			return halos.AHFCatalogue(yt_snap, filename=filename)
		else:
			raise Exception("Unimplemented finder: %s"%finder)


	def halo_cat_exists(self):
		import glob

		s = self.raw_snapshot()
		isRamses = (type(s) == pynbody.ramses.RamsesSnap)
		fname = self.tipsy_fname() if isRamses else self.path()
		ahf_files = glob.glob('%s.*.AHF_*'%fname)

		return (len(ahf_files) != 0)

	def halos_deprec(self, nmin_per_halo = 150, num_threads=16, configloc = True):
		import glob
		s = self.raw_snapshot()
		snap = self
		isTipsy = (type(s) == pynbody.tipsy.TipsySnap)

		fname = self.path() if isTipsy else self.tipsy_fname()

		if len(glob.glob('%s.*.AHF_*'%fname)) == 0:
			if os.path.exists(fname) is False and isTipsy is False: self.tipsy()
			if isTipsy:
				#Switch to the raw ramses output to run AHF
				ramses_path = os.path.abspath(os.path.join(self.path(), '../../../'))
				snap = load(ramses_path, self._simulation, self.output_number())
				s = snap.raw_snapshot()
			else:
				print 'Warning: Analysis of RAMSES output with tipsy halos can result in unexpected results...'

			print 'No AHF files found, running...'
			#Run AHF. Pynbody has some hooks for this I think, but for now will use existing code (Rok Roskar)
			lenunit, massunit, timeunit = snap.tipsy_units()

			l_unit = Unit('%f kpc'%lenunit)
			t_unit = Unit('%f Gyr'%timeunit)
			v_unit = l_unit/t_unit

			f = open('%s.AHF.input'%fname,'w')
			f.write('[AHF]\n')
			f.write('ic_filename = %s\n'%fname)
			f.write('ic_filetype = 90\n')
			f.write('outfile_prefix = %s\n'%fname)
			f.write('LgridDomain = 256\n')
			f.write('LgridMax = 2097152\n')
			f.write('NperDomCell = 5\n')
			f.write('NperRefCell = 5\n')
			f.write('VescTune = 1.0\n')
			f.write('NminPerHalo = %d\n'%nmin_per_halo)
			f.write('RhoVir = 0\n')
			f.write('Dvir = 200\n')
			f.write('MaxGatherRad = 1.0\n')
			f.write('[TIPSY]\n')
			f.write('TIPSY_BOXSIZE = %e\n'%(s.properties['boxsize'].in_units('Mpc')*s.properties['h']/s.properties['a']))
			f.write('TIPSY_MUNIT   = %e\n'%(massunit*s.properties['h']))#/s.properties['omegaM0']))
			f.write('TIPSY_OMEGA0  = %f\n'%s.properties['omegaM0'])
			f.write('TIPSY_LAMBDA0 = %f\n'%s.properties['omegaL0'])
			
		 #   velunit = Unit('%f cm'%s._info['unit_l'])/Unit('%f s'%s._info['unit_t'])
			
			f.write('TIPSY_VUNIT   = %e\n'%v_unit.ratio('km s^-1 a', **s.conversion_context()))			
		 
			# the thermal energy in K -> km^2/s^2
		 
			f.write('TIPSY_EUNIT   = %e\n'%((pynbody.units.k/pynbody.units.m_p).in_units('km^2 s^-2 K^-1')*5./3.))
			f.close()			

			os.environ['OMP_NUM_THREADS'] = str(num_threads)

			if not configloc:
				print "running AHF from your applications dir"
				os.system("~/apps/bin/AHF-v1.0-075 %s.AHF.input"%fname)
			else:
				print "running AHF from your compiled location"
				appstring = config.applications_dir + "/AHF-v1.0-084"
				filename = "%s.AHF.input" %fname
				print filename
				print appstring
				exe_string = appstring + " " + filename
				os.system(exe_string)	

		if isTipsy is False:
			s = self.tipsy().raw_snapshot()
		return s.halos()

class Type:
	GAS = 1
	STAR = 2
	DM = 3
