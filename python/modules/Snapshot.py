'''
Based on Snapshot.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import abc
import config
import numpy as np

#Abstract snapshot class for loading a ramses output

class Snapshot():
	'''
	Abstract class for loading a single ramses snapshot. Pymses/Pynbody specific loaders are based on this class.
	'''
	__metaclass__ = abc.ABCMeta

	def __init__(self, sim_type):
		'''
		Constructor
		'''
		self._type = sim_type

	@abc.abstractmethod
	def output_number(self):
		'''
		Return the output number for this snapshot
		'''
		return

	@abc.abstractmethod
	def path(self):
		'''
		Return the path to this snapshot
		'''
		return

	@abc.abstractmethod
	def raw_snapshot(self):
		'''
		Return the raw snapshot object
		'''
		return

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

	def tform(self, tform, rt=False):

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

		unit_t = self.info()['unit_t']

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

			#time = max( (time_simu - time)/(h0*1e5/3.08e24)/(365*24*3600*1e9), 0 )
			#output[j] = time

			output[j] = (time_simu - time)*unit_t/(365 * 24 * 3600 * 1e9)

		return output

	def integrate_friedman(self, aexp=None, store=False):
		from fortran import friedman as fm

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
