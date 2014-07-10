import loader
#from ramses_pp.modules import Simulation, Snapshot
from ramses_pp import config
import halomaker_utils
import os

class Halomaker():
	def __init__(self, Simulation, subvol=False, ncpu=1):
		self.Simulation = Simulation
		self._subvol=subvol
		self._ncpu=ncpu
		is_simulation = hasattr(Simulation, 'initial_conditions')
		if is_simulation == False:
			print "this is not a simulation object"
			return None

		is_simulation = hasattr(Simulation, '_halomaker_info')
		if is_simulation == False:
			print "No Halomaker data is stored. Returning"
			return None
	
		halomaker_dir = self.Simulation.data_dir() + "/HaloMaker"

		# make the halomaker data dir if it does not exist
		if not os.path.exists(halomaker_dir):
			os.makedirs(halomaker_dir)

		print "loading halomaker location"
		application_dir = config.applications_dir
		data_dir = config.json_dir + '/' + str(self.Simulation.name()) + '/' + 'HaloMaker'
		halomaker_dir = application_dir + "/HaloMaker/f90"
		self._halomaker_exe = halomaker_dir + "/HaloFinder"
		self._halomaker_dir = halomaker_dir
		self._data_dir = data_dir
		print "loading initial variables"



		self._nsteps = "1"
		self._halomaker_stats = self.Simulation.halomaker_info()[0]
		self.halomaker_info()


	def halomaker_info(self):
		widthform = 14
		print "current halomaker variables as defined are"
		print "##################################"
		print "method".ljust(widthform), self._halomaker_stats['method']
		print "b".ljust(widthform), self._halomaker_stats['b']
		print "cdm".ljust(widthform), self._halomaker_stats['cdm']
		print "npart".ljust(widthform), self._halomaker_stats['npart']
		print "nsteps".ljust(widthform), str(self._nsteps).ljust(widthform), "(not editable)"
		print "nvoisins".ljust(widthform), self._halomaker_stats['adaptahop']['nvoisins']
		print "nhop".ljust(widthform), self._halomaker_stats['adaptahop']['nhop']
		print "rhot".ljust(widthform), self._halomaker_stats['adaptahop']['rhot']
		print "fudge".ljust(widthform), self._halomaker_stats['adaptahop']['fudge']
		print "fudgepsilon".ljust(widthform), self._halomaker_stats['adaptahop']['fudgepsilon']
		print "alphap".ljust(widthform), self._halomaker_stats['adaptahop']['alphap']
		print "megaverbose".ljust(widthform), self._halomaker_stats['adaptahop']['megaverbose']
		print "verbose".ljust(widthform), self._halomaker_stats['verbose']
		print "##################################"
		print ""
		print "use .edit_halomaker_info(variable=new) to update any of this info"
		print "also: if ever aexp final is greater than 1 in sim_info, it is then set to 1"
		return

	def sim_info(self):
		'''
		returns relevent infomation for this simulation in the context of HaloMaker. See the readme
		'''

		initial_conditions = self.Simulation.initial_conditions()
		af = initial_conditions['faexp'] #expansion factor at redshift 0
		if af > 1.0:
			af = 1.0
		lbox = self.Simulation.box_size() / initial_conditions['h']
		H0 = initial_conditions['H0']
		omega_m_0 = initial_conditions['omega_m_0']
		omega_l_0 = initial_conditions['omega_l_0']
		omega_b_0 = initial_conditions['omega_b_0']
		periodic = self.Simulation.is_periodic()
		if periodic:
			FlagPeriod = 1
		else:
			FlagPeriod = 0
		sim_info = {
			'faexp':af,				# final expansion parameter
			'lbox(Mpc)':lbox,		# !NB in MPC, not Mpc/h  size of the box at redshift 0, [lbox] =  MPc
			'H0':H0, 					
			'omega_m_0':omega_m_0,  # Omega matter at redshift 0
			'omega_l_0':omega_l_0,  # Omega lambda at redshift 0
			'omega_b_0':omega_b_0,
			'FlagPeriod':FlagPeriod,  # 1 for periodic boundary conditions, 0 else
		}
		return sim_info

	def run_halomaker(self,parallel=False,ncpu=1):
		'''
		executes halomaker across all snapshots
		'''	

		# gather number of snapshots
		num = self.Simulation.num_snapshots()
		# mpi version eventually
		if parallel==False:
			print "running HaloMaker in serial"
			for ioutput in range(1,num+1):
				# make the directory
				direct = self._data_dir + '/' + str(ioutput)
				print "working with " + direct
				if not os.path.exists(direct):
					os.makedirs(direct)
				# make input files
				halomaker_utils.make_input(direct, self.sim_info(), self._halomaker_stats, self._nsteps)
				snapno = "%05d" % ioutput
				snaploc = self.Simulation._path + "/" + "output_" + snapno + '/'
				nsteps = self._nsteps
				halomaker_utils.make_inputfiles(direct, nsteps, snapno, snaploc, simtype="Ra3")
				
		return
							
		
		

	def executable_loc(self):
		'''
		returns the executable location of halomaker.
		'''
		return self._halomaker_exe

	def halomaker_dir(self):
		'''
		returns the directory halomaker is stored in.
		'''
		return self._halomaker_dir

	def output_dir(self):
		'''
		returns the director halomaker outputs are stored in
		'''
		return self._data_dir


	

## edit function coming soon when I can be arsed to write it
#	def edit_halomaker_info(self,method=None,b=None,cdm=None,adaptahop=None,nvoisins=None,nhop=None,rhot=None,fudge=None,fudgeepsilon=None,alphap=None,megaverbose=None,verbose=None):
		







		

