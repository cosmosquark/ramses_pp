from ramses_pp import config
import os


def write_string(var, val):
	widthform = 20
	string = var.ljust(widthform) + " = " + str(val) + '\n'
	return string

def write_string_final(var, val):
	widthform = 20
	string = var.ljust(widthform) + " = " + str(val)
	return string


def make_inputfiles(direct, nsteps, snapno, snaploc, simtype="Ra3"):
	filename = direct + '/inputfiles_HaloMaker.dat'
	if os.path.exists(filename):
		os.remove(filename)	

	f = open(filename, 'w')
	string = "'" + snaploc + " '  " + simtype + "   " + nsteps + "  " + snapno
	f.write(string)
	f.close()

def make_input(direct, sim_info, halomaker_stats, nsteps):
	filename = direct + '/input_HaloMaker.dat'
	# check if file already exists
	if os.path.exists(filename):
		os.remove(filename)	

	# write the file
	f = open(filename, 'w')
	f.write(write_string('af',sim_info['faexp']))
	f.write(write_string('lbox',sim_info['lbox(Mpc)']))
	f.write(write_string('H_f',sim_info['H0']))
	f.write(write_string('omega_f',sim_info['omega_m_0']))
	f.write(write_string('omega_b',sim_info['omega_b_0']))
	f.write(write_string('lambda_f',sim_info['omega_l_0']))
	f.write(write_string('FlagPeriod',sim_info['FlagPeriod']))
	f.write(write_string('npart',halomaker_stats['npart']))
	f.write(write_string('method',halomaker_stats['method']))
	f.write(write_string('cdm',halomaker_stats['cdm']))
	f.write(write_string('b',halomaker_stats['b']))
	f.write(write_string('nsteps',nsteps))
	f.write(write_string('nvoisins',halomaker_stats['adaptahop']['nvoisins']))
	f.write(write_string('nhop',halomaker_stats['adaptahop']['nhop']))
	f.write(write_string('rhot',halomaker_stats['adaptahop']['rhot']))
	f.write(write_string('fudge',halomaker_stats['adaptahop']['fudge']))
	f.write(write_string('fudgepsilon',halomaker_stats['adaptahop']['fudgepsilon']))
	f.write(write_string('alphap',halomaker_stats['adaptahop']['alphap']))
	f.write(write_string('verbose',halomaker_stats['verbose']))
	f.write(write_string_final('megaverbose',halomaker_stats['adaptahop']['megaverbose']))

	f.close()
	return

class HaloMaker():
	def __init__(self, Simulation, subvol=False, ncpu=1):
		'''
		Store initial properties and simulation infomation
		'''
		self.Simulation = Simulation
		self._subvol=subvol
		self._ncpu=ncpu
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
			 } 

		is_simulation = hasattr(Simulation, 'initial_conditions')
		if is_simulation == False:
			print "this is not a simulation object"
			return None
	
		data_dir = self.Simulation.path() + "/HaloMaker"

		# make the halomaker data dir if it does not exist
		if not os.path.exists(data_dir):
			os.makedirs(data_dir)

		print "loading halomaker location"
		application_dir = config.applications_dir
		halomaker_dir = application_dir + "/HaloMaker/f90"
		post_analysis_dir = application_dir + "/HaloMaker/post_analysis"
		TreeMaker_dir = application_dir + "/HaloMaker/TreeMaker"
		GalaxyMaker_dir = application_dir + "/HaloMaker/GalaxyMaker"
		self._halomaker_exe = halomaker_dir + "/HaloFinder"
		self._halomaker_dir = halomaker_dir
		self._post_analysis_dir = post_analysis_dir
		self._TreeMaker_dir = TreeMaker_dir
		self._GalaxyMaker_dir = GalaxyMaker_dir
		self._data_dir = data_dir
		print "loading initial variables"



		self._nsteps = "1"


	def halomaker_info(self):
		widthform = 14
		print "current halomaker variables as defined are"
		print "##################################"
		print "method".ljust(widthform) +  str(self._halomaker_info['method'])
		print "b".ljust(widthform) + str(self._halomaker_info['b'])
		print "cdm".ljust(widthform) + str(self._halomaker_info['cdm'])
		print "npart".ljust(widthform) + str(self._halomaker_info['npart'])
		print "nsteps".ljust(widthform) +  str(self._nsteps).ljust(widthform), "(not editable)"
		print "nvoisins".ljust(widthform) + str(self._halomaker_info['adaptahop']['nvoisins'])
		print "nhop".ljust(widthform) + str(self._halomaker_info['adaptahop']['nhop'])
		print "rhot".ljust(widthform) + str(self._halomaker_info['adaptahop']['rhot'])
		print "fudge".ljust(widthform) + str(self._halomaker_info['adaptahop']['fudge'])
		print "fudgepsilon".ljust(widthform) + str(self._halomaker_info['adaptahop']['fudgepsilon'])
		print "alphap".ljust(widthform) + str(self._halomaker_info['adaptahop']['alphap'])
		print "megaverbose".ljust(widthform) + str(self._halomaker_info['adaptahop']['megaverbose'])
		print "verbose".ljust(widthform) + str(self._halomaker_info['verbose'])
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

## part 1, halomaker

	def check_halofinder_snap(self,ioutput):
		direct = self._data_dir + '/' + str(ioutput)
		# does the directory exist? if not run it
		if not os.path.exists(direct):
			os.makedirs(direct)
		# check if run is not already done
		# first by checking if there is any out files
		props = "%03d" % ioutput
		fileprops = direct + "/" + "haloProps." + props
		if os.path.exists(fileprops):
			print "halomaker already run on " + str(ioutput)
			return True
		# then by checking the log file
		logfile = direct + "/log.out"
		if os.path.exists(logfile):
			run = False
			lines = [line.strip() for line in open(logfile)]
			for i in range(0,len(lines)):
				if lines[i] == "End of HaloMaker":
					run = True		
			if run == True:
				print "halomaker already run on " + str(ioutput)
				return True
		return False
		


	def run_halomaker_snap(self,ioutput):
		'''
		runs halomaker and treemaker preperations on a snapshot
		'''
		snapno = "%05d" % ioutput
		snaploc = self.Simulation._path + "/" + "output_" + snapno + '/'
		nsteps = self._nsteps
		direct = self._data_dir + '/' + str(ioutput)

		check = self.check_halofinder_snap(ioutput)
		if check == False: # if false, run halomaker
		# make the directory
						
			print "working with " + direct

		# make input files
			make_input(direct, self.sim_info(), self._halomaker_info, self._nsteps)
			make_inputfiles(direct, nsteps, snapno, snaploc, simtype="Ra3")
			os.chdir(direct)
			command = "ln -s " + self._halomaker_exe + " ."
			os.system(command)
			command = "chmod u+x HaloFinder"
			os.system(command)
			command = "./HaloFinder > log.out"
			os.system(command)
			command = "rm HaloFinder"
			os.system(command)
		else:
			"halomaker already run, skipping"

		print "preparing " + str(ioutput) + " for treemaker"
		
		os.chdir(direct)  # this directory must exist to get past halomaker check
		treeprepexe = self._post_analysis_dir + "/treebrick2ascii"
		naploc = self.Simulation._path + "/" + "output_" + snapno
		outnod = "halomaker_" + str(snapno) + ".nod"
		outnei = "halomaker_" + str(snapno) + ".nei"
		outbricks = "tree_bricks" + str("%03d" % ioutput)
		treecommand = treeprepexe + " -inp " + naploc + " " + outnod + " " + outnei + " " + outbricks
		print treecommand
		os.system(treecommand)
		print "treemaker done, assuming there are halos"


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
				self.run_halomaker_snap(ioutput)
		
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
		


