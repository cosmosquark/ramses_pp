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
	f.write(write_string('omega_b',sim_info['omega_l_0']))
	f.write(write_string('lambda_f',sim_info['omega_b_0']))
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
