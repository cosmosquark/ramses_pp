import pymses, sys, os
from pynbody.array import SimArray
import numpy as np

def critical_density(omega_m0, h):
	G = 6.67384E-11
	#Hubble paramater in SI units
	H = (h*100) * 1000/3.08E22

	return (3*H**2)/(8*np.pi*G)

def main(omega_m0, h, boxsize, ngrid):

	rho_crit0 = critical_density(omega_m0, h)
	mean_rho0 = SimArray(omega_m0*rho_crit0, "kg m**-3")
	mean_rho0.convert_units('Msol Mpc**-3')

	print 'Mean density: %e'%mean_rho0, ' Msol Mpc**-3'

	mean_M0 = mean_rho0 * (boxsize / h)**3 * h

	print 'Mean mass: %e'%mean_M0, ' Msol h**-1'

	mass_resolution = mean_M0/ngrid

	print 'Mass resolution: %e'%mass_resolution, ' Msol h**-1'
	print 'Halo mass with 50 particles: %e'%(mass_resolution*50), ' Msol h**-1'
	print 'Halo mass with 150 particles: %e'%(mass_resolution*150), ' Msol h**-1'

def jeans(n_star, T, mu=0.76):
	kB = 1.38e-23
	mH = 1.67e-27
	G = 6.67384E-11

	#Assume n_star in units H/cc
	n_star_SI = n_star * 1e6 #H/m^3
	rho_star_SI = n_star_SI * mH # kg/m^3

	r_J = ((15 * kB * T)/(4 * np.pi * G * mu * mH * rho_star_SI))**0.5
	print 'For T = %e, r_{J} = %e m, %e pc'%(T, r_J, r_J/3.08e16)

	return r_J

def n_jeans(delta_x):
	#See Teyssier, R et al. 2010. delta_x in pc, returns n_J in H/cc
	n_J = 6*((delta_x/100)**(-4./3.))
	return n_J

def min_temp(n):
	#See Teyssier, R et al. 2010. n is number density in units H/cc
	T_eq = 10**4 * (n/0.3)**(-0.5)
	return T_eq

if __name__ == "__main__":

	if (len(sys.argv) == 1):
		print 'Case 1 Usage: python %s <simdir> <ioutput> <nstar>'%sys.argv[0]
		print 'Case 2 Usage: python %s <omega_m0> <h> <L_box (Mpc/h)> <ell_min> <ell_max> <nstar>'%sys.argv[0]
		sys.exit(1)

	elif os.path.isdir(sys.argv[1]):

		ro = pymses.RamsesOutput(sys.argv[1], int(sys.argv[2]))
		omega_m0 = ro.info['omega_m']
		omega_b = ro.info['omega_b']
		h = ro.info['H0'] / 100
		aexp = ro.info['aexp']
		boxsize = ro.info['unit_length'].val / 3.08e22 / aexp * h # Mpc/h at z=0
		lgridmin = ro.info['levelmin']
		lgridmax = ro.info['levelmax']
		n_star = float(sys.argv[3]) # H/cc
	else:
		omega_m0 = float(sys.argv[1])
		omega_b = 0.045
		h = float(sys.argv[2])
		boxsize = float(sys.argv[3]) # Mpc/h
		lgridmin = int(sys.argv[4])
		lgridmax = int(sys.argv[5])
		n_star = float(sys.argv[6]) # H/cc

	ngrid_min = (2**lgridmin)**3 # Assume 3D!
	ngrid_max = (2**lgridmax)**3
	cell_size_coarse = SimArray((boxsize / 2**lgridmin), "Mpc h**-1")
	cell_size_fine = SimArray((boxsize / 2**lgridmax), "Mpc h**-1")
	unit_coarse = "kpc h**-1"
	unit_fine = "pc h**-1"
	print cell_size_coarse, cell_size_fine
	print 'boxsize = %e'%boxsize, 'Mpc /h'
	print 'levelmin = %s, levelmax = %s'%(lgridmin, lgridmax)
	print 'ngrid_min = %s, ngrid_max = %s'%(ngrid_min**(1./3.), ngrid_max**(1./3.))
	print 'omega_m0 = %s, omega_b = %s, h=%s kms**-1 Mpc**-1, boxsize = %e Mpc h**-1'%(omega_m0, omega_b, h, boxsize)
	print 'coarse resolution = ', cell_size_coarse.in_units(unit_coarse), '%s, fine resolution = '%unit_coarse, cell_size_fine.in_units(unit_fine), unit_fine

	print '---------------------------------------------'
	n_J = n_jeans(cell_size_fine.in_units(unit_fine))
	print 'Jeans density = %f H/cc. n_star = %f H/cc'%(n_J, n_star)
	print 'Suggested n_star (based on Tomassetti et al. 2014 and Jeans refinement with 4 cells) = %f'%((2./3.)*n_J)
	print 'NOTE: n_star should be significantly smaller than this value to properly resolve the star forming part of the gas PDF'
	if (n_star > n_J):
		print "ERROR: The chosen value of n_star is greater than the Jeans density! Exiting."
		sys.exit(1)
	T_eq = min_temp(n_J)
	print 'Minimum equillibrium temperature at the Jeans density (neglecting self-shielding), T_eq = %f K'%T_eq

	print '---------------------------------------------'
	print 'STAR FORMATION:'
	T = 100
	r_J = jeans(n_star, T)
	print 'Jeans length resolved with ~ %f cells prior to reionization'%((2 * (r_J/3.08e16)) / (cell_size_fine.in_units(unit_fine)/h))

	print '---------------------------------------------'
	T = 10**4
	r_J = jeans(n_star, T)
	print 'Jeans length resolved with ~ %f cells after reionization'%((2 * (r_J/3.08e16)) / (cell_size_fine.in_units(unit_fine)/h))

	#Compute resolution for gas and dm separately
	print '---------------------------------------------'
	print 'MASS RESOLUTION:'
	print 'Dark Matter coarse grid mass resolution (with hydro):'
	main(omega_m0-omega_b, h, boxsize, ngrid_min)
	print '---------------------------------------------'
	print 'Total Matter coarse grid mass resolution (nbody only):'
	main(omega_m0, h, boxsize, ngrid_min)
	print '---------------------------------------------'
	print 'Done'
