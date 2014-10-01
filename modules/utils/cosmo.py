#Virial overdensity as a function of redshift
import numpy as np

def omega_z(omega0, z):
	return omega0/(omega0 + (1./(1. + z))**3 * (1 - omega0))

def Hubble_z(z, omega_m, H0):
	aexp = 1./(1. + z)
	return ((omega_m/aexp**3) + (1-omega_m))**0.5 * H0

def del_c(z, omega_m=0.3):
	a = 1./(1 + z)
	x = -( (1 - omega_m)*a**3) / (omega_m + (1 - omega_m)*a**3)
	del_c = (178 + 82*x - 39*(x**2))/(1 + x)
	print 'del_c(%f) = %1.2e'%(z, del_c)

	return del_c

#Hoeft et al. 2006 characteristic mass
def hoeft_Mc(z, omega_m=0.3):
	if z < 0: z = 0
	#Tau encodes evolution of min. virial temp.
	tau_z = 0.73 * (z + 1)**0.18 * np.exp(-(0.25 * z)**2.1)
	del_c_z = del_c(z, omega_m=omega_m)
	del_c_0 = del_c(0, omega_m=omega_m)

	Mc_equ6 = (tau_z * (1./(1+z)))**(3./2.) * (del_c_0/del_c_z)**0.5
	print 'Equ 6 = ', Mc_equ6
	Mc = Mc_equ6 * 10**10
	return Mc

def okamoto_Mc(z, Vc, H0, omega_m=0.27):
	import constants as C
	omega_m = omega_z(omega_m, z)
	del_vir = del_c(z, omega_m)
	Vc_3 = Vc**3
	A = 1./np.sqrt(0.5 * del_vir * omega_m)
	B = (1. + z)**(-3./2.)
	return ((Vc_3 * A * B) / (C.G * H0)) / C.Msun # Msun

#Okamoto 2008 Virial Temperature

def T_vir(M, z, omega_m, H0, mu=0.59):
	import constants as C
	delta_vir = del_c(z, omega_m)
	#delta_vir = 1000
	M *= C.Msun
	return 0.5*(mu*C.mp/C.kb) * (delta_vir * omega_m / 2.)**(1./3.) * (1 + z) * (C.G * M * H0)**(2./3.)

def V_c(T_vir = None, mu=0.59, M = None, Rvir = None):
	import constants as C
	if T_vir:
		V_circ = np.sqrt((2 * C.kb * T_vir) / (mu * C.mp))
	elif (M is not None) and (Rvir is not None):
		V_circ = np.sqrt(C.G * M / Rvir)
	else:
		raise Exception("Must supply Virial Temperature or Mass and Virial Radius!")
	return V_circ
