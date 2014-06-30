#Virial overdensity as a function of redshift
import numpy as np

def del_c(z, omega_m=0.3):
	a = 1./(1 + z)
	x = -( (1 - omega_m)*a**3) / (omega_m + (1 - omega_m)*a**3)
	del_c = (178 + 82*x - 39*(x**2))/(1 + x)
	print 'del_c(%f) = %1.2e'%(z, del_c)

	return del_c

#Hoeft et al. 2006 characteristic mass
def hoeft_Mc(z):
	#Tau encodes evolution of min. virial temp.
	tau_z = 0.73 * (z + 1)**0.18 * np.exp(-(0.25 * z)**2.1)
	del_c_z = del_c(z)
	del_c_0 = del_c(0)

	Mc_equ6 = (tau_z * (1./(1+z)))**(3./2.) * (del_c_0/del_c_z)**0.5
	print 'Equ 6 = ', Mc_equ6
	Mc = Mc_equ6 * 10**10
	return Mc