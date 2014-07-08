import numpy as np
from ramses_pp.fortran import test
from ramses_pp.modules import Simulation


sim = Simulation.load('4Mpc_aton_256')
snapshot = sim.snapshot(215)

cosmo = snapshot.cosmology()

omega_m_0 = cosmo['omega_m_0']
omega_l_0 = cosmo['omega_l_0']
omega_k_0 = cosmo['omega_k_0']
aexp = cosmo['aexp']
h0 = cosmo['h']*100

alpha=1e-6
axpmin=1e-3
ntable=1000
age_tot=0


axp_out, hexp_out, tau_out, t_out, age_tot = test.friedman(omega_m_0, omega_l_0, omega_k_0, alpha, axpmin, ntable)
print age_tot

#Find neighbouring expansion factor
i = 1
while ( (axp_out[i] > aexp) and (i < ntable) ):
	i+=1


#Interpolate time
time_simu = t_out[i] * (aexp - axp_out[i-1])/(axp_out[i]-axp_out[i-1]) + \
	t_out[i-1]*(aexp - axp_out[i])/(axp_out[i-1] - axp_out[i])

print 'time_simu = ', time_simu

age_simu = (time_simu+age_tot)/(h0*1e5/3.08e24)/(365*24*3600*1e9)
print 'Age of Universe= ', age_simu
