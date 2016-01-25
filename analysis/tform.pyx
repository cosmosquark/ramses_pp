""" 
Taken from getstarlist.f90 in ramses/utils/f90 by R. Teysier
Python code written by D. Sullivan
converted into Cython and adapted for YT by B. Thompson
"""

import numpy as np
cimport numpy as np


def tform(np.ndarray tform, np.ndarray tau_out, np.ndarray t_out, double time_simu, double hubble_constant):

    """
    hubble constant in units of cm/s
    """

    cdef np.ndarray[dtype="double", ndim=1] output = np.zeros((len(tform)))
    cdef int i
    cdef int j
    cdef double time
    cdef int ntable
    
    ntable = len(tau_out)-1

    for j in range(0, len(tform)-1):

        if tform[j] != 0.0:

            i = 1
            while ( (tau_out[i] > tform[j]) and (i < ntable) ):
                i+=1

            #Interpolate time
            time = t_out[i] * (tform[j] - tau_out[i-1]) / (tau_out[i] - tau_out[i-1]) + \
                 t_out[i-1] * (tform[j] - tau_out[i]) / (tau_out[i-1] - tau_out[i])

            time = np.max( (time_simu - time)/(hubble_constant), 0 )
            output[j] = time
       
        else:
            output[j] = 0.0

    return output

def friedman(double omega_matter, double omega_lambda, double omega_curvature,
    double alpha, double axp_min, int ntable):

    """
    This function assumes that axp = 1 at z = 0 (today) 
    and that t and tau = 0 at z = 0 (today).              
    axp is the expansion factor, hexp the Hubble constant 
    defined as hexp=1/axp*daxp/dtau, tau the conformal    
    time, and t the look-back time, both in unit of 1/H0. 
    alpha is the required accuracy and axp_min is the     
    starting expansion factor of the look-up table.       
    ntable is the required size of the look-up table.     
    """

    print "integrating friedman for RAMSES"
    # defining output variables

    cdef double age_tot
    cdef np.ndarray[dtype="double", ndim=1] axp_out = np.zeros((ntable))
    cdef np.ndarray[dtype="double", ndim=1] hexp_out = np.zeros((ntable))
    cdef np.ndarray[dtype="double", ndim=1] tau_out = np.zeros((ntable))
    cdef np.ndarray[dtype="double", ndim=1] t_out = np.zeros((ntable))

    # defining local variables

    cdef double axp_tau
    cdef double axp_t
    cdef double axp_tau_pre
    cdef double axp_t_pre
    cdef double dtau 
    cdef double dt
    cdef double tau
    cdef double t
    cdef int nstep
    cdef int nout
    cdef int nskip

 #   if (omega_matter + omega_lambda + omega_curvature != 1.0):
 #       print('Error: non-physical cosmological constants')
 #       print("omega_matter, omega_lambda, omega_curvature = ", 
 #                omega_matter, omega_lambda, omega_curvature)
 #       print("The sum mus equal 1.0, but")
 #       print("omega_matter + omega_lambda + omega_curvature = ", 
 #                omega_matter + omega_lambda + omega_curvature)

    axp_tau = 1.0
    axp_t = 1.0
    tau = 0.0
    t = 0.0
    nstep = 0

    print('Cosmology: ', omega_matter, omega_lambda, omega_curvature)
    print('Params: ', alpha,axp_min)

    while (( axp_tau >= axp_min ) or ( axp_t >= axp_min )):
        nstep += 1
        dtau = alpha * axp_tau / dadtau(axp_tau,omega_matter,
                                     omega_lambda,omega_curvature)
        axp_tau_pre = axp_tau - (dadtau(axp_tau,omega_matter,
                                     omega_lambda,omega_curvature) * dtau / 2.0 )
        axp_tau = axp_tau - (dadtau(axp_tau_pre,omega_matter,
                                     omega_lambda,omega_curvature) * dtau)
        tau = tau - dtau
     
        dt = alpha * axp_t / dadt(axp_t,omega_matter,
                                     omega_lambda,omega_curvature)
        axp_t_pre = axp_t - ( dadt(axp_t,omega_matter,
                                     omega_lambda,omega_curvature) * dt / 2.0 )
        axp_t = axp_t - ( dadt(axp_t_pre,omega_matter,
                                     omega_lambda,omega_curvature) * dt)
        t = t - dt

    age_tot = -t

    nskip = nstep / ntable
  
    axp_t = 1.0
    t = 0.0
    axp_tau = 1.0
    tau = 0.0
    nstep = 0
    nout=0
    t_out[nout] = t
    tau_out[nout] =tau
    axp_out[nout] = axp_tau
    hexp_out[nout] = dadtau(axp_tau,omega_matter,omega_lambda,omega_curvature) / axp_tau

    while ( (axp_tau >= axp_min) or (axp_t >= axp_min) ):
        nstep += 1
        dtau = alpha * axp_tau / dadtau(axp_tau,omega_matter,
                               omega_lambda,omega_curvature)
        axp_tau_pre = axp_tau - (dadtau(axp_tau,omega_matter,
                               omega_lambda,omega_curvature) * dtau / 2.0)
        axp_tau = axp_tau - (dadtau(axp_tau_pre,omega_matter,
                               omega_lambda,omega_curvature)*dtau)
        tau = tau - dtau

        dt = alpha * axp_t / dadt(axp_t,omega_matter,
                               omega_lambda,omega_curvature)
        axp_t_pre = axp_t - (dadt(axp_t,omega_matter,
                               omega_lambda,omega_curvature) * dt / 2.0)
        axp_t = axp_t - (dadt(axp_t_pre,omega_matter,
                               omega_lambda,omega_curvature) * dt)
        t = t - dt
     
        if np.mod(nstep,nskip) == 0:
            t_out[nout] = t
            tau_out[nout] = tau
            axp_out[nout] = axp_tau
            hexp_out[nout] = dadtau(axp_tau,omega_matter,
                               omega_lambda,omega_curvature)/axp_tau
            nout += 1    

    t_out[ntable-1] = t
    tau_out[ntable-1] = tau
    axp_out[ntable-1] = axp_tau
    hexp_out[ntable-1] = dadtau(axp_tau,omega_matter,
                               omega_lambda,omega_curvature)/axp_tau

    return axp_out, hexp_out, tau_out, t_out, age_tot




def dadtau(double axp_tau, double omega_matter,double omega_lambda,
              double omega_curvature):
    cdef double dadtau 

    dadtau = axp_tau*axp_tau*axp_tau *  \
           ( omega_matter + \
             omega_lambda * axp_tau * axp_tau * axp_tau + \
             omega_curvature * axp_tau )

    dadtau = np.sqrt(dadtau)
    return dadtau


def dadt(double axp_t,double omega_matter,double omega_lambda,
            double omega_curvature):

    cdef double dadt

    dadt   = (1.0 / axp_t) * \
          ( omega_matter + \
            (omega_lambda * axp_t * axp_t * axp_t) + \
            omega_curvature * axp_t )

    dadt = np.sqrt(dadt)
    return dadt

