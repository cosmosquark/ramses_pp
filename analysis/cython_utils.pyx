import numpy as np
cimport numpy as np


# loop over all the stars for supernova yields


def compute_mbins(np.ndarray rad_fields, np.ndarray r_bins, np.ndarray mass_fields, np.ndarray m_bins):

	cdef int i
	cdef int j

	j = 0
	for i in range(0,len(mass_fields)):
		if rad_fields[i] > r_bins[j]:
			while (rad_fields[i] > r_bins[j]) and j < (len(mass_fields) - 1):
				j = j + 1
				m_bins[j] = m_bins[j-1]
	#			print j
		m_bins[j] = m_bins[j] + mass_fields[i]

	return m_bins
	

def cic_sample(np.ndarray ihx, np.ndarray y, np.ndarray x, np.ndarray x_bins, int nbins):

	cdef np.ndarray[dtype="double", ndim=1] pdf = np.zeros((nbins+1))
	cdef int i

	for i in range(0,len(y)):
		pdf[ihx[i]] = pdf[ihx[i]] + y[i] * np.power((x_bins[ihx[i]+1]-x[i]) / (x_bins[ihx[i]+1] - x_bins[ihx[i]]),2.0)
		pdf[ihx[i]+1] = pdf[ihx[i]+1] + y[i] * np.power((x[i] - x_bins[ihx[i]]) / (x_bins[ihx[i]+1] - x_bins[ihx[i]]),2.0)
	return pdf

def supernova_yield_loop(double[:] formation_time, double[:] metallicity, double[:] mass_formation, 
			int nstar, int nbins, double dy, double dx_1, double z_min, double time_min,
			np.ndarray yieldtab_zstar, np.ndarray yieldtab_astar, np.ndarray dt,
			np.ndarray yieldtab_NSN1, 
			np.ndarray yieldtab_NSN2,):

	# define variables here
	

	cdef double yy, lbt, lbt0, xx
	cdef int i, j, ihx, ihy 
	cdef double wya, wyb, wxa, wxb, yieldval
	cdef np.ndarray[dtype="double", ndim=1] nsn1 = np.zeros((nbins))
	cdef np.ndarray[dtype="double", ndim=1] nsn2 = np.zeros((nbins))

	print nsn1
	

	for j in range(0, nstar):
		lbt = formation_time[j]
	
		# the Z of this star
		yy = metallicity[j]

		# star Z
		if yy > 0.0000000000000000000:
			ihy = np.int(np.floor(dy * (np.log10(yy) - z_min))) + 1

			if (yy < yieldtab_zstar[1]):  # fix to minimum
				ihy = 0

		else: # if 0, then np.log10(yy) is infinity).. cython needs to be more specific
			ihy = 0

		# find the weightings
		wya = (yy-yieldtab_zstar[ihy])  /(yieldtab_zstar[ihy+1]-yieldtab_zstar[ihy])
		wyb = (yieldtab_zstar[(ihy+1)]-yy)/(yieldtab_zstar[ihy+1]-yieldtab_zstar[ihy])

		for i in range(0,nbins): # looping through metallicity
			lbt0 = dt[i] # lbt0 is the time in the cosmic table
			xx = lbt0-lbt # table time - cosmic time of formation = the amount of the table this star contributes at the given time (NB not actually sampling the late end of the table)
			if (xx >= 0.0 ):

				ihx = np.array(np.floor(dx_1  * (np.log10(xx) - time_min)), dtype="int")
				if (ihx >= 0):
					wxa = (xx-yieldtab_astar[ihx]) / (yieldtab_astar[ihx+1] - yieldtab_astar[ihx])
					wxb = (yieldtab_astar[ihx+1] - xx)/(yieldtab_astar[ihx+1] - yieldtab_astar[ihx])


					if (yieldtab_NSN1[ihx,ihy] >= 0.0):
						yieldval = (yieldtab_NSN1[ihx,ihy] * wxb * wyb) \
								+ (yieldtab_NSN1[ihx,ihy+1] * wxb * wya) \
								+ (yieldtab_NSN1[ihx+1,ihy] * wxa * wyb) \
								+ (yieldtab_NSN1[ihx+1,ihy+1] * wxa * wya)

						if (yieldval >= 0.): 
							nsn1[i] = nsn1[i] + yieldval * mass_formation[j]
         				
					if (yieldtab_NSN2[ihx,ihy] >= 0.0):
						yieldval = (yieldtab_NSN2[ihx,ihy] * wxb * wyb) \
								+ (yieldtab_NSN2[ihx,ihy+1] * wxb * wya) \
								+ (yieldtab_NSN2[ihx+1,ihy] * wxa * wyb) \
								+ (yieldtab_NSN2[ihx+1,ihy+1] * wxa * wya)

						if (yieldval >= 0.):
							nsn2[i] = nsn2[i] + yieldval * mass_formation[j]

	return nsn1, nsn2
