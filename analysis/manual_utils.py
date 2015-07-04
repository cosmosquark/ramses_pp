# these are analysis utilities
# that are like utils.py
# but do NOT use YT data sources (i.e can take in any numpy array as a dataset)
from ramses_pp.analysis import read_utils, utils
import numpy as np
import scipy as sp
import yt
from ramses_pp.analysis.plot_utils import flatten_line, kill_first_nan
from ramses_pp.analysis.read_utils import read_yield_table
from matplotlib import pyplot as plt
from ramses_pp.analysis.cython_utils import supernova_yield_loop, cic_sample
from yt.utilities.math_utils import ortho_find
from ramses_pp.analysis.visualisation_utils import gen_data_source

def sfr(formation_time, mass_formation, n_bins=80, time_range = [-1.0,14.0], flip_axis = False):

	"""
	Inputs:
		Mass : Solar Mass
		Formation time: Gyr

	Returns:
		Time: binned time array
		SFR: Star formation rate (Msun/ Year)
	"""

	print "TIMES"
	print formation_time.min(), formation_time.max()

	time_range = [formation_time.min(),formation_time.max()]
	hist, bins = np.histogram(formation_time, bins=n_bins, range=time_range,)
	inds = np.digitize(formation_time, bins=bins)
	inds = inds - 1  # shifts everything down one.

	time = (bins[:-1] + bins[1:]) / 2.0

	sfr = np.array([mass_formation[inds == j].sum()/((bins[j+1]-bins[j])*(10**9)) for j in range(0,len(time))])  # conversion from Gyr to Yr
	sfr[sfr == 0] = np.nan


	if flip_axis == True:
		age = np.fliplr([time])[0] # just need to flip the axis
	else:
		age = time


	return age, sfr 


def sfrhist(formation_time,mass_formation,no_bins=100):

	t_min = np.log10(formation_time.min())
	t_max = np.log10(formation_time.max())
	delt = float(no_bins-1)/(t_max-t_min)
	t = np.zeros((no_bins+1)) # want the extra bin for ihx
	pdf = np.zeros((no_bins))

	for i in range(0,no_bins+1):
		t[i]=np.power(10, (t_min + ( float(i-1)*(t_max-t_min)/float(no_bins-1) ) ) )
	
	# log10 spacing
	ihx = np.array(np.floor(delt  * (np.log10(formation_time) - t_min)), dtype="int")

	# cic sampling
	pdf = cic_sample(ihx, mass_formation, formation_time, t, no_bins)

	# convert into Msun / yr
	pdf = pdf / (10. ** 9.0)		
	pdf = pdf[np.where(np.log10(t) <= t_max)]
	t = t[np.where(np.log10(t) <= t_max)]

	# get the pdf on the right axis
#	if flip_axis == True:
		

	return t, pdf


def supernova(yield_file,formation_time,mass_formation,metallicity, no_bins = 100):
	# load the yield table	
	yield_table = read_yield_table(yield_file)

	# COMPLICATED GUBBINS FOR CALCULATING THE SNR as Dr. Few would say
	nmetal = yield_table.shape[0]    #n22
	nsteps = len(yield_table[0,1]["time"])  #n11
	yieldtab_astar = yield_table[0,1]["time"]
	yieldtab_zstar = yield_table[:,0]

	yieldtab_NSN1 = np.zeros((nsteps,nmetal))
	yieldtab_NSN2 = np.zeros((nsteps,nmetal))

	# populate arrays for fast acess

	for i in range(0,nmetal):
		yieldtab_NSN1[:,i] = yield_table[:,1][i]["N_SNIa"]
		yieldtab_NSN2[:,i] = yield_table[:,1][i]["N_SNII"]

	time_min = np.log10(yieldtab_astar.min()) 
	time_max = np.log10(yieldtab_astar.max())

	# size of x bins
	dx_1 = np.float(nsteps-1)/(time_max - time_min) # width of bin

	# convert into cummulatitive format

	for j in range(0,len(yieldtab_astar)):
		yieldtab_astar[j] = np.power(10, (time_min + (float(j) * (time_max-time_min))/float(nsteps-1)) )

	z_min = np.log10(yieldtab_zstar[1])
	z_max = np.log10(yieldtab_zstar.max())

	# size of y bins
	dy = np.float(nmetal-2)/(z_max-z_min)
	
	# dt is actually the time bins at the CENTER of each dx_1

	# nbins = n11 or nsteps -1 
	nbins = nsteps - 1
	
	# initialise arrays

	dt = np.zeros((nbins))
	nsn1 = np.zeros((nbins))
	nsn2 = np.zeros((nbins))

	for j in range(0,nbins):
		dt[j] = (yieldtab_astar[j+1] + yieldtab_astar[j])/2.0

	# This block is converting the table to a non cumulative format

	for i in range(0,nmetal):
		tmp = np.zeros((nbins))
		for j in range(0,nbins):
			delta=yieldtab_astar[j+1] - yieldtab_astar[j]
			tmp[j]=(yieldtab_NSN1[j+1,i]-yieldtab_NSN1[j,i])/delta

		yieldtab_NSN1[:nbins,i]=tmp[:]

		tmp = np.zeros((nbins))
		for j in range(0,nbins):
			delta=yieldtab_astar[j+1] - yieldtab_astar[j]
			tmp[j]=(yieldtab_NSN2[j+1,i]-yieldtab_NSN2[j,i])/delta

		yieldtab_NSN2[:nbins,i]=tmp[:]

	# find SNR for stellar

	nstar = len(formation_time)

	# debug

	nsn1, nsn2 = supernova_yield_loop(abs(formation_time), metallicity, mass_formation,
			nstar, nbins, dy, dx_1, z_min, time_min,
			yieldtab_zstar, yieldtab_astar, dt,
			yieldtab_NSN1, yieldtab_NSN2)

	# convert to solar masses per year

	nsn1 = nsn1 * 100 / (10 ** 9)
	nsn2 = nsn2 * 100 / (10 ** 9)

	# we now have the number of SN as a function of time for a number of
	# steps equal to that in the yield table. the next step is to bin these
	# values into a more managable set of linear bins (since the times are log in the table)
    	

	jjmax = no_bins # i.e the output bins
	dx = np.zeros((jjmax))
	wgt = dx

	# define new age bins
	t_min = np.log10(0.1)
	t_max = np.log10(formation_time.max())

	dx_1 = np.float((jjmax - 1) / (t_max - t_min))

	for i in range(0,jjmax):
		dx[i] = np.power(10, (t_min + (float(i) * (t_max-t_min))/float(jjmax-1.0)) )

	massive = np.array(np.floor(dx_1 * (np.log10(max(dt)) - t_min ) ), dtype="int") - 1  #??

	tab = np.zeros((massive,2)) # tab needs initialising properly, it's the rate
	wgt = np.zeros((massive))


	for i in range(0,nbins): # cyclle over the table

		dd1 = dt[i]

		if dd1 < dx[0]:
			dd1 = dx[0] # fix the earliest up so we account for them

		if dd1 <= dx[jjmax-1]:
			
			ihx = np.array(np.floor( dx_1 * (np.log10(dd1) - t_min) ), dtype="int" )
			dd = nsn1[i]

			tab[ihx,0] = tab[ihx,0] + dd * (dx[ihx+1] - dd1) / (dx[ihx+1] - dx[ihx])  # cloud in cell sampling
			tab[ihx+1,0] = tab[ihx+1,0] + dd * (dd1 - dx[ihx]) / (dx[ihx+1] - dx[ihx])

			wgt[ihx] = wgt[ihx] + ( dx[ihx+1] - dd1) / (dx[ihx+1] - dx[ihx])
			wgt[ihx+1] = wgt[ihx+1] + ( dd1 - dx[ihx] ) / (dx[ihx+1] - dx[ihx]) # weight


	for i in range(0,jjmax):
		if (wgt[i] > 0.0):
			tab[i,0] = (tab[i,0] / wgt[i]) * 5.0 # averaging part

	# same again for SNII

	wgt = np.zeros((massive))

	for i in range(0,nbins):
					
		dd1 = dt[i]

		if dd1 < dx[0]:
			dd1 = dx[0] # fix the earliest up so we account for them

		if dd1 <= dx[jjmax-1]:
			
			ihx = np.array(np.floor( dx_1 * (np.log10(dd1) - t_min) ), dtype="int") 
			dd = nsn2[i]
			
			tab[ihx,1] = tab[ihx,1] + dd * (dx[ihx+1] - dd1) / (dx[ihx+1] - dx[ihx])  # cloud in cell sampling
			tab[ihx+1,1] = tab[ihx+1,1] + dd * (dd1 - dx[ihx]) / (dx[ihx+1] - dx[ihx])

			wgt[ihx] = wgt[ihx] + ( dx[ihx+1] - dd1) / (dx[ihx+1] - dx[ihx])
			wgt[ihx+1] = wgt[ihx+1] + ( dd1 - dx[ihx] ) / (dx[ihx+1] - dx[ihx]) # weight

	for i in range(0,jjmax):
		if (wgt[i] >= 0.0):
			tab[i,1] = tab[i,1] / wgt[i] # averaging part


	return dx, tab[:jjmax,0], tab[:jjmax,1]
	


def mstar_mhalo(ytsnap):
	"""
	inputs:
		stars in Msun
		dark in Msun

	outputs:
		log10 of mass of DM M_dm
		log10 of M*/Mdm

	Returns the result of the mass ratio of the DM to the stellar matter and the corresponding x value
	"""
	dm_mass = ytsnap["dark","particle_mass"].sum().in_units("Msun")
	stellar_mass = ytsnap["stars","particle_mass"].sum().in_units("Msun")

	mass_ratio = stellar_mass / dm_mass
	mass_log = np.log10(mass_ratio)

	print mass_ratio, "mass_ratio"
	print mass_log, "mass_log"
	print np.log10(dm_mass.value), "dm mass log"

	x = np.log10(dm_mass.value)
	return x, mass_log


def ks_law(ytdataset,raw_snapshot,width,max_age,image_width=800,field="stars",n_bins=100, r_max = None):

	"""
	Use this to compute the ks_law for each individual pixel within a galaxy
	"""




	weight_field = None

	if r_max == None:
		r_max = ytdataset.ds.arr(25,"kpc")
	
	print "Computing KS law"
	print "ks 1) Gathering initial projections"

	width = ((width.v,width.units),(width.v,width.units))
	image_width = (image_width,image_width) 

	# slice of the star particle mass
	normal = ytdataset.get_field_parameter("normal")
	data_source = gen_data_source(0,ytdataset,raw_snapshot,width,width[0],"kpc")

	
	plot = yt.OffAxisProjectionPlot(data_source.ds,normal,[("gas","density"),("deposit","%s_density" % field)],center=ytdataset.center,width=width,depth=width[0],weight_field=weight_field, north_vector=normal)
	# set the units
	plot.set_axes_unit("kpc")
	plot.set_unit(("gas","density"),"Msun/pc**2")
	plot.set_unit(("deposit","%s_density" % field),"Msun/kpc**2")
	images = plot.frb
	gas_image = images[('gas', 'density')]
	star_image = images[('deposit', '%s_density' % field)]

	print "ks 2) collecting data"

	# filtering things withing the first 15 kpc does the trick usuall


	n_bins = gas_image.shape[0]
	bin_width = width[0][0] / float(n_bins)
	left = - width[0][0] / 2.0
	right = width[0][0] / 2.0
	length = np.arange(left,right,bin_width)
	shift_thing = (length[2] - length[1]) / 2.0
	length = length + shift_thing

	# to make radial array, we need to calculate the radial distance of each pixel in a 2d array. Namely.

	length_squared = np.power(length,2.0) # power of 2
	radius_squared = np.add.outer(length_squared, length_squared) # adds as an outer product
	radius = np.sqrt(radius_squared)

	# now we need to generate some radial bins

	radius_bin_width = radius.max() / float(n_bins)
	radius_bins = np.arange(0,radius.max(), radius_bin_width)

	# flatten the arrays

	flat_radius = radius.flatten()

	# now need to find the indices in which radius is maximum
	radius_filter = (flat_radius < r_max.in_units("kpc").v)

	# ok, now get the data
	gas_density = np.array(gas_image.in_units("Msun/pc**2").value).flatten()
	star_density = np.array(star_image.in_units("Msun/kpc**2").value).flatten()

	# filter by radius
	gas_density = np.log10(gas_density[radius_filter])
	star_density = np.log10(star_density[radius_filter] / max_age )

	return gas_density, star_density




def vcirc(ytdataset, fields=["all","stars","dark","gas"],r_min=None,r_max=None,dr_factor=2):
	"""
	Inputs:
		Ytdataset (either particles or not).. with "stars", "dark", and "gas" defined (or not)
		fields = list containing ["all","stars","dark","gas"

	returns

		r_bins = list of numpy arrays corresponding to the field order
		v_circ = same as above except with the vcirc results	
			
	"""


	if r_max == None:
		r_max = ytdataset["particle_position_cylindrical_radius"].max().in_units("kpc")
#		r_max = sphere.ds.arr(50,"kpc")

	if r_min == None:
		r_min = ytdataset.ds.arr(0,"kpc")

	# 1) get radial values
	r_indices, r_bins, r_filter, r_truths = utils.radial_indices('particle_position_spherical_radius',ytdataset,r_min=None,r_max=r_max,dr_factor=dr_factor)


	n_bins = len(r_bins)
	r_bins = ytdataset.ds.arr(r_bins,"cm")

	vcirc_list = []
	r_bins_list = []

	for field in fields:

		if field == "all":
			m_bins = utils.mass_enclosed_bins(ytdataset,ytdataset,r_bins,shape="sphere",type=field,sim_object_type="ad", AMR=True)
			vcirc = utils.vcirc(r_bins,m_bins,ytdataset)

			vcirc_list.append(vcirc.in_units("km/s").value)
			r_bins_list.append(r_bins.in_units("kpc").value)

		if field == "stars":
			m_bins = utils.mass_enclosed_bins(ytdataset,ytdataset,r_bins,shape="sphere",type=field,sim_object_type="ad", AMR=True)
			vcirc = utils.vcirc(r_bins,m_bins,ytdataset)

			vcirc_list.append(vcirc.in_units("km/s").value)
			r_bins_list.append(r_bins.in_units("kpc").value)

		if field == "dark":
			m_bins = utils.mass_enclosed_bins(ytdataset,ytdataset,r_bins,shape="sphere",type=field,sim_object_type="ad", AMR = True)
			vcirc = utils.vcirc(r_bins,m_bins,ytdataset)

			vcirc_list.append(vcirc.in_units("km/s").value)
			r_bins_list.append(r_bins.in_units("kpc").value)


		if field == "gas":
			m_bins = utils.mass_enclosed_bins(ytdataset,ytdataset,r_bins,shape="sphere",type=field,sim_object_type="ad", AMR=True)
			vcirc = utils.vcirc(r_bins,m_bins,ytdataset)

			vcirc_list.append(vcirc.in_units("km/s").value)
			r_bins_list.append(r_bins.in_units("kpc").value)


	return r_bins_list, vcirc_list


def vrot(ytdataset, fields=["stars","young_stars","disk_stars","gas","cold_gas"],n_bins=50, AMR=True):


	"""
	Input

		ytdataset: A ytdataset with gas cells or gas particles.
			particulally of the region you want to work with (e.g a disk)
		fields: list of fields to compute vrot for

	output
		r_bins_list = list of r)bins in kpc
		
	"""

	vrot_list = []
	rbins_list = []

	r_max = ytdataset["particle_position_cylindrical_radius"].max().in_units("kpc")
	r_indices, r_bins, r_filter, r_truths = utils.radial_indices('particle_position_spherical_radius',ytdataset,r_min=None,r_max=r_max,dr_factor=2)
	for field in fields:

		if field == "gas" and AMR == True:
			vrot = yt.ProfilePlot(ytdataset,'cylindrical_r',["velocity_cylindrical_theta"],weight_field="cell_mass",n_bins=n_bins,x_log=False,y_log={"velocity_cylindrical_theta":False})
			r_binned = vrot.profiles[0].x.in_units("kpc").value
			vrot_binned = np.abs(vrot.profiles[0]["velocity_cylindrical_theta"]).in_units("km/s").value # this overwrites the yt plot object
		

		elif field == "cold_gas" and AMR == True:

			cold_disk = ytdataset.cut_region(["obj['temperature'] < 1e4"])
			vrot = yt.ProfilePlot(cold_disk,'cylindrical_r',["velocity_cylindrical_theta"],weight_field="cell_mass",n_bins=n_bins,x_log=False,y_log={"velocity_cylindrical_theta":False})
			r_binned = vrot.profiles[0].x.in_units("kpc").value
			vrot_binned = np.abs(vrot.profiles[0]["velocity_cylindrical_theta"]).in_units("km/s").value	

		else:

			if field != "disk_stars":
				vrot = yt.ProfilePlot(ytdataset,(field,'particle_position_cylindrical_radius'),[(field,'particle_velocity_cylindrical_theta')],weight_field=(field,'particle_mass'),n_bins=n_bins,x_log=False,y_log={(field,'particle_velocity_cylindrical_theta'):False})	
				r_binned = vrot.profiles[0].x.in_units("kpc").value
				vrot_binned = np.abs(vrot.profiles[0][(field,'particle_velocity_cylindrical_theta')]).in_units("km/s").value
			
			else:
				# manual profile plot
				hist, bins = np.histogram(ytdataset[field,'particle_position_cylindrical_radius'].in_units("kpc").value, bins=n_bins)
				inds = np.digitize(ytdataset[field,'particle_position_cylindrical_radius'].in_units("kpc").value, bins=bins)
				width = (bins[:-1] + bins[1:]) / 2.0
      				  #time = time[:i1]
				# want to compute the mass weight ... vrot * mass_i / mass_mean in each bin
				#  average of sum of (vrot * mass_i)  / mass_tot

				

				vrot_thing = np.array([ np.sum(ytdataset[field,"particle_velocity_cylindrical_theta"][inds == j].in_units("km/s").value * ytdataset[field,"particle_mass"][inds == j].in_units("Msun").value) for j in range(len(bins))])
				vrot_binned = vrot_thing / np.array([ (ytdataset[field,"particle_mass"][inds == j].in_units("Msun").value).sum()  for j in range(len(bins))])

				#vrot_binned[vrot_binned == 0] = np.nan
				vrot_binned[vrot_binned == np.nan] = 0.0
				vrot_binned = abs(vrot_binned)
			
				r_binned = bins
				

			# make things pretty
	#	r_binned, vrot_binned = flatten_line(r_binned,vrot_binned,no_nan=True,append_max=True,double_zero=True,extra_x = r_bins.in_units("kpc").max())
		
		vrot_list.append(vrot_binned)
		rbins_list.append(r_binned)
								
	
	return rbins_list, vrot_list


def velocity_dispersion(ytdataset, fields=["stars","disk_stars","solar_stars"], max_age = 14.5, bins = 30, AMR=True, time_field="particle_birth_epoch"):


	"""
	Requirements
		ytdataset
		field list

	returns

		for particles
			age_bins
			vel_dispersion
			vel_disp_err

		for gas

			sig_u
			sig_v
			sig_w
	"""

	age_fields = []
	vel_disp_fields = []
	vel_err_fields = []


	for field in fields:

		if field == "gas" and AMR == True:
			u = ytdataset["velocity_cylindrical_radius"].in_units("km/s").value
			v = ytdataset["velocity_cylindrical_theta"].in_units("km/s").value
			w = ytdataset["velocity_cylindrical_z"].in_units("km/s").value
			age_fields.append(u)
			vel_disp_fields.append(v)
			vel_err_fields.append(w)

		elif field == "cold_gas" and AMR == True:
			cold_cylinder = ytdataset.cut_region(["obj['temperature'] < 1.0e4"])
			u = cold_cylinder["velocity_cylindrical_radius"].in_units("km/s").value
			v = cold_cylinder["velocity_cylindrical_theta"].in_units("km/s").value
			w = cold_cylinder["velocity_cylindrical_z"].in_units("km/s").value

			age_fields.append(u)
			vel_disp_fields.append(v)
			vel_err_fields.append(w)

		else:
			u = ytdataset[field,"particle_velocity_cylindrical_radius"].in_units("km/s").value
			v = ytdataset[field,"particle_velocity_cylindrical_theta"].in_units("km/s").value
			w = ytdataset[field,"particle_velocity_cylindrical_z"].in_units("km/s").value
			age = abs(ytdataset[field,time_field].in_units("Gyr").value)

			sig_u = np.std(u)
			sig_v = np.std(v)
			sig_w = np.std(w)

			bin_things = np.linspace(0.0,max_age,bins)
			hist, bin_edges = np.histogram(age,bins=bin_things) # gets the distribution
			inds = np.digitize(age, bins=bin_edges) # put the ages into the right bins
			inds = inds - 1.0
	
			bin_edges = np.delete(bin_edges,len(bin_edges)-1) # deletes the last bin edge elemtn

			sig_u_age = np.zeros((bins-1))
			sig_v_age = np.zeros((bins-1))
			sig_w_age = np.zeros((bins-1))

			for i in range(0,len(sig_u_age)):
				ages_inds = np.where(inds == i) # find the index locations of stars of this gage
				sig_u_age[i] = np.std(u[ages_inds[0]])
				sig_v_age[i] = np.std(v[ages_inds[0]])
				sig_w_age[i] = np.std(w[ages_inds[0]])


			#	temp_age = age
			sig_age_temp = sig_u_age
			age_bins, sig_u_age = kill_first_nan(bin_edges,sig_u_age)
			age_bins, sig_v_age = kill_first_nan(bin_edges,sig_v_age)
			age_bins, sig_w_age = kill_first_nan(bin_edges,sig_w_age)
			
			#	freq_hist, dud = plot_utils.kill_first_nan(hist, sig_age_temp)
			freq_hist = hist

			# errors
			sig_u_age_err = sig_u_age / np.sqrt(freq_hist)
			sig_v_age_err = sig_v_age / np.sqrt(freq_hist)
			sig_w_age_err = sig_w_age / np.sqrt(freq_hist)

			vel_dispersion = np.sqrt(np.power(sig_u_age,2) + np.power(sig_v_age,2) + np.power(sig_w_age,2))

			# power errors
			# http://www.rit.edu/cos/uphysics/uncertainties/Uncertaintiespart2.html
			u_err_power = (sig_u_age_err * sig_u_age * 2.0) 
			v_err_power = (sig_v_age_err * sig_v_age * 2.0)
			w_err_power = (sig_w_age_err * sig_w_age * 2.0)

			uvw_err_sq = np.sqrt(np.power(sig_u_age_err,2) + np.power(sig_v_age_err,2) + np.power(sig_w_age,2))
			vel_disp_err = (uvw_err_sq * (0.5)) / vel_dispersion

			# add the errors
			vel_disp_err = vel_dispersion / np.sqrt(freq_hist)


			age_fields.append(age_bins)
			vel_disp_fields.append(vel_dispersion)
			vel_err_fields.append(vel_disp_err)
		
	return age_fields, vel_disp_fields, vel_err_fields



def velocity_dispersion_r(ytdataset, fields=["stars","disk_stars","solar_stars"], max_radius = 30, bins = 30, AMR=True):


	"""
	Requirements
		ytdataset
		field list

	returns

		for particles
			rad_bins
			vel_dispersion
			vel_disp_err

		for gas

			sig_u
			sig_v
			sig_w
	"""

	rad_fields = []
	vel_disp_fields = []
	vel_err_fields = []


	for field in fields:

		if field == "gas" and AMR == True:
			u = ytdataset["velocity_cylindrical_radius"].in_units("km/s").value
			v = ytdataset["velocity_cylindrical_theta"].in_units("km/s").value
			w = ytdataset["velocity_cylindrical_z"].in_units("km/s").value
			rad_fields.append(u)
			vel_disp_fields.append(v)
			vel_err_fields.append(w)

		elif field == "cold_gas" and AMR == True:
			cold_cylinder = ytdataset.cut_region(["obj['temperature'] < 1.0e4"])
			u = cold_cylinder["velocity_cylindrical_radius"].in_units("km/s").value
			v = cold_cylinder["velocity_cylindrical_theta"].in_units("km/s").value
			w = cold_cylinder["velocity_cylindrical_z"].in_units("km/s").value

			rad_fields.append(u)
			vel_disp_fields.append(v)
			vel_err_fields.append(w)

		else:
			u = ytdataset[field,"particle_velocity_cylindrical_radius"].in_units("km/s").value
			v = ytdataset[field,"particle_velocity_cylindrical_theta"].in_units("km/s").value
			w = ytdataset[field,"particle_velocity_cylindrical_z"].in_units("km/s").value
			rad = ytdataset[field,"particle_position_cylindrical_radius"].in_units("kpc").value

			sig_u = np.std(u)
			sig_v = np.std(v)
			sig_w = np.std(w)

			bin_things = np.linspace(0.0,max_rad,bins)
			hist, bin_edges = np.histogram(rad,bins=bin_things) # gets the distribution
			inds = np.digitize(rad, bins=bin_edges) # put the rads into the right bins
			inds = inds - 1.0
	
			bin_edges = np.delete(bin_edges,len(bin_edges)-1) # deletes the last bin edge elemtn

			sig_u_rad = np.zeros((bins-1))
			sig_v_rad = np.zeros((bins-1))
			sig_w_rad = np.zeros((bins-1))

			for i in range(0,len(sig_u_rad)):
				rad_inds = np.where(inds == i) # find the index locations of stars of this rad
				sig_u_rad[i] = np.std(u[rad_inds[0]])
				sig_v_rad[i] = np.std(v[rad_inds[0]])
				sig_w_rad[i] = np.std(w[rad_inds[0]])


			#	temp_rad = rad
			sig_rad_temp = sig_u_rad
			rad_bins, sig_u_rad = kill_first_nan(bin_edges,sig_u_rad)
			rad_bins, sig_v_rad = kill_first_nan(bin_edges,sig_v_rad)
			rad_bins, sig_w_rad = kill_first_nan(bin_edges,sig_w_rad)
			
			#	freq_hist, dud = plot_utils.kill_first_nan(hist, sig_age_temp)
			freq_hist = hist

			# errors
			sig_u_rad_err = sig_u_rad / np.sqrt(freq_hist)
			sig_v_rad_err = sig_v_rad / np.sqrt(freq_hist)
			sig_w_rad_err = sig_w_rad / np.sqrt(freq_hist)

			vel_dispersion = np.sqrt(np.power(sig_u_rad,2) + np.power(sig_v_rad,2) + np.power(sig_w_rad,2))

			# power errors
			# http://www.rit.edu/cos/uphysics/uncertainties/Uncertaintiespart2.html
			u_err_power = (sig_u_rad_err * sig_u_rad * 2.0) 
			v_err_power = (sig_v_rad_err * sig_v_rad * 2.0)
			w_err_power = (sig_w_rad_err * sig_w_rad * 2.0)

			uvw_err_sq = np.sqrt(np.power(sig_u_rad_err,2) + np.power(sig_v_rad_err,2) + np.power(sig_w_rad,2))
			vel_disp_err = (uvw_err_sq * (0.5)) / vel_dispersion

			# add the errors
			vel_disp_err = vel_dispersion / np.sqrt(freq_hist)


			rad_fields.append(rad_bins)
			vel_disp_fields.append(vel_dispersion)
			vel_err_fields.append(vel_disp_err)
		
	return rad_fields, vel_disp_fields, vel_err_fields



