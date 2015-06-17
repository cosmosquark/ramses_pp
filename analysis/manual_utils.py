# these are analysis utilities
# that are like utils.py
# but do NOT use YT data sources (i.e can take in any numpy array as a dataset)
from ramses_pp.analysis import read_utils, utils
import numpy as np
import scipy as sp
import yt
from ramses_pp.analysis.plot_utils import flatten_line, kill_first_nan

def sfr(mass, formation_time, n_bins=50, time_range = [0.0,14.0]):

	"""
	Inputs:
		Mass : Solar Mass
		Formation time: Gyr

	Returns:
		Time: binned time array
		SFR: Star formation rate (Msun/ Year)
	"""

	time_range = [0.0,14.0] # Giga Years

        #print formation_time, "this"
	hist, bins = np.histogram(formation_time, bins=n_bins, range=time_range,)
	inds = np.digitize(formation_time, bins=bins)

        #ime = max(time) - time
	time = (bins[:-1] + bins[1:]) / 2.0
        #time = time[:i1]


	sfr = np.array([mass[inds == j].sum()/((bins[j+1]-bins[j])*(10**9)) for j in range(len(time))])  # conversion from Gyr to Yr
	sfr[sfr == 0] = np.nan
	
	age = max(time) - time  # is this not a convoluted time = time[:-1]  quick and dirty.. really we want to be calling the year from the snapshot time

	return age, sfr 


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
			print ytdataset[field,'particle_velocity_cylindrical_theta']
			print len(ytdataset[field,'particle_velocity_cylindrical_theta']), field

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


def velocity_dispersion(ytdataset, fields=["stars","disk_stars","solar_stars"], max_age = 14.5, bins = 30 ):


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

	if "Temperature" in ytdataset.ds.derived_field_list:
		AMR = True
	else:
		AMR = False

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
			age = - ytdataset[field,"particle_birth_epoch"].in_units("Gyr").value

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


