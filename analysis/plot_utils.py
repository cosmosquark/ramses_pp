import yt
from yt import derived_field
from ramses_pp.modules import Simulation
from yt.data_objects.particle_filters import add_particle_filter
from ramses_pp.modules.yt.YT import star_filter, dark_filter, young_star_filter
from yt.units.yt_array import YTArray
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from yt.utilities.physical_constants import G
import shelve
from ramses_pp import config

# ramses_pp modules

from ramses_pp.analysis import read_utils
from ramses_pp.modules import Simulation


plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('text', usetex=True)

"""
this code simply takes an input x/y/z field, andd plots things
"""

def subplot_plot(plotno,x,y,lab,title,x_axis,y_axis,alpha=0.5):

	axis_range = x_axis + y_axis
	plt.subplot(2,2,plotno)
	plt.scatter(x,y,s=10,c="r",alpha=alpha)
	plt.title(title)
	plt.xlabel(lab)
	plt.ylabel(lab)
	plt.axis(axis_range)

def plot_position(data,x_field,y_field,z_field,width, depth,filter=None, filter_name=None,type="particle", quantity="particle_mass", galaxy_name=None, dir_name = None):
	"""
	takes in a YT xyz field and plots the position
	in xy xz and yz
	and also does a 3d plot
	"""

	if width[1] != depth[1]:
		print "at least use freaking consistant units"
		return None

	if type != None and type != "particle":
		if filter != None and filter_name != None:
			x = data[type,x_field][filter].in_units(width[1]).value
			y = data[type,y_field][filter].in_units(width[1]).value
			z = data[type,z_field][filter].in_units(width[1]).value
		else:
			x = data[type,x_field].in_units(width[1]).value
			y = data[type,y_field].in_units(width[1]).value
			z = data[type,z_field].in_units(width[1]).value

	else:
		if filter != None and filter_name != None:
			x = data[x_field][filter].in_units(width[1]).value
			y = data[y_field][filter].in_units(width[1]).value
			z = data[z_field][filter].in_units(width[1]).value
		else:
			x = data[x_field].in_units(width[1]).value
			y = data[y_field].in_units(width[1]).value
			z = data[z_field].in_units(width[1]).value

	x_range =  [-width[0],width[0]]
	y_range = [-width[0],width[0]]
	z_range = [-depth[0],depth[0]]

	plt.figure(dpi=1200)

	subplot_plot(0,x,y,"kpc","x-y",x_range, y_range)
	subplot_plot(1,x,z,"kpc","y-z",x_range, z_range)
	subplot_plot(2,y,z,"kpc","x-z",x_range, z_range)

	if type:
		type_string = "_" + str(type)
	else:
		type_string = ""
	if filter != None and filter_name != None:
		filter_string = "_" +  str(filter_name)
	else:
		filter_string = ""

	if galaxy_name:
		galaxy_string = "_" + galaxy_name
	else:
		galaxy_string = ""

	print "plotting position" + type_string + filter_string + ".png"
	if dir_name != None:
		plt.savefig(dir_name + "position" + galaxy_name + type_string + filter_string + ".png")
	else:
		plt.savefig("position" + galaxy_name + type_string + filter_string + ".png")
	plt.close()


def plot_profile(data,x_field,y_field,units=["Mpc","km/s"],filter=None,bins=100,weight_profile="ones",abs=False):
	"""
	get some x data.. and some y data.. then plot it.. simples
	"""
	if filter != None:
		x = data[x_field][filter].in_units(units[0]).value
		y = data[y_field][filter].in_units(units[1]).value
	else:
		x = data[x_field].in_units(units[0]).value
		y = data[y_field].in_units(units[1]).value

	if abs:
		x = np.absolute(x)
		y = np.absolute(y)

	x_range = [x.min(),x.max()]
	y_range = [y.min(),y.max()]

	print "test"
	print x_range, y_range
	nbin = bins + 1

	x_bins_temp, step = np.linspace(x_range[0],x_range[1],nbin, retstep = True)

	x_bins = x_bins_temp + (step / 2.0)
	x_bins = np.delete(x_bins,(len(x_bins)-1)) # pop the last element of the array

	print x_bins
	val_array, xedges = np.histogram(x,bins=x_bins_temp,weights=y)
	print val_array

	if weight_profile=="ones":
		freq_array, fedges = np.histogram(x,bins=x_bins_temp)
	else:
		freq_array, fedges = np.histogram(x,bins=x_bins_temp,weights=weight_profile)

	# normalise
	if weight_profile == None:
		val_array = val_array
	else:
		val_array = val_array / freq_array
	return x_bins, val_array

def plot_manual_profile(x_field,y_field,units=["Mpc","km/s"],filter=None,bins=100,weight_profile="ones",abs=False):
	"""
	works for particles kinda... but super deprecated
	"""	
	if filter != None:
		x = x_field[filter].in_units(units[0]).value
		y = y_field[filter].in_units(units[1]).value
	else:
		x = x_field.in_units(units[0]).value
		y = y_field.in_units(units[1]).value

	if abs:
		x = np.absolute(x)
		y = np.absolute(y)

	x_range = [x.min(),x.max()]
	y_range = [y.min(),y.max()]

	print "test"
	print x_range, y_range

	

	nbin = bins + 1

	x_bins_temp, step = np.linspace(x_range[0],x_range[1],nbin, retstep = True)

	x_bins = x_bins_temp + (step / 2.0)
	x_bins = np.delete(x_bins,(len(x_bins)-1)) # pop the last element of the array

	print x_bins

	val_array, xedges = np.histogram(x,bins=x_bins_temp,weights=y)
	print val_array

	if weight_profile=="ones":
		freq_array, fedges = np.histogram(x,bins=x_bins_temp)
	else:
		freq_array, fedges = np.histogram(x,bins=x_bins_temp,weights=weight_profile)

	# normalise
	if weight_profile == None:
		val_array = val_array
	else:
		val_array = val_array / freq_array
	return x_bins, val_array




def plot_projection_deprecated(data,x_field,y_field,q_field,width, q_unit,filter=None, filter_name=None,type="particle", bins=100, r_bins=10, extra="mean", density_scale=1.0):
	"""
	takes in a YT xyz field and plots the position
	in xy xz and yz
	and also does a 3d plot
	assumes the z field has already beein filtered

	DEPRECATED
	"""


	if type != None and type != "particle":
		if filter != None and filter_name != None:
			x = data[type,x_field][filter].in_units(width.units).value
			y = data[type,y_field][filter].in_units(width.units).value
			q = data[type,q_field][filter].in_units(q_unit).value
		else:
			x = data[type,x_field].in_units(width.units).value
			y = data[type,y_field].in_units(width.units).value
			q = data[type,q_field].in_units(q_unit).value

	else:
		if filter != None and filter_name != None:
			x = data[x_field][filter].in_units(width[1]).value
			y = data[y_field][filter].in_units(width[1]).value
			q = data[q_field][filter].in_units(q_unit).value
		else:
			x = data[x_field].in_units(width[1]).value
			y = data[y_field].in_units(width[1]).value
			q = data[q_field].in_units(q_unit).value

	x_range =  [-width.value,width.value]
	y_range = [-width.value,width.value]

	nbin = bins + 1

	x_bin_temp = np.linspace(x_range[0],x_range[1],nbin)
	y_bin_temp = np.linspace(y_range[0],y_range[1],nbin)
	x_bin_width = (x_bin_temp[1] - x_bin_temp[0]) / 2.0
	y_bin_width = (x_bin_temp[1] - x_bin_temp[0]) / 2.0

	# find the quantity of the weight in each bin
	data_array, xedges, yedges = np.histogram2d(x,y,bins=[x_bin_temp,y_bin_temp],weights=q)

	x_bin_temp = x_bin_temp + x_bin_width
	y_bin_temp = y_bin_temp + y_bin_width

	x_bins = x_bin_temp[:(bins)]
	y_bins = y_bin_temp[:(bins)]

	x_size = abs(x_bins[1] - x_bins[0])
	y_size = abs(y_bins[1] - y_bins[0])

        # multiply this
	x_bins_r = np.power(x_bins,2)
#       print radius_array      
	y_bins_r = np.power(y_bins,2)
	r_array = np.add.outer(x_bins_r,y_bins_r)
	r_array = np.sqrt(r_array)


	r_bins_temp, step = np.linspace(r_array.min(), r_array.max(), num=(r_bins+1), retstep = True)
	r_bins = r_bins_temp + (step / 2.0)
	r_bins = np.delete(r_bins,(len(r_bins)-1)) # pop the last element of the array

	# sum the quantity within each radial bin
	val_array, val_width = np.histogram(r_array,bins=r_bins_temp,weights=data_array)

	# calculate the total number of points in each bin
	freq_array, val_width = np.histogram(r_array,bins=r_bins_temp) # how many pixels per radial bin?


        # density
	if extra == "density":
		# here we care about the physical densitiy/dimensions of everything		
		val_array = val_array / (freq_array * (x_size * density_scale * y_size * density_scale)) # mass divided by the area of each radial bin

	if extra == "mean":
		# here we do not care about the physical size of everything
		val_array = val_array / freq_array

	
        return r_bins, val_array, 


def plot_sfr(sphere, plot_name="selene", n_bins=50, dir_name=None):

        #TODO improve this with the more recent histogram plotting

        # now compute the SFR (mass formed at epoch time / bin size in years)

        # warning of ambiguity
        # particle_age = coformal time
        # particle_birth_epoch = cosmic time things are born
        # particle_total_time ... don't use this... related to feedback

	masses = sphere["stars","particle_mass"].in_units("Msun").value
	formation_time = - sphere["stars","particle_birth_epoch"].in_units("Gyr").value
	metallicity = sphere["stars","particle_metallicity"].value
	age = -sphere["stars","particle_birth_epoch"].in_units("Gyr").value
#	age_max = age.max()
#	age_min = age.min()
#	age = abs(age - age_max) + age_min 

	time_range = [0.0,14.0] # Giga Years

        #print formation_time, "this"
	hist, bins = np.histogram(formation_time, bins=n_bins, range=time_range,)
	inds = np.digitize(formation_time, bins=bins)
	# include the max values in the last bin
	inds[np.argmax(inds)] = inds.max() - 1

        #ime = max(time) - time
	time = (bins[:-1] + bins[1:]) / 2.0
        #time = time[:i1]


	sfr = np.array([masses[inds == j].sum()/((bins[j+1]-bins[j])*(10**9)) for j in range(len(time))])
	sfr[sfr == 0] = np.nan
	
	time = max(time) - time  # is this not a convoluted time = time[:-1]  quick and dirty.. really we want to be calling the year from the snapshot time

	print time
	print sfr

	if plot_name != None:
	        ## age abundance
		plt.plot(time, sfr)
		plt.xlabel('Time  [Gyr]')
		plt.ylabel('SFR [Msun/Year]')
		if dir_name != None:
			plot_name = dir_name + plot_name
		plt.savefig(("%s_sfr.png" % plot_name))
		plt.savefig(("%s_sfr.pdf" % plot_name))	
		plt.close()

	return time, sfr


def flatten_line(x,y,extra_x=None,no_nan=True,append_max=True, double_zero = True):
	"""
	This will flatten a curve to 0.0 as of when needed.
	It will also append an extra x value if needed.
	in which the corresponding y value for this will be zero.
	Mainly used to make plots look a bit more pretty
	"""

	if no_nan:
		# this kills any curve from the first NaN, inf etc
		idx = np.argmax(y, axis=None)
		y[idx:] = 0.0
	else:
		idx = len(y) - 1

	# lets set the next bin element to be zero
	bin_width = (x[1] - x[0]) / 2.0
	if append_max:
		new_element = x.max() + bin_width
		x = np.append(x,new_element)
		y = np.append(y,0.0)

	if extra_x != None:
		x = np.append(x,extra_x)
		y = np.append(y,0.0)

	if double_zero == True:
		# check if the first element is (0,0)
		if x[0] != 0 and y[0] != 0:
			x = np.insert(x,0,0.0)
			y = np.insert(y,0,0.0)
		# maybe the length of X is greater than Y.. 
		# if so, then we already have the extra value in X, so just append Y
		if x[0] == 0.0 and len(x) == len(y) + 1:
			y = np.insert(y,0,0.0)


		
	return x, y

def kill_first_nan(x,y):
	"""
	This will delete array values from and beyond the first NaN found
	"""

	# TODO have this cancel the loop and just delete from the minimum value found
	# instead of looping over each individual element.
	delete = False
	for i in range(0,len(x)):
		if delete == True:
			np.delete(x,i)
			np.delete(y,i)

		if x[i] == np.nan:
			np.delete(x,i)
			np.delete(y,i)
			delete = True
			
	return x, y
	
def plot_pdf(x, plotname, x_lab, y_lab="PDF", x_min=-1.0, x_max=1.5, nbins=100, spline=False):
	"""
	This function plots the PDF of a distribution x, usually for Jz/Jcirc cuts
	"""

	# filter x between the min and max values

	filter = np.logical_and(x >= x_min, x <= x_max)
	x = x[filter]
	weights = np.ones_like(x)/len(x)
	plt.hist(x,bins=nbins,weights=weights, histtype="step")
	hist, bins = np.histogram(x,bins=nbins,weights=weights)

	if spline:
		from scipy.interpolate import UnivariateSpline
		x_vals = bins[:-1] + ((bins[1] - bins[0]) / 2.0 )
		f = UnivariateSpline(x_vals, hist, s=nbins)
		plt.plot(x_vals,f(x_vals))
	plt.xlabel(x_lab)
	plt.ylabel(y_lab)
	plt.axis([x_min, x_max,0.0,(max(hist) + (max(hist) * 0.1))])
	plt.savefig(plotname)
	plt.close()


def plot_standard_vcirc(sphere,cylinder,snap,galaxy_name,r_min=None,r_max=None,dr_factor=2,vrot=True):
	"""
	Your bog standard Vcirc plot, with the option to append rotation velocities to the plot.
	Requires stars, dark matter to be pre-filtered first
	"""

	from ramses_pp.analysis import utils, filter_utils

	if r_max == None:
		r_max = sphere["particle_positions_cylindrical_radius"].max().in_units("kpc")
#		r_max = sphere.ds.arr(50,"kpc")

	if r_min == None:
		r_min = sphere.ds.arr(0,"kpc")

	# 1) get radial values
	r_indices, r_bins, r_filter, r_truths = utils.radial_indices('particle_position_spherical_radius',sphere,r_min=None,r_max=r_max,dr_factor=2)


	ytsnap = snap.raw_snapshot()
	n_bins = len(r_bins)
	r_bins = ytsnap.arr(r_bins,"cmcm")
	m_bins_all = utils.mass_enclosed_bins(sphere,snap,r_bins,shape="sphere",type="all")
	m_bins_dark = utils.mass_enclosed_bins(sphere,snap,r_bins,shape="sphere",type="dark")
	m_bins_star = utils.mass_enclosed_bins(sphere,snap,r_bins,shape="sphere",type="stars")
	m_bins_gas = utils.mass_enclosed_bins(sphere,snap,r_bins,shape="sphere",type="gas")

	vcirc_all = utils.vcirc(r_bins,m_bins_all,sphere)
	vcirc_dark = utils.vcirc(r_bins,m_bins_dark,sphere)
	vcirc_star = utils.vcirc(r_bins,m_bins_star,sphere)
	vcirc_gas = utils.vcirc(r_bins,m_bins_gas,sphere)

	plt.plot(r_bins.in_units("kpc"),vcirc_all.in_units("km/s"),color="black",label=r"$V_{circ,all}$")
	plt.plot(r_bins.in_units("kpc"),vcirc_dark.in_units("km/s"),color="green",label=r"$V_{circ,dark}$")
	plt.plot(r_bins.in_units("kpc"),vcirc_star.in_units("km/s"),color="blue",label=r"$V_{circ,stars}$")
	plt.plot(r_bins.in_units("kpc"),vcirc_gas.in_units("km/s"),color="red",label=r"$V_{circ,gas}$")

	plt.xlabel("Radius [kpc]")
	plt.ylabel("Rotation Speed [km/s]")


	if vrot:
		# THE DISK SHOULD ALREADY BE FILTERED SPATIALLY KINDA
		print "doing stuff with vrot"
		min_z = sphere.ds.arr(-3.0,"kpc")
		max_z = sphere.ds.arr(3.0,"kpc")
	
		cold_disk = cylinder.cut_region(["obj['temperature'] < 1e4"])
		
#		L = sphere.quantities.angular_momentum_vector(use_gas=True,use_particles=False)
#		bv = sphere.get_field_parameter("bulk_velocity")
#		L_mag = np.sqrt(np.power(L[0],2.0) + np.power(L[1],2.0) + np.power(L[2],2.0))
#		L_norm = np.zeros(3)
#		L_norm[0] = L[0].value/L_mag
#		L_norm[1] = L[1].value/L_mag
#		L_norm[2] = L[2].value/L_mag
		
#		sphere.set_field_parameter("bulk_velocity",bv)
#		sphere.set_field_parameter("normal",sphere.ds.arr(L_norm,"code_length"))

#		cold_sphere.set_field_parameter("bulk_velocity",bv)
#		cold_sphere.set_field_parameter("normal",sphere.ds.arr(L_norm,"code_length"))
		

		add_particle_filter("young_stars", function=young_star_filter, filtered_type="all", requires=["particle_age"])
		cylinder.ds.add_particle_filter("young_stars")


		# cold gas
		
		print "filtering spatially"
		spatial_filter_tot = filter_utils.min_max_filter(cylinder,"particle_position_relative_z",min_z,max_z)
		spatial_filter_dark = filter_utils.min_max_filter(cylinder,('dark','particle_position_relative_z'),min_z,max_z)
		spatial_filter_stars = filter_utils.min_max_filter(cylinder,('stars','particle_position_relative_z'),min_z,max_z)
		spatial_filter_young_stars = filter_utils.min_max_filter(cylinder,('young_stars','particle_position_relative_z'),min_z,max_z)

		r_tot = np.sqrt(np.power(cylinder["particle_position_relative_x"][spatial_filter_tot],2) + np.power(cylinder["particle_position_relative_y"][spatial_filter_tot],2))
		r_stars = np.sqrt(np.power(cylinder["stars","particle_position_relative_x"][spatial_filter_stars],2) +  np.power(cylinder["stars","particle_position_relative_y"][spatial_filter_stars],2))
		r_young_stars = np.sqrt(np.power(cylinder["young_stars","particle_position_relative_x"][spatial_filter_young_stars],2) +  np.power(cylinder["young_stars","particle_position_relative_y"][spatial_filter_young_stars],2))
		r_dark = np.sqrt(np.power(cylinder["dark","particle_position_relative_x"][spatial_filter_dark],2) + np.power(cylinder["dark","particle_position_relative_y"][spatial_filter_dark],2))
		
#		r_gas = sphere["cylindrical_r"]
#		r_cold_gas = cold_sphere["cylindrical_r"]

		print "computing vrots"
#		vrot_tot = utils.manual_vrot(cylinder,type="all",r_filter = spatial_filter_tot)

#		vrot_stars = utils.manual_vrot(cylinder,type="stars",r_filter = spatial_filter_stars)
#		vrot_young_stars = utils.manual_vrot(cylinder,type="young_stars",r_filter = spatial_filter_young_stars)
#		vrot_dark = utils.manual_vrot(cylinder,type="dark",r_filter = spatial_filter_dark)

		vrot_stars = yt.ProfilePlot(cylinder,('stars','particle_position_cylindrical_radius'),[('stars','particle_velocity_cylindrical_theta')],weight_field=('stars','particle_mass'),n_bins=n_bins,x_log=False,y_log={('stars','particle_velocity_cylindrical_theta'):False})		
		vrot_young_stars = yt.ProfilePlot(cylinder,('young_stars','particle_position_cylindrical_radius'),[('young_stars','particle_velocity_cylindrical_theta')],weight_field=('young_stars','particle_mass'),n_bins=n_bins,x_log=False,y_log={('young_stars','particle_velocity_cylindrical_theta'):False})		
		vrot_dark = yt.ProfilePlot(cylinder,('dark','particle_position_cylindrical_radius'),[('dark','particle_velocity_cylindrical_theta')],weight_field=('dark','particle_mass'),n_bins=n_bins,x_log=False,y_log={('dark','particle_velocity_cylindrical_theta'):False})		
		
		r_stars_binned = vrot_stars.profiles[0].x.in_units("kpc").value
		r_young_stars_binned = vrot_young_stars.profiles[0].x.in_units("kpc").value
		vrot_stars_binned = np.abs(vrot_stars.profiles[0][('stars','particle_velocity_cylindrical_theta')]).in_units("km/s").value
		vrot_young_stars_binned = np.abs(vrot_young_stars.profiles[0][('young_stars','particle_velocity_cylindrical_theta')]).in_units("km/s").value
		


		# this is not filtering in height yet

		vrot_gas = yt.ProfilePlot(cylinder,'cylindrical_r',["velocity_cylindrical_theta"],weight_field="cell_mass",n_bins=n_bins,x_log=False,y_log={"velocity_cylindrical_theta":False})
		r_gas = vrot_gas.profiles[0].x
		vrot_gas = np.abs(vrot_gas.profiles[0]["velocity_cylindrical_theta"]) # this overwrites the yt plot object


		vrot_cold_gas = yt.ProfilePlot(cold_disk,'cylindrical_r',["velocity_cylindrical_theta"],weight_field="cell_mass",n_bins=n_bins,x_log=False,y_log={"velocity_cylindrical_theta":False})
		r_cold_gas = vrot_cold_gas.profiles[0].x
		vrot_cold_gas = np.abs(vrot_cold_gas.profiles[0]["velocity_cylindrical_theta"])

		# bin the particle data
		print "binning the data"
	#	r_tot_binned, vrot_tot_binned = plot_manual_profile(r_tot,vrot_tot,units=["kpc","km/s"],filter=None,bins=30,abs=True,weight_profile="ones")
	#	r_stars_binned, vrot_stars_binned = plot_manual_profile(r_stars,vrot_stars,units=["kpc","km/s"],filter=None,bins=30,abs=True,weight_profile="ones")
	#	r_dark_binned, vrot_dark_binned = plot_manual_profile(r_dark,vrot_dark,units=["kpc","km/s"],filter=None,bins=30,abs=True,weight_profile="ones")
	#	r_young_stars_binned, vrot_young_stars_binned = plot_manual_profile(r_young_stars,vrot_young_stars,units=["kpc","km/s"],filter=None,bins=30,abs=True,weight_profile="ones")

		# make the plots look a little more pretty
		r_gas, vrot_gas = flatten_line(r_gas.in_units("kpc").value,vrot_gas.in_units("km/s").value,no_nan=False,extra_x=r_bins.in_units("kpc").max())
		# make the plots look a little more pretty
		r_cold_gas, vrot_cold_gas = flatten_line(r_cold_gas.in_units("kpc").value,vrot_cold_gas.in_units("km/s").value,no_nan=False,extra_x=r_bins.in_units("kpc").max())
		
		r_stars_binned, vrot_stars_binned = flatten_line(r_stars_binned,vrot_stars_binned,no_nan=True,append_max=True,double_zero=True,extra_x = r_bins.in_units("kpc").max())
		r_young_stars_binned, vrot_young_stars_binned = flatten_line(r_young_stars_binned,vrot_young_stars_binned,no_nan=True,append_max=True,double_zero=True,extra_x = r_bins.in_units("kpc").max())
		
		# finally do the plots
		print "plotting the data"
#		plt.plot(r_dark,vrot_dark,label="vrot dark")
		plt.plot(r_stars_binned,vrot_stars_binned,label=r"$V_{rot,stars}$",color="blue",linestyle='dashed')
		plt.plot(r_young_stars_binned,vrot_young_stars_binned,label=r"$V_{rot,young stars}$",color="purple",linestyle='dashed')
		plt.plot(r_gas, vrot_gas, label=r"$V_{rot,gas}$",color="red",linestyle='dashed')
		plt.plot(r_cold_gas, vrot_cold_gas, label=r"$V_{rot,cold gas}$",color="orange",linestyle='dashed')

	plt.legend()
	plt.savefig(("%s_rot_circ.png" % galaxy_name))
	plt.savefig(("%s_rot_circ.pdf" % galaxy_name))
	plt.close()

def plot_xy(x,y,x_lab,y_lab,filename,type="line",x_min=None,x_max=None,y_min=None,y_max=None,label=None,y_err=None):
	"""
	Use this function to plot x vs y super lazily
	"""

	#TODO make it so this function accepts a list of x and y 
	#TODO get min and max working nicely
	

	if type == "line":
		plt.plot(x,y)
	elif type == "bar":
		plt.bar(x,y)
	elif type == "step":
		plt.step(x,y,where="mid")

	plt.xlabel(x_lab)
	plt.ylabel(y_lab)

	if y_err != None:
		plt.errorbar(x,y,yerr=y_err)

	if x_min != None and x_max != None:
		plt.xlim([x_min,x_max])

	if y_min != None and y_max != None:
		plt.ylim([y_min,y_max]) 

	print filename
	plt.savefig(("%s.png" % filename))
	plt.close()
	


	return


def plot_frb_profile(image,width,y_units,x_lab,y_lab,filename,n_bins=50,ylog=True,extrema=[None,None,None,None]):
	"""
	Use this method to plot the 2d radial profile of an FRB image generated from YT
	"""


	n_bins = image.shape[0]
	bin_width = width.v / float(n_bins)
	left = - width.v / 2.0
	right = width.v / 2.0
	length = np.arange(left,right,bin_width)
	shift_thing = (length[2] - length[1]) / 2.0
	length = length + shift_thing

	# to make radial array, we need to calculate the radial distance of each pixel in a 2d array. Namely.

	length_squared = np.power(length,2.0) # power of 2
	radius_squared = np.add.outer(length_squared, length_squared) # adds as an outer product
	radius = np.sqrt(radius_squared)

	# now we need to generate some radial bins to histogram.

	radius_bin_width = radius.max() / float(n_bins)
	radius_bins = np.arange(0,radius.max(), radius_bin_width)

	# flatten the arrays

	flat_radius = radius.flatten()
	flat_data = np.array(image.in_units(y_units)).flatten()

	n, _ = np.histogram(flat_radius, bins=radius_bins)
	sy, _ = np.histogram(flat_radius, bins=radius_bins, weights=flat_data)
	sy2, _ = np.histogram(flat_radius, bins=radius_bins, weights=flat_data*flat_data)
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)

	radius_bins = radius_bins + (radius_bin_width / 2.0)
	radius_bins = np.delete(radius_bins,(len(radius_bins)-1))


	plt.plot(radius_bins,mean)
	if ylog == True:
		plt.yscale('log')
	plt.xlabel(x_lab)
	plt.ylabel(y_lab)

	if extrema[0] != None and extrema[1] != None:
		plt.xlim(extrema[0],extrema[1])
	
	if extrema[2] != None and extrema[3] != None:
		plt.ylim(extrema[2],extrema[3])


	plt.savefig(("%s.png" % filename))
	plt.close()

def plot_particle_projection(x_field, y_field, z_field, weight_field="mean", density=False):
	"""
	This method takes in a series of particle data in a 3D volume
	Then projects it into a 2D image
	And then can be followed on with some further routines

	Data must be in YT arrays

	if weight field = "mean"... it will just calculate the mean
	if density = True, it will also take into account the physical dimensions of each bin too.
	"""
	return None

	







