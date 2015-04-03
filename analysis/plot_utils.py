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

from ramses_pp.analysis import read_utils, utils
from ramses_pp.modules import Simulation


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

def plot_position(data,x_field,y_field,z_field,width, depth,filter=None, filter_name=None,type="particle", quantity="particle_mass"):
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

	print "plotting position" + type_string + filter_string + ".png"
	plt.savefig("position" + type_string + filter_string + ".png")
	plt.close()


def plot_profile(data,x_field,y_field,units=["Mpc","km/s"],filter=None,bins=100,weight_profile="ones",abs=False):
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




def plot_projection(data,x_field,y_field,q_field,width, q_unit,filter=None, filter_name=None,type="particle", bins=100, r_bins=10, density=False, density_scale=1000.0):
	"""
	takes in a YT xyz field and plots the position
	in xy xz and yz
	and also does a 3d plot
	assumes the z field has already beein filtered
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
	if density:
    	     val_array = val_array / (freq_array * (x_size * density_scale * y_size * density_scale)) # mass divided by the area of each radial bin
        return r_bins, val_array, 


def plot_sfr(sphere, plot_name="selene", n_bins=50):

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

	time_range = [0.0,14.0] # Giga Years

        #print formation_time, "this"
	hist, bins = np.histogram(formation_time, bins=n_bins, range=time_range,)
	inds = np.digitize(formation_time, bins=bins)

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
		plt.ylabel('SFR')
		plt.savefig(("%s_sfr.png" % plot_name))
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
	
def plot_pdf(x, plotname, x_lab, y_lab="PDF", x_min=-1.0, x_max=1.5, nbins=100):
	"""
	This function plots the PDF of a distribution x
	"""

	# filter x between the min and max values

	filter = np.logical_and(x >= x_min, x <= x_max)
	x = x[filter]
	weights = np.ones_like(x)/len(x)
	plt.hist(ratio,bins=nbins,weights=weights, histtype="step")
	hist, bins = np.histogram(ratio,bins=nbins,weights=weights)

	from scipy.interpolate import UnivariateSpline
	x_vals = bins[:-1] + ((bins[1] - bins[0]) / 2.0 )
	f = UnivariateSpline(x_vals, hist, s=nbins)
	plt.plot(x_vals,f(x_vals))
	plt.xlabel(x_lab)
	plt.ylabel(y_lab)
	plt.axis([x_min, x_max,0.0,(max(hist) + (max(hist) * 0.1))])
	plt.savefig(plotname)
	plt.close()
