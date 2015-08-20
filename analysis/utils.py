###
# A list of functions, routines etc that do very useful things 
# with ramses_pp and YT data
#
#
#
#
#

import yt
from yt import derived_field
from ramses_pp.modules import Simulation
from yt.data_objects.particle_filters import add_particle_filter
from ramses_pp.modules.yt.YT import star_filter, dark_filter, young_star_filter
from yt.units.yt_array import YTArray
import numpy as np
from matplotlib import pyplot as plt
from yt.utilities.physical_constants import G
import shelve
from ramses_pp.analysis import read_utils, filter_utils, plot_utils
import scipy as sp


import weakref, copy
import abc
import re, uuid
import abc


def compute_star_met(object, ref=None, type="stars"):
	"""
	Computes the metallicities from either a YT dataset
	Or a dictionary of YT Arrays
	"""


	## Anders & Grevesse 1989 as default values
	sol_abund = {"H":0.706, "He":0.275, "C":3.03e-3, "N":1.11e-3, "O":9.59e-3, "Ne":0.00000001, "Mg":5.15e-4, "Si":6.53e-4, "Fe":1.17e-3,"Z":0.019}

	if ref=="asplund":
		sol_abund = {"H":0.715, "He":0.270, "C": 0.0024, "N":0.000728 , "O":0.006, "Ne":0.0013, "Mg":0.000742, "Si":0.0007, "Fe":0.00135,  "Z":0.0143}

	if type != None:

		FeH = np.log10(object[type,"particle_Fe_fraction"].value / object[type,"particle_H_fraction"].value) - np.log10(sol_abund["Fe"] / sol_abund["H"])
		OFe = np.log10(object[type,"particle_O_fraction"].value / object[type,"particle_Fe_fraction"].value) - np.log10(sol_abund["O"] / sol_abund["Fe"])
		MgFe = np.log10(object[type,"particle_Mg_fraction"].value / object[type,"particle_Fe_fraction"].value) - np.log10(sol_abund["Mg"] / sol_abund["Fe"])
		CFe = np.log10(object[type,"particle_C_fraction"].value / object[type,"particle_Fe_fraction"].value) - np.log10(sol_abund["C"] / sol_abund["Fe"])
		NFe = np.log10(object[type,"particle_N_fraction"].value / object[type,"particle_Fe_fraction"].value) - np.log10(sol_abund["N"] / sol_abund["Fe"])
		NeFe = np.log10(object[type,"particle_Ne_fraction"].value / object[type,"particle_Fe_fraction"].value) - np.log10(sol_abund["Ne"] / sol_abund["Fe"])
		SiFe = np.log10(object[type,"particle_Si_fraction"].value / object[type,"particle_Fe_fraction"].value) - np.log10(sol_abund["Si"] / sol_abund["Fe"])
	
	else:

		FeH = np.log10(object["particle_Fe_fraction"].value / object["particle_H_fraction"].value) - np.log10(sol_abund["Fe"] / sol_abund["H"])
		OFe = np.log10(object["particle_O_fraction"].value / object["particle_Fe_fraction"].value) - np.log10(sol_abund["O"] / sol_abund["Fe"])
		MgFe = np.log10(object["particle_Mg_fraction"].value / object["particle_Fe_fraction"].value) - np.log10(sol_abund["Mg"] / sol_abund["Fe"])
		CFe = np.log10(object["particle_C_fraction"].value / object["particle_Fe_fraction"].value) - np.log10(sol_abund["C"] / sol_abund["Fe"])
		NFe = np.log10(object["particle_N_fraction"].value / object["particle_Fe_fraction"].value) - np.log10(sol_abund["N"] / sol_abund["Fe"])
		NeFe = np.log10(object["particle_Ne_fraction"].value / object["particle_Fe_fraction"].value) - np.log10(sol_abund["Ne"] / sol_abund["Fe"])
		SiFe = np.log10(object["particle_Si_fraction"].value / object["particle_Fe_fraction"].value) - np.log10(sol_abund["Si"] / sol_abund["Fe"])


	data = {
		"FeH": FeH,
		"OFe": OFe,
		"MgFe": MgFe,
		"CFe": CFe,
		"NFe": NFe,
		"NeFe": NeFe,
		"SiFe": SiFe,
	}

	return data
	
def radial_indices(field,data,r_min=None,r_max=None,dr_factor=1):
	"""
	use this function to create radial bins of particles based on the smallest possible cell size
	and bin increments of the smallest cell size
	of a radial field of your choice

	field = e.g data["stars","particle_spherical_position_radius"]
	data = a data object so you can do data.ds
	correction = a slight correction so that sphere too small errors do not appear

	the dr_factor will allow you to reduce the number of bins
	a dr_factor of 2 will produce radial bins for min_radii x2 bin widths
	"""


	r_min_used = True
	r_max_used = True
	if r_min == None:
		r_min_used = False
		if data[field].min() > data.ds.index.get_smallest_dx():
			r_min = data[field].min() 
		else:
			r_min = data.ds.index.get_smallest_dx() #+ data.ds.arr(0.0001,"code_length")
	else:
		r_min_temp = r_min
		if r_min_temp.in_units("code_length").value > data.ds.index.get_smallest_dx().in_units("code_length").value:
			r_min = r_min_temp
		else:
			r_min = data.ds.index.get_smallest_dx() #+ data.ds.arr(0.0001,"code_length")
	if r_max == None:
		r_max_used = False
		r_max = data[field].max()
	
	# filter in the selection range
	if isinstance(field,tuple):
		filter_r = field[1]
		type = field[0]
	else:
		filter_r = field
		type = None

	if r_min_used == True or r_max_used == True:
		r_truths = filter_utils.min_max_truths(data,filter_r,r_min,r_max,type=type)
		r_filter = filter_utils.min_max_filter(data,filter_r,r_min,r_max,extra_filter=None,type=type)

	else:
		# basically returns a length of arrays that
		r_truths = np.ones(len(data[field]), dtype=bool)
		r_filter = np.arange(0,len(data[field]))

		# with that filter added, this array is not going to be the same length as the origonal data
		# hence the need to return the filter
		# r_indices is of the same length as the data array *after* applying r_filter,
		# r_indicies contains the bin indexes of the filtered field for whatever quantity that field should belong to 
		# r_truths retains the context of the origional array size and states whether a field is filtered or not

		# tl'dr 
		# r_indices contains the bin indicies of each filtered field for which property in y belongs to bin x
		# e.g the particle indicies belonging to each bin
		# f_filter contains additional filter
		# r_truths is what you need to apply to the origional array and maintains the same size as the origional array


	print "len of the truths and filter arrays"
	print len(r_truths), len(r_filter)
	r_indices = np.floor(data[field][r_filter].in_units("code_length").value / (dr_factor * r_min.in_units("code_length").value)).astype(int)
	r_bins = data.ds.arr(np.zeros(int(r_indices.max() + 1)),"code_length")
	for i in range(0,len(r_bins)):
		r_bins[i] = ( (i * dr_factor) + 1) * r_min.in_units("code_length")
	try:
		r_bins = r_bins.in_units("cmcm")
	except yt.units.unit_object.UnitParseError:
		r_bins = r_bins.in_units("cm")
	return r_indices, r_bins, r_filter, r_truths


def mass_enclosed_bins(object,snap,r_bins,shape="sphere",type="all",sim_object_type=None, AMR=True):
	"""
	extending radial_indices, this function can find the enclosed mass of
	of a particular field (stars, dm, gas, total or custom)
	returns the m_bins to accompany the r_bins

	tl;dr you feed this r_bins, and you get m_bins out
	"""

	# overide the need for snap


	m_bins = object.ds.arr(np.zeros(len(r_bins)),"g")
	invalid = False
	if sim_object_type == None:
		print "ramses_pp snapshot"
		dataset = snap.raw_snapshot()
	elif sim_object_type == "ad":
		print "ytsnapshot"
		dataset = snap.ds
	else:
		print "using object as snapshot"
		dataset = object.ds

	try:
		test = r_bins.in_units("cmcm")
		min_unit = "cmcm"
	except yt.units.unit_object.UnitParseError:
		min_unit = "cm"


	for i in range(0,len(r_bins)):
		if type == "all":
			#print object.ds.field_list
			if AMR == True:
				try:
					m_bins[i] = dataset.sphere(object.center,object.ds.arr(r_bins[i],min_unit)).quantities.total_quantity(["particle_mass"]) + dataset.sphere(object.center,object.ds.arr(r_bins[i],min_unit)).quantities.total_quantity(["cell_mass"])
				except:
					m_bins[i] = object.ds.arr(0.0,"g")
			else:
				try:
					temp = dataset.sphere(object.center,object.ds.arr(r_bins[i],min_unit))
					m_bins[i] = temp['all','particle_mass'].sum()
				except:
					m_bins[i] = object.ds.arr(0.0,"g")
		elif type == "gas":
			if AMR == True:
				try:
					m_bins[i] = dataset.sphere(object.center,object.ds.arr(r_bins[i],min_unit)).quantities.total_quantity(["cell_mass"])
				except:
					m_bins[i] = object.ds.arr(0.0,"g")
			else:
				try:
					temp = dataset.sphere(object.center,object.ds.arr(r_bins[i],min_unit))
					m_bins[i] = temp['gas','particle_mass'].sum()
				except:
					m_bins[i] = object.ds.arr(0.0,"g")
		elif type == "stars":
			try:
				temp = dataset.sphere(object.center,object.ds.arr(r_bins[i],min_unit))
				m_bins[i] = temp['stars','particle_mass'].sum()
			except:
				m_bins[i] = object.ds.arr(0.0,"g")
		elif type == "dark":
			try:
				temp = dataset.sphere(object.center,object.ds.arr(r_bins[i],min_unit))
				m_bins[i] = temp['dark','particle_mass'].sum()
			except:
				m_bins[i] = object.ds.arr(0.0,"g")
		else:
			invalid = True
			break
	if invalid == True:
		return None
	else:
		return m_bins


def vcirc(r_bins,m_bins,data):
	"""
	computes vcirc for the enclosed material from a given m_bins
	"""

	vcirc = data.ds.arr(np.zeros(len(m_bins)))
	vcirc = (np.sqrt(G * (m_bins) / r_bins))
	return vcirc

def jcirc_bins(r_bins,v_circ,m_bins,data):
	""" this computes j_circ for each individual r_bin, v_circ bin and m_bin"""
	j_circ = data.ds.arr(np.zeros(len(m_bins)))
	j_circ = r_bins * v_circ
	return j_circ

def jcirc(data, r_indices,v_circ,r_filter,type="all", AMR = True):
	""" this computes j_circ for each individual r_bin, v_circ bin and m_bin"""
	if AMR == True and type=="gas":
		print r_filter
		print data["index","spherical_radius"]
		print r_indices
		print v_circ[r_indices]
		print r_filter.shape, data["index","spherical_radius"].shape
		print data["index","spherical_radius"][r_filter] 
		if r_filter != None:
			j_circ = (data["index","spherical_radius"][r_filter] * (v_circ[r_indices].astype(int)))
		else:
			j_circ = (data["index","spherical_radius"] * (v_circ[r_indices].astype(int)))
	else:
		if r_filter != None:
			j_circ = (data[type, "particle_position_spherical_radius"][r_filter] * (v_circ[r_indices].astype(int)))
		else:
			j_circ = (data[type, "particle_position_spherical_radius"] * (v_circ[r_indices].astype(int)))
	# j_circ is of the length of data[field][filter]

	return j_circ



def vrot(data,type="all",r_filter = None):
	"""
	returns vrot (it's the spherical theta componant
	"""

	# this is a bit of a hack for now
	# but if there is a bulk velocity, with the way YT caches the fields
	# we want to recompute the field with the bulk velocity applied

	if data.has_field_parameter("bulk_velocity"):
		data.field_data.pop((type, 'particle_velocity_spherical_radius'))
		data.field_data.pop((type, 'particle_velocity_spherical_theta'))

	if r_filter != None:
		vrot = data[type,"particle_velocity_cylindrical_theta"][r_filter]
		rrot = data[type,"particle_position_cylindrical_radius"][r_filter]
		
	else:
		vrot = data[type,"particle_velocity_cylindrical_theta"]
		rrot = data[type,"particle_position_cylindrical_radius"]
	return vrot

def manual_vrot(data,type="all",r_filter = None,AMR = True):

	# vx * y - yv * x
#	if data.has_field_parameter("bulk_velocity"):
#		data.field_data.pop((type, 'particle_position_relative_x'))
#		data.field_data.pop((type, 'particle_position_relative_y'))
#		data.field_data.pop((type, 'particle_velocity_relative_x'))
#		data.field_data.pop((type, 'particle_velocity_relative_y'))

	if r_filter != None:
		r = np.sqrt(np.power(data[type,"particle_position_relative_x"][r_filter],2) + np.power(data[type,"particle_position_relative_y"][r_filter],2))
		vrot = - (data[type,"particle_position_relative_y"][r_filter] * data[type,"particle_velocity_relative_x"][r_filter])
		vrot = vrot + (data[type,"particle_position_relative_x"][r_filter] * data[type,"particle_velocity_relative_y"][r_filter])
		vrot = vrot / r


	else:
		r = np.sqrt(np.power(data[type,"particle_position_relative_x"],2) + np.power(data[type,"particle_position_relative_y"],2))
		vrot = (data[type,"particle_position_relative_y"] * data[type,"particle_velocity_relative_x"])
		vrot = vrot - (data[type,"particle_position_relative_x"] * data[type,"particle_velocity_relative_y"])
		vrot = vrot / r

	return vrot


def jz(data, type="all", r_filter = None, reverse_orientation = False,AMR = True):
	"""
	returns jz (its the specific angular momentum in the z axis
	assuming your object has a centre and a normal vector
	"""
	try:
		if data.has_field_parameter("bulk_velocity"):
			data.field_data.pop((type, 'particle_specific_angular_momentum_z'))
	except:
		pass
		

	try:
		test = data[type,"particle_specific_angular_momentum_z"].in_units("kpcc**2/s")
		min_unit = "kpccm**2/s"
	except yt.units.unit_object.UnitParseError:
		min_unit = "kpc**2/s"


	if r_filter != None:
		if reverse_orientation == True:
			jz = -data[type,"particle_specific_angular_momentum_z"][r_filter].in_units(min_unit)
		else:
			jz = data[type,"particle_specific_angular_momentum_z"][r_filter].in_units(min_unit)
	else:
		if reverse_orientation == True:
			jz = -data[type,"particle_specific_angular_momentum_z"].in_units(min_unit)
		else:
			jz = data[type,"particle_specific_angular_momentum_z"].in_units(min_unit)
	return jz

def manual_jz(data,type="all",r_filter = None,AMR = True):

	# vx * y - yv * x
#	if data.has_field_parameter("bulk_velocity"):
#		data.field_data.pop((type, 'particle_position_relative_x'))
#		data.field_data.pop((type, 'particle_position_relative_y'))
#		data.field_data.pop((type, 'particle_velocity_relative_x'))
#		data.field_data.pop((type, 'particle_velocity_relative_y'))

	if AMR == True and type == "gas":
		from yt.utilities.math_utils import modify_reference_frame
		from yt.units.yt_array import ucross

		if data.has_field_parameter("normal"):
			normal = data.get_field_parameter("normal")
		else:
			normal = data.ds.arr([0.0,0.0,1.0],"code_length") # default to simulation axis
		bv = data.get_field_parameter("bulk_velocity")
		pos = data.ds.arr([data["%s" % ax] for ax in "xyz"]).T
		vel = data.ds.arr([data["gas", "velocity_%s" % ax] - bv[iax] for iax, ax in enumerate("xyz")]).T
		center = data.get_field_parameter('center')
		L, r_vec, v_vec = modify_reference_frame(center, normal, P=pos, V=vel)
		ang_momentum = -ucross(r_vec, v_vec, registry=data.ds.unit_registry)
		jz = ang_momentum[:, "2"]
		if r_filter != None:
			jz = jz[r_filter]

	else:
		if r_filter != None:
			vrot = (data[type,"particle_position_relative_x"][r_filter] * data[type,"particle_velocity_relative_y"][r_filter])
			vrot = vrot - (data[type,"particle_position_relative_y"][r_filter] * data[type,"particle_velocity_relative_x"][r_filter])
		else:
			vrot = (data[type,"particle_position_relative_x"] * data[type,"particle_velocity_relative_y"])
			vrot = vrot - (data[type,"particle_position_relative_y"] * data[type,"particle_velocity_relative_x"])

		jz = vrot
	return jz


def decomp_stars(data,snap,disk_min=0.8, disk_max=1.1,plot=False,r_min=None,r_max=None,spline=False,sim_object_type=None, manual_jz_flag = False, reverse_orientation=True, AMR = True, type="stars"):
	"""
	At the moment, this technique can only filter for disk stars
	This is computed by calculating the ratio of Jz/Jcirc
	See this paper http://arxiv.org/pdf/1210.2578.pdf
	returns the indices, id's and the ratio
	"""
	print "starting decomp stars"
	try:
		data["stars","particle_index"]
	except: # stars field does not exist, lets add it
		print "adding stars field to data"
		add_particle_filter("stars", function=star_filter, filtered_type="all", requires=["particle_age"])
       		data.ds.add_particle_filter("stars")
	print "decomposing stars into disk and non-disk, please wait"
	if AMR == True and type=="gas":
		print "step on the gas"
		print "field len", len(data["index","spherical_radius"])
		r_indices, r_bins, r_filter, r_truths = radial_indices(("index","spherical_radius"),data,r_min,r_max)
	else:
		print "stars"
		print "field len", len(data[type,"particle_spherical_position_radius"])
		r_indices, r_bins, r_filter, r_truths = radial_indices((type,"particle_spherical_position_radius"),data,r_min,r_max)
	print "LENGTHS"
	print len(r_indices), len(r_bins), len(r_filter), len(r_truths)
	m_bins = mass_enclosed_bins(data,snap,r_bins,shape="sphere",type="all",sim_object_type=sim_object_type, AMR = AMR)
	v_circ = vcirc(r_bins,m_bins,data)
	j_circ = jcirc(data,r_indices,v_circ,r_filter,type=type, AMR = AMR)

	print "FLAGS", manual_jz_flag, reverse_orientation

	if manual_jz_flag == False:
		j_z = jz(data,r_filter=r_filter,type=type, reverse_orientation=reverse_orientation, AMR=AMR)
	else:
		j_z = manual_jz(data,r_filter=r_filter,type=type, AMR=AMR)	


	try:
		test = r_bins.in_units("cmcm")
		min_unit = "cmcm"
	except yt.units.unit_object.UnitParseError:
		min_unit = "cm"

	print "test lengths"
	print len(j_z), "j_z"
	print len(j_circ), "j_circ"
	print len(r_truths), "r_truths"
	print len(r_indices), "r_indices"

	if min_unit == "cmcm":
		ratio = data.ds.arr(j_z.in_units("kpccm**2/s").value / j_circ.in_units("kpccm**2/s").value, "dimensionless")
	else:
		ratio = data.ds.arr(j_z.in_units("kpc**2/s").value / j_circ.in_units("kpc**2/s").value, "dimensionless")
	# get indices of stars which satisfy the conditions

	# expand ratio into context of origonal array
	if AMR == True and type=="gas":
		ratio_expanded = np.zeros((len(data["gas","cell_mass"])))	
	else:
		ratio_expanded = np.zeros((len(data[type,"particle_spherical_position_radius"])))
	print len(ratio_expanded), "ratio_expanded", AMR, type
	j = 0
	for i in range(0, len(r_truths)):
		# if false, insert a dummy value that will never get through the filter
		if r_truths[i] == False:
			ratio_expanded[i] = -99999
		else:
			ratio_expanded[i] = ratio[j]
			j += 1

	x = np.logical_and(ratio_expanded >= disk_min, ratio_expanded <= disk_max)
	y = np.where(x == True)
	disk_star_indices = ratio_expanded[y]

	if AMR == True and type=="gas":
		disk_star = data["gas","cell_mass"][y]
	else:
		disk_star = data[type,"particle_index"][y]
		
	

	print "NUMBER OF DISK", len(disk_star)

	if plot != None:

		plot_utils.plot_pdf(ratio,(plot + "_jzjcirc_dist.png"), "jz/jcirc","PDF",x_min=-1.0,x_max=1.5,nbins=100,spline=spline)
#		weights = np.ones_like(ratio)/len(ratio)
#		from scipy.interpolate import UnivariateSpline
#		x_vals = bins[:-1] + ((bins[1] - bins[0]) / 2.0 )
#		f = UnivariateSpline(x_vals, hist, s=100)
#		plt.plot(x_vals,f(x_vals),"r-")
#		plt.xlabel('jz/jcirc')
#		plt.ylabel('Distribution')
#		plt.axis([-1.0, 1.5,0.0,(max(hist) + (max(hist) * 0.1))])
#		plt.savefig(plot + "_jzjcirc_dist.png")
#		plt.close()

		# attempt 2
		#from scipy.stats import rv_continuous

		

	#	x_plot = np.linspace(-1.0,1.5,1000)
#		fig = plt.figure()
#		ax = fig.add_subplot(1,1,1)
#		ax.hist(ratio, bins=100, normed=True)
#		ax.plot(x_plot, rv_continuous.pdf(ratio), "r-", label="pdf")
#		ax.legend(loc="best")
#		
#		plt.xlabel('jz/jcirc')
##		plt.ylabel('Distribution')
#		plt.axis([-1.0, 1.5,0.0,0.1])
#		plt.savefig(plot + "_jzjcirc_dist2.png")
#		plt.close()


	
#		hist, bin_edges = np.histogram(ratio,bins=300)
#
#
#		real_bins_temp = np.resize(bin_edges,len(bin_edges)-1)
#		bins = real_bins_temp + 0.5 * np.diff(bin_edges)
#		y = np.zeros(len(hist))
##                inverse_density_function = sp.interpolate.interp1d(
#		for i in range(0,len(y)):
#			y[i] = hist[i].astype(float) / 
#
#		l = plt.plot(bins,y,"r",linewidth=1, marker="o")
#		plt.xlabel('jz/jcirc')
#		plt.ylabel('Distribution')
#		plt.axis([-1.0, 1.5,0.0,(y.max() + 0.02)])
#		plt.savefig(plot + "jzjcirc_dist.png")
#		plt.close()
		
	return disk_star_indices, disk_star, ratio, x

def gradient(x,y,col="r",facecolor="blue",label="Fit Line",filename="name",significance=False,plt=None):
	"""
	plots the gradient ( dy/dx ) of y as a function of x
	if a plt object is supplied, this data is plotted onto the plot object
	"""
	# sort the values
#	print x, y
	x_sorted_temp = np.array(x)
	y_sorted_temp = np.array(y)
	indices = x_sorted_temp.argsort()
	x_sorted = np.sort(x_sorted_temp)
	y_sorted = y_sorted_temp[indices]
#	print "final"
#	print x_sorted, y_sorted

                # predicted values from fitting model
	slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x_sorted,y_sorted)

	if plt != None:
		plt.plot(x_sorted, (slope*x_sorted + intercept), "r", label="Fitted line")

#
# Calculate some additional outputs
#predict_y = intercept + slope * x
#pred_error = y - predict_y
#degrees_of_freedom = len(x) - 2
#residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)

# min and max lin regression
# References: 
# - Statistics in Geography by David Ebdon (ISBN: 978-0631136880)
# - Reliability Engineering Resource Website: 
# - http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
# - University of Glascow, Department of Statistics:
# - http://www.stats.gla.ac.uk/steps/glossary/confidence_intervals.html#conflim
# https://www.students.ncl.ac.uk/tom.holderness/software/pythonlinearfit


#  # Std. deviation of an individual measurement (Bevington, eq. 6.15)  
#  N=numpy.size(x)  
#  sd=1./(N-2.)* numpy.sum((y-a*x-b)**2); sd=numpy.sqrt(sd)  


# more importantly
# http://nbviewer.ipython.org/url/bagrow.com/dsv/LEC10_not_2014-02-13.ipynb
# http://bagrow.com/dsv/
	if significance:
		conf = 0.95
		alpha = 1.0 - conf # significance
		n = x_sorted.size
		linefit = slope*x_sorted + intercept
		x_test = np.linspace(x_sorted.min(),x_sorted.max(),len(x_sorted))

		sd = 1.0 /(n-2.0) * np.sum((y_sorted-slope*x_sorted-intercept)**2) # variance
		sd = np.sqrt(sd) # standard deviation
		sxd = np.sum((x_sorted - x_sorted.mean())**2)
		sx = ((x_test - x_sorted.mean())**2)

	#       t distribution for p=1-alpha/2
		q = sp.stats.t.ppf(1.-alpha/2, n-2)

                # get the upper and lower CI:
		dy = q*sd*np.sqrt( 1./n + sx/sxd )
		yl = linefit - dy
		yu = linefit + dy

                # finally plot


		plt.fill_between(x_sorted, yl, yu, alpha=0.3, facecolor='blue',edgecolor='none')

	if filename:
		plt.savefig(filename)
        return slope, intercept, r_value, p_value, std_err


def sig_uvw(cylinder,galaxy_name,snap,solar=True,disk=True,cold=True):
	"""
	Calculates the sigma_uvw for the galaxy and distributes it by age
	"""

	# gas is really simple
	
	print "###################"
	print "sig_U, sig_V, sig_W"
	if cold == True:
		cold_cylinder = cylinder.cut_region(["obj['temperature'] < 1.5e4"])
		u = cold_cylinder["velocity_cylindrical_radius"]
		v = cold_cylinder["velocity_cylindrical_theta"]
		w = cold_cylinder["velocity_cylindrical_z"]

		sig_u = np.std(u.in_units("km/s").value)
		sig_v = np.std(v.in_units("km/s").value)
		sig_w = np.std(w.in_units("km/s").value)

		print "cold gas sigma"
		print sig_u, sig_v, sig_w

	print "stars"

	# gather all the indices
	filter = np.arange(0,len(cylinder["stars","particle_velocity_cylindrical_radius"]))
	plot_name = galaxy_name + "_uvw"

	if disk == True:
		disk_star_indices, disk_star, ratio, filter = decomp_stars(cylinder,snap,disk_min=0.7, disk_max=1.1,plot=None,r_min=cylinder.ds.arr(2.0,"kpccm"))
		plot_name = plot_name + "_disk"

	if solar == True:

		r_min = cylinder.ds.arr(6.0,"kpc")
		r_max = cylinder.ds.arr(10.0,"kpc")
		z_min = cylinder.ds.arr(-3.0,"kpc")
		z_max = cylinder.ds.arr(3.0,"kpc")

		print filter, len(filter)

		theta_min = cylinder.ds.arr(0.0,"dimensionless")
		theta_max = cylinder.ds.arr(np.radians(359.90),"dimensionless")	

		new_filter = filter_utils.triple_filter_function(cylinder,filter_r="particle_position_cylindrical_radius",filter_z="particle_position_cylindrical_z",
		filter_theta="particle_position_cylindrical_theta",r_min = r_min, r_max = r_max, z_min = z_min, z_max = z_max,
		theta_min = theta_min, theta_max = theta_max, type="stars",extra_filter = filter )
		
		if disk == True:
			filter = new_filter
		else:
			filter = np.where(new_filter == 1)[0]

		print filter, len(filter)
		plot_name = plot_name + "_solar"

	# just consider all the stars
	u = cylinder["stars","particle_velocity_cylindrical_radius"][filter]
	v = cylinder["stars","particle_velocity_cylindrical_theta"][filter]
	w = cylinder["stars","particle_velocity_cylindrical_z"][filter]

	sig_u = np.std(u.in_units("km/s").value)
	sig_v = np.std(v.in_units("km/s").value)
	sig_w = np.std(w.in_units("km/s").value)

	print sig_u, sig_v, sig_w

		# sort out the age distribution
	age = -cylinder["stars","particle_birth_epoch"][filter].in_units("Gyr").value
		# histogram ages

	bin_things = np.linspace(0.0,14.5,30)
	hist, bin_edges = np.histogram(age,bins=bin_things) # gets the distribution
	inds = np.digitize(age, bins=bin_edges) # put the ages into the right bins
	inds = inds - 1.0

	bin_edges = np.delete(bin_edges,len(bin_edges)-1) # deletes the last bin edge elemtn

	# really now we want to shift the bin edges up by the bin width (0.25) for it to be symetric around the bin widths.

	sig_u_age = np.zeros(29)
	sig_v_age = np.zeros(29)
	sig_w_age = np.zeros(29)

	for i in range(0,len(sig_u_age)):
		ages_inds = np.where(inds == i) # find the index locations of stars of this gage
		sig_u_age[i] = np.std(u[ages_inds[0]].in_units("km/s").value)
		sig_v_age[i] = np.std(v[ages_inds[0]].in_units("km/s").value)
		sig_w_age[i] = np.std(w[ages_inds[0]].in_units("km/s").value)



#	temp_age = age
	sig_age_temp = sig_u_age
	age_bins, sig_u_age = plot_utils.kill_first_nan(bin_edges,sig_u_age)
	age_bins, sig_v_age = plot_utils.kill_first_nan(bin_edges,sig_v_age)
	age_bins, sig_w_age = plot_utils.kill_first_nan(bin_edges,sig_w_age)


#	freq_hist, dud = plot_utils.kill_first_nan(hist, sig_age_temp)
	freq_hist = hist

	# errors

	sig_u_age_err = sig_u_age / np.sqrt(freq_hist)
	sig_v_age_err = sig_v_age / np.sqrt(freq_hist)
	sig_w_age_err = sig_w_age / np.sqrt(freq_hist)


	# so inds[id] will return the hist bin for that star.. essentially giving you a filtered distribution
	# and any filter numpy things can filter by value

#	for i in range(0, len(bin_edges)):
	

	# TODO Latex this part

	plot_utils.plot_xy(age_bins,sig_u_age,"Age [Gyr]",r'$\sigma_{U}$ [km/s]',(plot_name + "_u"),type="line",x_min=0.0,x_max=12.5,y_err = sig_u_age_err)
	plot_utils.plot_xy(age_bins,sig_v_age,"Age [Gyr]",r'$\sigma_{V}$ [km/s]',(plot_name + "_v"),type="line", x_min=0.0,x_max=12.5,y_err = sig_v_age_err)
	plot_utils.plot_xy(age_bins,sig_w_age,"Age [Gyr]",r'$\sigma_{W}$ [km/s]',(plot_name + "_w"),type="line", x_min=0.0,x_max=12.5,y_err = sig_w_age_err)
	

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

	plot_utils.plot_xy(age_bins,vel_dispersion,"Age [Gyr]","Velocity Dispersion [km/s]",(plot_name + "_vel_disp"),type="line", x_min=0.0,x_max=12.5,y_err = vel_disp_err)
	
	
	return	

		


def mstar_mhalo(sphere):
	"""
	Returns the result of the mass ratio of the DM to the stellar matter and the corresponding x value
	"""
	dm_mass = sphere["dark","particle_mass"].sum().in_units("Msun")
	stellar_mass = sphere["stars","particle_mass"].sum().in_units("Msun")

	mass_ratio = stellar_mass / dm_mass
	mass_log = np.log10(mass_ratio)

	print mass_ratio, "mass_ratio"
	print mass_log, "mass_log"
	print np.log10(dm_mass.value), "dm mass log"

	x = np.log10(dm_mass.value)
	return x, mass_log


if __name__ == "__main__":
	
	# mainly for testing / an example code
	sim_name = "real_selene_no_delay"
	snapno = None
	redshift = None
	sim = Simulation.load(sim_name)
	sim_patch = "ch"
	if redshift == None:
		if snapno:
			snap = sim.snapshot(snapno, patch=sim_patch)
		else:
			snap = sim.snapshot(sim.num_snapshots(), patch=sim_patch)
	else:
		snap = sim.snapshot(sim.redshift(redshift), patch=sim_patch)
	
	cylinder, disks  = read_utils.load_disk_data(snap, sim_name, sim_patch, halo_id = 0, snapno = None, n_disks=1, disk_h=(5,"kpc"), disk_w=(50,"kpc"), cylinder_w=(50,"kpc"), cylinder_h=(5,"kpc"),extra=None,overwrite=False)

	disk_star_indices, disk_star, ratio, logical = decomp_stars(cylinder,snap,disk_min=0.8, disk_max=1.1)

	# plot the ratio distribution

	# the histogram of the data
	hist, bin_edges = np.histogram(ratio,bins=300)
	real_bins_temp = np.resize(bin_edges,len(bin_edges)-1)
	bins = real_bins_temp + 0.5 * np.diff(bin_edges)
	y = np.zeros(len(hist))
	for i in range(0,len(y)):
		y[i] = hist[i].astype(float) / hist.sum()

	l = plt.plot(bins,y,"r",linewidth=1, marker="o")
	plt.xlabel('jz/jcirc')
	plt.ylabel('Distribution')
	plt.axis([-1.0, 1.5,0.0,(y.max() + 0.02)])
	plt.savefig("jzjcirc.png")
	plt.close()
	
