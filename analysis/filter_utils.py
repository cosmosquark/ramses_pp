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

def triple_filter_function(data,filter_r,filter_z,filter_theta,r_min,r_max,z_min,z_max,theta_min,theta_max,extra_filter=None, type=None):
	"""
	filters a function by two seperate boundaries
	"""
	if type != None:
		if extra_filter == None:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value ) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value ) & \
				(data[type,filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[type,filter_z].in_units(z_max.units).value <= z_max.value) & \
				(data[type,filter_theta].in_units(theta_min.units).value >= theta_min.value) & \
				(data[type,filter_theta].in_units(theta_max.units).value <= theta_max.value)

		else:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value) & \
				(data[type,filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[type,filter_z].in_units(z_max.units).value <= z_max.value) & \
				(data[type,filter_theta].in_units(theta_min.units).value >= theta_min.value) & \
				(data[type,filter_theta].in_units(theta_max.units).value <= theta_max.value) & \
				extra_filter

	else:
		if extra_filter == None:
			filter = (data[filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[filter_r].in_units(r_max.units) <= r_max.value) & \
				(data[filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[filter_z].in_units(z_max.units).value <= z_max.value) & \
				(data[filter_theta].in_units(theta_min.units).value >= theta_min.value) & \
				(data[filter_theta].in_units(theta_max.units).value <= theta_max.value)
		else:
			filter = (data[filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[filter_r].in_units(r_max.units) <= r_max.value) & \
				(data[filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[filter_z].in_units(z_max.units).value <= z_max.value) & \
				(data[filter_theta].in_units(theta_min.units).value >= theta_min.value) & \
				(data[filter_theta].in_units(theta_max.units).value <= theta_max.value) & \
				extra_filter

	return filter





def filter_function(data,filter_r,filter_z,r_min,r_max,z_min,z_max,extra_filter=None, type=None):
	"""
	filters a function by two seperate boundaries
	"""
	print "lol"
	print data[type,filter_r].in_units(r_min.units)
	print r_min, r_max, z_min, z_max
	print extra_filter	
	if type != None:
		if extra_filter == None:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value ) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value ) & \
				(data[type,filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[type,filter_z].in_units(z_max.units).value <= z_max.value)
		else:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value) & \
				(data[type,filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[type,filter_z].in_units(z_max.units).value <= z_max.value) & \
				extra_filter

	else:
		if extra_filter == None:
			filter = (data[filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[filter_r].in_units(r_max.units) <= r_max.value) & \
				(data[filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[filter_z].in_units(z_max.units).value <= z_max.value)
		else:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value) & \
				(data[type,filter_z].in_units(z_min.units).value >= z_min.value) & \
				(data[type,filter_z].in_units(z_max.units).value <= z_max.value) & \
				extra_filter

	return filter

def min_max_filter(data,filter_r,r_min,r_max,extra_filter=None,type=None):
	
	if type != None:
		if extra_filter == None:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value ) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value ) 
		else:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value) & \
				extra_filter

	else:
		if extra_filter == None:
			filter = (data[filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[filter_r].in_units(r_max.units) <= r_max.value) 
		else:
			filter = (data[type,filter_r].in_units(r_min.units).value >= r_min.value) & \
				(data[type,filter_r].in_units(r_max.units) <= r_max.value) & \
				extra_filter

	return filter

def min_max_truths(data,filter_r,r_min,r_max,type=None):
	
	if type != None:
		filter = np.logical_and(data[type,filter_r].in_units(r_min.units).value >= r_min.value, data[type,filter_r].in_units(r_max.units) <= r_max.value )

	else:
		filter = np.logical_and(data[filter_r].in_units(r_min.units).value >= r_min.value, data[filter_r].in_units(r_max.units) <= r_max.value) 

	return filter


