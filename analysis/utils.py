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
from ramses_pp.analysis import read_utils


def radial_indices(field,data,correction=YTArray(1,"pc")):
	"""
	use this function to create radial bins of particles based on the smallest possible cell size
	and bin increments of the smallest cell size
	of a radial field of your choice

	field = e.g data["stars","particle_spherical_position_radius"]
	data = a data object so you can do data.ds
	correction = a slight correction so that sphere too small errors do not appear
	"""
	
	# computes how many bins will be needed depending on the size of the object
	r_indices = np.floor(data[field].in_units("code_length").value / data.ds.index.get_smallest_dx().value).astype(int) 
	r_bins = YTArray(np.zeros(int(r_indices.max() + 1)),"cm")
	for i in range(0,len(r_bins)):
		r_bins[i] = (i + 1) * data.ds.index.get_smallest_dx().in_units("cm")
	r_bins[0] += correction
	print r_indices, "the first"
	return r_indices, r_bins


def mass_enclosed_bins(object,snap,r_bins,shape="sphere",type="all"):
	"""
	extending radial_indices, this function can find the enclosed mass of
	of a particular field (stars, dm, gas, total or custom)
	returns the m_bins to accompany the r_bins
	"""

	m_bins = YTArray(np.zeros(len(r_bins)),"g")
	invalid = False
	for i in range(0,len(r_bins)):
		if type == "all":
			m_bins[i] = snap.raw_snapshot().sphere(cylinder.center,r_bins[i]).quantities.total_quantity(["particle_mass"]) + snap.raw_snapshot().sphere(cylinder.center,r_bins[i]).quantities.total_quantity(["cell_mass"])
		elif type == "gas":
			m_bins[i] = snap.raw_snapshot().sphere(cylinder.center,r_bins[i]).quantities.total_quantity(["cell_mass"])
		elif type == "stars":
			 m_bins[i] = snap.raw_snapshot().sphere(cylinder.center,r_bins[i]).quantities.total_quantity(["stars","particle_mass"])
		elif type == "dark":
			 m_bins[i] = snap.raw_snapshot().sphere(cylinder.center,r_bins[i]).quantities.total_quantity(["dark","particle_mass"])
		else:
			print "invalid type"
			invalid = True
			break
	if invalid == True:
		return None
	else:
		return m_bins


def vcirc(r_bins,m_bins):
	"""
	computes vcirc for the enclosed material
	"""

	vcirc = YTArray(np.zeros(len(m_bins)))
	vcirc = (np.sqrt(G * m_bins / r_bins))
#	for i in range(0,len(m_bins)):
#		vcirc[i] = (np.sqrt(G * m_bins[i] / r_bins[i]))
	print vcirc
	return vcirc

def jcirc_bins(r_bins,v_circ,m_bins):
	""" this computes j_circ for each individual r_bin, v_circ bin and m_bin"""
	j_circ = YTArray(np.zeros(len(m_bins)))
	j_circ = (r_bins * v_circ)
	return j_circ

def jcirc(data, r_indices,v_circ):
	""" this computes j_circ for each individual r_bin, v_circ bin and m_bin"""
	star_count = len(data['particle_age'] != 0)
	print r_indices, "r_indices"
	print v_circ, "v_circ"
	j_circ = YTArray(np.zeros(star_count))
	j_circ = (data["stars", "particle_position_spherical_radius"] * v_circ[r_indices.astype(int)])
	#j_circ = data["stars", "particle_position_spherical_radius"]) * np.sqrt(G * m_bins[r_indices.astype(int)] / data["stars", "particle_position_spherical_radius"])
	print j_circ, "j_circ"
	return j_circ



def vrot(data):
	"""
	returns vrot (it's the spherical theta componant
	"""
	vrot = data["stars","particle_velocity_spherical_theta"]
	rrot = data["stars","particle_position_spherical_radius"]
	return vrot

def jz(data):
	"""
	returns jz (its the specific angular momentum in the z axis
	assuming your object has a centre and a normal vector
	"""	
	jz = data["stars","particle_specific_angular_momentum_z"]
	return jz

def decomp_stars(data,snap,disk_min=0.8, disk_max=1.1):
	"""
	At the moment, this technique can only filter for disk stars
	This is computed by calculating the ratio of Jz/Jcirc
	See this paper http://arxiv.org/pdf/1210.2578.pdf
	returns the indices, id's and the ratio
	"""

	
	add_particle_filter("stars", function=star_filter, filtered_type="all", requires=["particle_age"])
        data.ds.add_particle_filter("stars")
	j_z = jz(data)
	r_indices, r_bins = radial_indices(("stars","particle_spherical_position_radius"),data)
	print r_indices, "the second"
	m_bins = mass_enclosed_bins(data,snap,r_bins,shape="sphere",type="all")
	print r_indices, "the third"
	v_circ = vcirc(r_bins,m_bins)
	print r_indices, "the fourth"
	j_circ = jcirc(data,r_indices,v_circ)

	print j_z
	print j_circ
	ratio = j_z / j_circ

	print ratio
	# get indices of stars which satisfy the conditions

	x = np.logical_and(ratio >= disk_min, ratio <= disk_max)
	y = np.where(x == True)
	print x
	print y
	disk_star_indices = ratio[y]
	disk_star = data["stars","particle_index"][y]

	return disk_star_indices, disk_star, ratio

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

	disk_star_indices, disk_star, ratio = decomp_stars(cylinder,snap,disk_min=0.8, disk_max=1.1)

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
	
