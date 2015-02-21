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
from ramses_pp import config

object_storage = config.yt_object_dir

def load_disk_data(snap, sim_name, sim_patch, halo_id = 0, snapno = None, n_disks=1, disk_h=None, disk_w=None, cylinder_w=None, cylinder_h=None,extra=None,overwrite=False,save=True, center = None, normal = None):
	"""
	Efficient way of loading disks and cylinder
	designed to take either a redshift or a snapshot number
	use this for effective loops
	processes can be speeded up by reloading data if it already exists
	"""

	if disk_h == None:
		disk_h = snap.raw_snapshot().arr(5,"kpccm")

	if disk_w == None:
		disk_w = snap.raw_snapshot().arr(50,"kpccm")

	if cylinder_w == None:
		cylinder_w = snap.raw_snapshot().arr(50,"kpccm")

	if cylinder_h == None:
		cylinder_h = snap.raw_snapshot().arr(5,"kpccm")

	if overwrite == True:
		halos = snap.halos()
		print halos[halo_id]["pos"].in_units("kpccm"), "halo_location"
		disks, cylinder = halos[halo_id].galaxy_disk(n_disks=n_disks,disk_h=disk_h,disk_w=disk_w,cylinder_w=cylinder_w,cylinder_h=cylinder_h, center=center,normal=normal)
		
		if extra == "young": # lets load disks based on the position of the youngest stars
			add_particle_filter("young_stars", function=young_star_filter, filtered_type="all", requires=["particle_age"])
			sphere.ds.add_particle_filter("young_stars")
			disk_length = cylinder_all["young_stars","particle_position_cylindrical_radius"].in_units("kpc").max() 
			#from the cylinder, grab the distance of the furthest young star
			# and then enter some recursion
			return load_disk_data(snap, sim_name, sim_patch, halo_id = halo_id, n_disks=n_disks, disk_h=disk_h, disk_w=disk_length, cylinder_w=disk_length, cylinder_h=cylinder_h,extra=None,overwrite=True,save=save,center=center,normal=normal)

		elif extra == "density":
			# find the radius at which the stellar density is at 10 Msun/pc
			return load_disk_data(snap, sim_name, sim_patch, halo_id = halo_id, n_disks=n_disks, disk_h=disk_h, disk_w=disk_w, cylinder_w=disk_w, cylinder_h=cylinder_h,extra=None,overwrite=True,save=save,center=center,normal=normal)
		else:
			# will fix it for disks later
#			for i in range(0,n_disks):
#				disks[i].save_object((str(halo_id) + "_disk_" + str(i)),(object_storage + "_" + sim_name + "_" + str(snap.output_number) + ".cpkl"  ) )
			print object_storage +  sim_name + "_" + str(snap.output_number()) + ".cpkl"
			if save:
				cylinder.save_object((str(halo_id) + "_cylinder"),(object_storage +  sim_name + "_" + str(snap.output_number()) + ".cpkl" ) )
			return cylinder, disks
	else:
		# might be worth storing the number of disks somewhere safe
		# lets reload these objects rather than spend time recomputing them

		ds = snap.raw_snapshot()
		disks = []
		try:
			data_file = shelve.open((object_storage +  sim_name + "_" + str(snap.output_number()) + ".cpkl"))
		except:
			print "data object does not exist, going to make it"
			return load_disk_data(snap, sim_name, sim_patch, halo_id = halo_id, n_disks=n_disks, disk_h=disk_h, disk_w=disk_w, cylinder_w=cylinder_w, cylinder_h=cylinder_h,extra=extra,overwrite=True,save=save,center=center, normal=normal)
#		for i in range(0,n_disks):
#			ds, disk = data_file[(str(halo_id) + "_disk_" + str(i))]
#			disks.append(disk)
		ds, cylinder = data_file[(str(halo_id) + "_cylinder")]
		return cylinder, disks
