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



def read_yield_table(filename):

	# first find the number of "metallicity" values

	initial = np.genfromtxt(filename, skiprows=0, usecols=(0), comments=None, dtype="string")

	# find the "metallicity" values
	metal_indices = np.char.startswith(initial,"#")

	# filters them
	metal_table = initial[metal_indices]
	# convert metals into floats
	split = np.char.partition(metal_table,"#")

	# we want the numbers... i.e the values AFTER the partition
	# and convert into floats
	metal_table = np.array(split[:,2], dtype="float")
	metal_length = len(metal_table)
	# now extract the size of each of the yield table info for each metal... this is really a 3D array
	
	# now how many rows are there in each element
	yield_length = len(initial) - metal_length
	yield_tab_length = int(yield_length / metal_length)


	# array of dictionarys
	# first value being the metallicity Z
	# rest being a dictionary of values
	yield_table = np.zeros((metal_length,2))
	#data = np.loadtxt(filename, skiprows=1)
	# lets start to populate the yield table

	# this bit is easy
	yield_table[:,0] = metal_table[:]
	yield_table = yield_table.tolist() # so we can actually stick in dictionary data
	# i.e populating the first column with the metal values


	# now for the harder part. since we want to load yield_tab_length rows at a time]
	## well we can simply use # as comments

	data = np.loadtxt(filename,comments="#")


	# seperate data into columns

	#;time, H-2, He, C-4, N, O-6, Ne, Mg-8, Si, Fe-10, M Ejecta-11, MZ-12, Nb of SNII-13, Nb of SNIa-14, mass of SNII, mass of SNIa-16
	# time = Gyr
	# everything else mass of that elements ejected up to that time entry per initial solar mass
	for i in range(0,yield_length,yield_tab_length):
		first = i
		last = (i + yield_tab_length)		
		yield_no = int(np.divide(i,yield_tab_length))

		yield_info = { "time": data[first:last,0],
				"H": data[first:last,1],
				"He": data[first:last,2],
				"C": data[first:last,3],
				"N": data[first:last,4],
				"O": data[first:last,5],
				"Ne": data[first:last,6],
				"Mg": data[first:last,7],
				"Si": data[first:last,8],
				"Fe": data[first:last,9],
				"M_ejecta": data[first:last,10],
				"MZ": data[first:last,11],
				"N_SNII": data[first:last,12],
				"N_SNIa": data[first:last,13],
				"M_SNII": data[first:last,14],
				"M_SNIa": data[first:last,15],
				}
		yield_table[yield_no][1] = yield_info
		yield_table[yield_no][0] = np.float(yield_table[yield_no][0])

	return np.array(yield_table) # much better handled as a numpy array
	

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


def lazy_load_galaxy(sim_name,sim_patch="normal",snap_no=None,halo_id=None,return_snap=True,return_halos=True,filter_stars=True,filter_dark=True,disk_h=10,disk_w=50,disk_units="kpccm"):
	"""
	reads in the simulation and returns the galaxy, halo_sphere the ytsnapshot and the snap object
	since you find yourself doing this a lot, I figured it would make sense as a function in it'w own right
	"""
	sim = Simulation.load(sim_name)

	# super lazy mode activated
	if snap_no == None:
		snap = sim.snapshot(sim.num_snapshots(),module="yt",patch=sim_patch)
	else:
		snap = sim.snapshot(snap_no,module="yt",patch=sim_patch)
	if halo_id == None:
		halo_id = 0


	ytsnap = snap.raw_snapshot()
	all_data = ytsnap.all_data()

	halos = snap.halos()
	print halos
	print halos[0]["pos"], halos[0]["Rvir"]
	galaxy_halo = halos[halo_id]
	print galaxy_halo
	print "getting halo sphere"
	halo_sphere = ytsnap.sphere(galaxy_halo["pos"],galaxy_halo["Rvir"])
	add_particle_filter("stars", function=star_filter, filtered_type="all", requires=["particle_age"])
	halo_sphere.ds.add_particle_filter("stars")
	add_particle_filter("dark", function=dark_filter, filtered_type="all", requires=["particle_age"])
	halo_sphere.ds.add_particle_filter("dark")
	disk_height = ytsnap.arr(disk_h,disk_units)
	disk_width = ytsnap.arr(disk_w,disk_units)
	
	cylinder, disks = load_disk_data(snap, sim_name, sim_patch, halo_id = halo_id, snapno = None, n_disks=1, disk_h=disk_height, disk_w=disk_width, cylinder_w=disk_width, cylinder_h=disk_height,extra=None,overwrite=True,save=False, center = None, normal = None)


	add_particle_filter("stars", function=star_filter, filtered_type="all", requires=["particle_age"])
	cylinder.ds.add_particle_filter("stars")

	add_particle_filter("dark", function=dark_filter, filtered_type="all", requires=["particle_age"])
	cylinder.ds.add_particle_filter("dark")

	return cylinder, halo_sphere, ytsnap, snap
