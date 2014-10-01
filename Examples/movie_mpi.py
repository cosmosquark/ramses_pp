from mpi4py import MPI

#from ramses_pp.modules import Simulation
from ramses_pp.modules.pymses import Pymses
from ramses_pp.modules.pymses import PymsesProjection
from pymses.analysis.visualization import *
import matplotlib.pyplot as plt
from pymses.utils import constants as C

import sys, os

verbose = True
ncpu = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
comm = MPI.COMM_WORLD

def terminate(rank, ierr):
	print 'CPU %d called terminate with error-code %d. Shutting down'%(rank, ierr)
	MPI.COMM_WORLD.Abort(ierr)

def make_movie(fps, filename):
	print 'CPU on rank %d processing movie...'%rank
	os.system('ls _tmp_*.png | sort -n -t _ -k 3 > filelist')
	os.system("mencoder 'mf://@filelist' -mf type=png:fps=%d -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s.avi"%(fps, filename))

	#cleanup
	os.system("rm ./_tmp_*.png ./_tmp_*.h5")
	os.system("rm filelist")

def projection(sim_dir, ioutput, source_type, field, method='Slice', cmap='gist_stern_r'):
	if verbose:
		print 'CPU on rank %d taking projection of output %05d'%(rank, ioutput)
	#Load the snapshot
	snapshot = Pymses.load(sim_dir, ioutput)
	z = snapshot.current_redshift()
	aexp = 1/(1+z)

	info = snapshot.info()
	#scale = info["unit_length"].express(C.Mpc)
	scale = info["unit_length"].express(C.Mpc) / aexp
	#Temp hard-coding
	if field == 'rho':
		factor = info["unit_density"].express(C.H_cc)
	else: factor = None

	#Set our camera
	cam = PymsesProjection.camera(center=[0.5, 0.45, 0.58], region_size=[0.05, 0.05], log_sensitive=True)
	op = PymsesProjection.scalar_operator(field, factor)
	#cam = PymsesProjection.default_camera()

	#Create a PymsesProjection object which points to everything required for projections
	proj = PymsesProjection.load(snapshot, source_type, [field], camera=cam, operator=op)
	if hasattr(proj, method):
		#Call the desired function, if it exists
		func = getattr(proj, method)
		map = func()

		#Save the image
		fname = "_tmp_%d"%ioutput
		#PymsesProjection.save_HDF5(map, fname, PymsesProjection.default_camera(), scale, color_map=cmap)
		print 'CPU on rank %d saving file %s'%(rank, fname)
		h5fname = save_map_HDF5(map, cam, map_name=fname)

		# Save into PNG image file
		plt.contour(map)
		fig = save_HDF5_to_plot(h5fname, map_unit=("H/cc",factor), axis_unit=(r"Mpc  (1 + z)$^{-1}$", scale),
			 img_path="./", cmap=cmap, save_into_png=True)#, cmap_range=crange)
		
		#plt.contour(map)
		#plt.savefig('%s.png'%fname)
		#plt.clf()
		#save_HDF5_to_img(h5fname, img_path='./', cmap=cmap, adaptive_gaussian_blur=True, ramses_output=snapshot.raw_snapshot())#, cmap_range=crange)

	else:
		terminate(rank, 404)

if __name__ == "__main__":
	host = (rank == 0)

	if host:
		if len(sys.argv) == 1:
			print 'Usage: python %s <sim_dir> <i_min> <i_max> <field> <fps> <outname>'
			terminate(rank, 0)
		if os.path.isdir(sys.argv[1]):
			data = {'sim_dir':sys.argv[1],
				'i_min':int(sys.argv[2]),
				'i_max':int(sys.argv[3]),
				'field':sys.argv[4],
				'fps':int(sys.argv[5]),
				'outname':sys.argv[6]}
		else:
			terminate(rank, 404)
	else:
		data = None
	data = comm.bcast(data, root=0)
	if verbose and not host:
		print 'CPU on rank %d received data: '%rank, data

	irange = range(data['i_min'], data['i_max']+1)
	nfiles = len(irange)
	files_per_cpu = nfiles/ncpu
	if host and verbose:
		print 'ncpu: %d, nfiles: %d'%(ncpu, nfiles)
		print 'Files per cpu: %d'%files_per_cpu

	if ncpu == nfiles:
		print 'One MPI process per file...'
		ioutput = irange[rank]
		projection(data['sim_dir'], ioutput, Pymses.Type.AMR, data['field'])
	else:
		print 'More files than MPI processes...splitting'
			#Share out the files
		firstfile = (rank)*files_per_cpu
		lastfile = min(firstfile + files_per_cpu, nfiles)
		for i in range(firstfile, lastfile+1):
			ioutput = irange[i]
			projection(data['sim_dir'], ioutput, Pymses.Type.AMR, data['field'])

	#Sync
	MPI.COMM_WORLD.Barrier()
	if host:
		make_movie(data['fps'], data['outname'])
