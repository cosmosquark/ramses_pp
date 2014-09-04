from ramses_pp.modules import Simulation
from ramses_pp.modules.pymses import Pymses
from ramses_pp.modules.pymses import PymsesProjection
from pymses.analysis.visualization import *
import matplotlib.pyplot as plt
from pymses.utils import constants as C

import sys, os

verbose = True

def make_movie(fps, filename):
	os.system('ls ./movie/_tmp_*.png | sort -n -t _ -k 3 > filelist_movie')
	os.system("mencoder 'mf://@filelist_movie' -mf type=png:fps=%d -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s.avi"%(fps, filename))
	os.system("mencoder 'mf://@filelist_movie' -mf type=png:fps=%d -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o %s.mp4"%(fps, filename))
	# highest quality
	os.system("mencoder 'mf://@filelist_movie' -mf type=png:fps=%d -ovc x264 -x264encopts preset=slow:tune=film:crf=20 -of rawvideo -o %s.264"%(fps, filename))

	#cleanup
#	os.system("rm ./movie/_tmp_*.png ./movie/_tmp_*.h5")
#	os.system("rm filelist_movie")

def projection(snapshot, source_type, field, method='RayTracer', cmap='jet'):
	info = snapshot.info()
	scale = info["unit_length"].express(C.Mpc)
	#Temp hard-coding
	if field == 'rho':
		#Hydrogens per cubic centimetre
		factor = info["unit_density"].express(C.H_cc)
	else: factor = 1

	#Set our camera
	cam = PymsesProjection.camera(center=[0.5, 0.5, 0.5], region_size=[0.1, 0.1])
	op = PymsesProjection.scalar_operator(field, factor)

	#cam = PymsesProjection.default_camera()

	#Create a PymsesProjection object which points to everything required for projections
	proj = PymsesProjection.load(snapshot, source_type, [field], camera=cam, operator=op)
	if hasattr(proj, method):
		#Call the desired function, if it exists
		func = getattr(proj, method)
		map = func()

		#Save the image
		fname = "movie/_tmp_%d"%ioutput
		#PymsesProjection.save_HDF5(map, fname, PymsesProjection.default_camera(), scale, color_map=cmap)
		h5fname = save_map_HDF5(map, cam, map_name=fname)

		# Save into PNG image file
		fig = save_HDF5_to_plot(h5fname, map_unit=("H/cc",factor), axis_unit=("Mpc", scale), img_path="./",  cmap=cmap)#, cmap_range=crange)
		#plt.savefig('%s_mpl.png'%fname)
		save_HDF5_to_img(h5fname, cmap=cmap, img_path="./")#, cmap_range=crange)
		#mv_str = "mv _tmp_%d movie/"%ioutput
		#os.system(mv_str)

	else:
		raise Exception("Method %s not found"%method)

if __name__ == "__main__":

	if len(sys.argv) == 1:
		print 'Usage: python %s <sim_name> <field> <fps> <outname>'
		sys.exit(1)

	sim_name = sys.argv[1]
	field = sys.argv[2]
	fps = int(sys.argv[3])
	outname = sys.argv[4]
	start = (sys.argv[5])

	if start == None:
		start = 1
	else:
		start = int(start)

	sim = Simulation.load(sim_name)
	num_snapshots = sim.num_snapshots()

	for ioutput in range(start, num_snapshots+1):
		projection(sim.snapshot(ioutput, module='pymses'), Pymses.Type.AMR, field)

	make_movie(fps, outname)
