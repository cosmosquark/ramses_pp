from . import Pymses
from .. import Projection

from pymses.analysis.visualization import *

import matplotlib.pylab as plt
import os, binascii

def default_camera():
	center = [0.5, 0.5, 0.5]
        region_size = [1., 1.]

        cam = Camera(center=center, line_of_sight_axis='x', up_vector='z',
                                region_size=region_size, log_sensitive=True)
	return cam

def default_operator(field):
	func = lambda dset: dset[field]
	op = ScalarOperator(func)
	return op

def load(snapshot, source_type, fields, operator=None, camera=default_camera()):
	if (operator is None):
		#Default to ScalarOperator of first field
		operator = default_operator(fields[0])
	return PymsesProjection(snapshot, source_type, fields, operator, camera)

def save_HDF5(map, mapname, camera, scale=None, units=None, unit_range=None, path="./", color_map="jet"):
	'''
	Save a map to a HDF5 image. Note: scale should be something like ('Mpc', ro.info['unit_length'].express(C.Mpc))
	'''
	h5fname = save_map_HDF5(map, camera, map_name=mapname)

	# Save into PNG image file
	save_HDF5_to_plot(h5fname, map_unit=units, axis_unit=scale, img_path=path, cmap=color_map, cmap_range=unit_range)
	save_HDF5_to_img(h5fname, img_path=path, cmap=color_map, cmap_range=unit_range)

class PymsesProjection(Projection.Projection):
	def __init__(self, snapshot, source_type, fields, operator, camera):
		Projection.Projection.__init__(self, snapshot, camera, "Pymses")
		'''
		Constructor. Set the snapshot in the superclass and also set it's type. For Pymses this should be an AMR or Particle source
		TODO verify we were given a PymsesSnapshot
		'''
		ro = snapshot.raw_snapshot()
		if (source_type == Pymses.Type.AMR):
			self._source = ro.amr_source(fields)
		elif (source_type == Pymses.Type.PART):
			self._source = ro.particle_source(fields)
		else:
			raise Exception("Unrecognised source type: %s"%source_type)
		self._operator = operator
		self._camera = camera
		self._info = ro.info
		self._fields = fields

	def slice(self, depth=0.5):
		'''
		Returns a simple slice of the snapshot. Defaults depth to centre of box
		'''		
		return SliceMap(self._source, self._camera, self._operator, z=depth)

	def FFTProjection(self):
		'''
		Returns an FFT convolved 3D -> 2D projection of the snapshot. Type can be AMR or Particle
		'''
		mp = fft_projection.MapFFTProcessor(self._source, self._info)
		map = mp.process(self._operator, self._camera)
		return map

	def RayTracer(self):
		'''
		Returns a ray traced projection of the snapshot
		'''
		rt = raytracing.RayTracer(self._snapshot.raw_snapshot(), self._fields)
		map = rt.process(self._operator, self._camera)
		return map
