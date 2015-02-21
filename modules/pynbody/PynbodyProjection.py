#from . import Pynbody
from .. import Projection

import pynbody
from . import Pynbody
from ..pymses import PymsesProjection

def camera(center=[0.5, 0.5, 0.5], region_size=[1., 1.], los_axis='x',
		up_vector='z', log_sensitive=True, **kwargs):
	return PymsesProjection.camera(center=center, los_axis=los_axis, up_vector=up_vector,
                                region_size=region_size, log_sensitive=log_sensitive, **kwargs)

def load(snapshot, camera=camera()):
	return PynbodyProjection(snapshot, camera)

class PynbodyProjection(Projection.Projection):
	def __init__(self, snapshot, camera):
		'''
		Constructor. Set the snapshot and (pymses) camera in the superclass
		'''
		Projection.Projection.__init__(self, snapshot, camera, "Pynbody")
		if (not self._snapshot.has_key('centred') or self._snapshot.get('centered') is False):
			#Shift the domain to centre the image on the box
			s = self._snapshot.raw_snapshot()
			#Only shift if we haven't already
			s['pos'] -= 0.5
			self._snapshot.put('centered', True)

	def image(self, field, units, filt_func=None, width=None, axis_units='Mpc a h^-1', ret_im=True, **kwargs):
		'''
		Shortcut to pynbody.sph.image, which centres the domain
		'''
		s = self._snapshot.raw_snapshot()
		if width is None:
			#Get it from the camera
			from pynbody.array import SimArray
			boxsize = s.properties['boxsize']
			print boxsize
			arr = SimArray(self._camera.region_size, boxsize)
			arr.conversion_context = s.conversion_context
			arr = arr.in_units(axis_units)
			#width = float(arr[0])/2.
			width = float(arr[0])

		if not filt_func is None:
			s = filt_func(s)
		#s.physical_units()
		#s.original_units()
		#print s['pos']
		#s['pos'].convert_units(axis_units)
		s[field].in_units(units)
		print s[field]
		#im = pynbody.plot.sph.image(s, width='%f %s'%(width, axis_units), units=units, ret_im=ret_im, **kwargs)
		im = pynbody.plot.sph.image(s, width='%f %s'%(width, axis_units), ret_im=ret_im, **kwargs)
		return im

	def velocity_image(self, field, units, filt_func=None, mode='stream', width=None, axis_units='Mpc a h^-1', cmap='gist_stern_r', ret_im=True, **kwargs):
		'''
		Shortcut to pynbody.sph.image, which centres the domain
		'''
		s = self._snapshot.raw_snapshot()
		if width is None:
			#Get it from the camera
			from pynbody.array import SimArray
			boxsize = s.properties['boxsize']
			print boxsize
			arr = SimArray(self._camera.region_size, boxsize)
			arr.conversion_context = s.conversion_context
			arr = arr.in_units(axis_units)
			#width = float(arr[0])/2.
			width = float(arr[0])

		if not filt_func is None:
			s = filt_func(s)
		#s.physical_units()
		#s.original_units()
		#print s['pos']
		#s['pos'].convert_units(axis_units)
		s[field].in_units(units)
		print s[field]
		#im = pynbody.plot.sph.image(s, width='%f %s'%(width, axis_units), units=units, ret_im=ret_im, **kwargs)
		im = pynbody.plot.sph.velocity_image(s, mode=mode, width='%f %s'%(width, axis_units), ret_im=ret_im, **kwargs)
		return im
