
#########################################
#
# origional authors go to the yt Development Team. So those guys get all the credit
# copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import h5py

from distutils.version import LooseVersion

from matplotlib.patches import Circle
from matplotlib.colors import colorConverter

from yt.funcs import *
from yt.extern.six import add_metaclass
from yt.visualization._mpl_imports import *
from yt.utilities.physical_constants import \
	sec_per_Gyr, sec_per_Myr, \
	sec_per_kyr, sec_per_year, \
	sec_per_day, sec_per_hr
from yt.units.yt_array import YTQuantity, YTArray
from yt.visualization.image_writer import apply_colormap
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from yt.visualization import _MPL

callback_registry = {}




class RegisteredCallback(type):
	def __init__(cls, name, b, d):
		type.__init__(cls, name, b, d)
		callback_registry[name] = cls


@add_metaclass(RegisteredCallback)
class PlotCallback(object):
	def __init__(self, *args, **kwargs):
		pass

	def convert_to_plot(self, plot, coord, offset = True):
		# coord should be a 2 x ncoord array-like datatype.
		try:
			ncoord = np.array(coord).shape[1]
		except IndexError:
			ncoord = 1

		# Convert the data and plot limits to tiled numpy arrays so that
		# convert_to_plot is automatically vectorized.

		x0 = np.array(np.tile(plot.xlim[0],ncoord))
		x1 = np.array(np.tile(plot.xlim[1],ncoord))
		xx0 = np.tile(plot._axes.get_xlim()[0],ncoord)
		xx1 = np.tile(plot._axes.get_xlim()[1],ncoord)

		y0 = np.array(np.tile(plot.ylim[0],ncoord))
		y1 = np.array(np.tile(plot.ylim[1],ncoord))
		yy0 = np.tile(plot._axes.get_ylim()[0],ncoord)
		yy1 = np.tile(plot._axes.get_ylim()[1],ncoord)

		ccoord = np.array(coord)
		
		# We need a special case for when we are only given one coordinate.
		if ccoord.shape == (2,):
			return ((ccoord[0]-x0)/(x1-x0)*(xx1-xx0) + xx0,
			 (ccoord[1]-y0)/(y1-y0)*(yy1-yy0) + yy0)
		else:
			return ((ccoord[0][:]-x0)/(x1-x0)*(xx1-xx0) + xx0,
			 (ccoord[1][:]-y0)/(y1-y0)*(yy1-yy0) + yy0)


	def pixel_scale(self,plot):
		x0, x1 = np.array(plot.xlim)
		xx0, xx1 = plot._axes.get_xlim()
		dx = (xx1 - xx0)/(x1 - x0)

		y0, y1 = np.array(plot.ylim)
		yy0, yy1 = plot._axes.get_ylim()
		dy = (yy1 - yy0)/(y1 - y0)

		return (dx,dy)



class HaloCatalogCallback(PlotCallback):	
	"""
	annotate_halos(halo_catalog, circle_kwargs=None,
		width = None, annotate_field=False,
        font_kwargs = None, factor = 1.0)

	Plots circles at the locations of all the halos
	in a halo catalog with radii corresponding to the
	virial radius of each halo. 

	circle_kwargs: Contains the arguments controlling the
		appearance of the circles, supplied to the 
		Matplotlib patch Circle.

	width: the width over which to select halos to plot,
      		useful when overplotting to a slice plot. Accepts
       		a tuple in the form (1.0, 'Mpc').

	annotate_field: Accepts a field contained in the 
		halo catalog to add text to the plot near the halo.
		Example: annotate_field = 'particle_mass' will
		write the halo mass next to each halo.

	font_kwargs: Contains the arguments controlling the text
		appearance of the annotated field.
    
	factor: A number the virial radius is multiplied by for
		plotting the circles. Ex: factor = 2.0 will plot
		circles with twice the radius of each halo virial radius.
	"""	

	_type_name = 'halos'
	region = None
	_descriptor = None

	def __init__(self, halo_catalog, plot, circle_kwargs = None,
		width = None, annotate_field = False,
		font_kwargs = None, factor = 1.0):
		print "haha"
		print dir(plot.data_source.axis)
		print dir(plot.ds.coordinates)
		PlotCallback.__init__(self)
		self.halo_catalog = halo_catalog
		self.width = width
		self.annotate_field = annotate_field
		self.font_kwargs = font_kwargs
		self.factor = factor
		if circle_kwargs is None:
			circle_kwargs = {'edgecolor':'white', 'facecolor':'None'}
		self.circle_kwargs = circle_kwargs
		print "loading"

		print plot.xlim
		data = plot.data_source
 		x0, x1 = plot.xlim
		y0, y1 = plot.ylim
#		xx0, xx1 = plot._axes.get_xlim()
#		yy0, yy1 = plot._axes.get_ylim()
		xx0, xx1 = plot.xlim ## not right, translate for above
		yy0, yy1 = plot.ylim ### not right

		halo_data= self.halo_catalog
		axis_names = plot.ds.coordinates.axis_name
		xax = plot.ds.coordinates.x_axis[data.axis]
		yax = plot.ds.coordinates.y_axis[data.axis]
		field_x = "particle_position_%s" % axis_names[xax]
		field_y = "particle_position_%s" % axis_names[yax]
		field_z = "particle_position_%s" % axis_names[data.axis]
		print field_x,field_y,field_z,"fields x y z"
		plot._axes.hold(True)
		print "calling"
  	   	   # Set up scales for pixel size and original data
		pixel_scale = self.pixel_scale(plot)[0]
		data_scale = data.ds.length_unit
		units = data_scale.units

		# Convert halo positions to code units of the plotted data
		# and then to units of the plotted window
		px = halo_data._nhalo[field_x][:].in_units(units) / data_scale  
		py = halo_data._nhalo[field_y][:].in_units(units) / data_scale
		px, py = self.convert_to_plot(plot,[px,py])

		# Convert halo radii to a radius in pixels
		radius = halo_data._nhalo[:]['Rvir'].in_units(units)
		radius = np.array(radius*pixel_scale*self.factor/data_scale)

		if self.width:
			pz = halo_data._nhalo[field_z][:].in_units(units)/data_scale
			pz = data.ds.arr(pz, 'code_length')
			c = data.center[data.axis]

		# I should catch an error here if width isn't in this form
 		# but I dont really want to reimplement get_sanitized_width...
			width = data.ds.arr(self.width[0], self.width[1]).in_units('code_length')

			indices = np.where((pz > c-width) & (pz < c+width))

			px = px[indices]
			py = py[indices]
			radius = radius[indices]

		for x,y,r in zip(px, py, radius):
			plot._axes.add_artist(Circle(xy=(x,y),
				radius = r, **self.circle_kwargs))

		plot._axes.set_xlim(xx0,xx1)
		plot._axes.set_ylim(yy0,yy1)
		plot._axes.hold(False)

		if self.annotate_field:
			annotate_dat = halo_data[self.annotate_field]
			texts = ['{0}'.format(dat) for dat in annotate_dat]
			for pos_x, pos_y, t in zip(px, py, texts):
				plot._axes.text(pos_x, pos_y, t, **self.font_kwargs)
                                                   


