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
from yt.utilities.math_utils import ortho_find
import re
import yt.visualization.eps_writer as eps
from yt.extern.six import iteritems

from yt.funcs import get_image_suffix
import os


width_thing = YTArray(100,"kpc")
depth_thing = YTArray(10,"kpc")

from ramses_pp.analysis import plot_utils

def plot_frb_profile(image,width,y_units,x_lab,y_lab,filename,n_bins=50,ylog=True):
	return plot_utils.plot_frb_profile(image=image,width=width,y_units=y_units,x_lab=x_lab,y_lab=y_lab,filename=filename,n_bins=n_bins,ylog=ylog)


def plot_eps_image(plot, name, multi=False, mpl_kwargs = None, suffix=None):

	print name, "thign"
	eps_fig = eps.single_plot(plot)
	n = name[0].split(".")[0]
	n = n + ".pdf"
	eps_fig.save_fig(name, format='pdf')
	

def gen_data_source(axis,container,ytsnap,width,depth,axis_unit):

	if axis==0:
		try:
			x = container.center[0].in_units("code_length").value - ytsnap.arr(depth[0],axis_unit).in_units("code_length").value
			y = container.center[1].in_units("code_length").value - ytsnap.arr(width[0][0],axis_unit).in_units("code_length").value
			z = container.center[2].in_units("code_length").value - ytsnap.arr(width[1][0],axis_unit).in_units("code_length").value
			left = [x,y,z]

			x = container.center[0].in_units("code_length").value + ytsnap.arr(depth[0],axis_unit).in_units("code_length").value
			y = container.center[1].in_units("code_length").value + ytsnap.arr(width[0][0],axis_unit).in_units("code_length").value
			z = container.center[2].in_units("code_length").value + ytsnap.arr(width[1][0],axis_unit).in_units("code_length").value
			right = [x,y,z]
		except:
			x = container.center[0].in_units(axis_unit).value - ytsnap.arr(depth[0],axis_unit).in_units(axis_unit).value
			y = container.center[1].in_units(axis_unit).value - ytsnap.arr(width[0][0],axis_unit).in_units(axis_unit).value
			z = container.center[2].in_units(axis_unit).value - ytsnap.arr(width[1][0],axis_unit).in_units(axis_unit).value
			left = [x,y,z]

			x = container.center[0].in_units(axis_unit).value + ytsnap.arr(depth[0],axis_unit).in_units(axis_unit).value
			y = container.center[1].in_units(axis_unit).value + ytsnap.arr(width[0][0],axis_unit).in_units(axis_unit).value
			z = container.center[2].in_units(axis_unit).value + ytsnap.arr(width[1][0],axis_unit).in_units(axis_unit).value
			right = [x,y,z]			

	if axis==1:
	
		try:
			x = container.center[0].in_units("code_length").value - ytsnap.arr(width[1][0],axis_unit).in_units("code_length").value
			y = container.center[1].in_units("code_length").value - ytsnap.arr(depth[0],axis_unit).in_units("code_length").value
			z = container.center[2].in_units("code_length").value - ytsnap.arr(width[0][0],axis_unit).in_units("code_length").value
			left = [x,y,z]

			x = container.center[0].in_units("code_length").value + ytsnap.arr(width[1][0],axis_unit).in_units("code_length").value
			y = container.center[1].in_units("code_length").value + ytsnap.arr(depth[0],axis_unit).in_units("code_length").value
			z = container.center[2].in_units("code_length").value + ytsnap.arr(width[0][0],axis_unit).in_units("code_length").value
			right = [x,y,z]
		except:
			x = container.center[0].in_units(axis_unit).value - ytsnap.arr(width[1][0],axis_unit).in_units(axis_unit).value
			y = container.center[1].in_units(axis_unit).value - ytsnap.arr(depth[0],axis_unit).in_units(axis_unit).value
			z = container.center[2].in_units(axis_unit).value - ytsnap.arr(width[0][0],axis_unit).in_units(axis_unit).value
			left = [x,y,z]

			x = container.center[0].in_units(axis_unit).value + ytsnap.arr(width[1][0],axis_unit).in_units(axis_unit).value
			y = container.center[1].in_units(axis_unit).value + ytsnap.arr(depth[0],axis_unit).in_units(axis_unit).value
			z = container.center[2].in_units(axis_unit).value + ytsnap.arr(width[0][0],axis_unit).in_units(axis_unit).value
			right = [x,y,z]
	if axis==2:

		try:
			x = container.center[0].in_units("code_length").value - ytsnap.arr(width[0][0],axis_unit).in_units("code_length").value
			y = container.center[1].in_units("code_length").value - ytsnap.arr(width[1][0],axis_unit).in_units("code_length").value
			z = container.center[2].in_units("code_length").value - ytsnap.arr(depth[0],axis_unit).in_units("code_length").value
			left = [x,y,z]

			x = container.center[0].in_units("code_length").value + ytsnap.arr(width[0][0],axis_unit).in_units("code_length").value
			y = container.center[1].in_units("code_length").value + ytsnap.arr(width[1][0],axis_unit).in_units("code_length").value
			z = container.center[2].in_units("code_length").value + ytsnap.arr(depth[0],axis_unit).in_units("code_length").value
			right = [x,y,z]
		except:
			x = container.center[0].in_units(axis_unit).value - ytsnap.arr(width[0][0],axis_unit).in_units(axis_unit).value
			y = container.center[1].in_units(axis_unit).value - ytsnap.arr(width[1][0],axis_unit).in_units(axis_unit).value
			z = container.center[2].in_units(axis_unit).value - ytsnap.arr(depth[0],axis_unit).in_units(axis_unit).value
			left = [x,y,z]

			x = container.center[0].in_units(axis_unit).value + ytsnap.arr(width[0][0],axis_unit).in_units(axis_unit).value
			y = container.center[1].in_units(axis_unit).value + ytsnap.arr(width[1][0],axis_unit).in_units(axis_unit).value
			z = container.center[2].in_units(axis_unit).value + ytsnap.arr(depth[0],axis_unit).in_units(axis_unit).value
			right = [x,y,z]

	data_source = ytsnap.region(container.center, left, right)
	return data_source



def plot_viz(plot,name,suffix):
	plot.save(name)
	if suffix != None:
		for suf in suffix:
			plot.save(name, suffix = suf)
		


def visualisation(viz_type, container, raw_snapshot, module=config.default_module, gas=True, stars=False, dark=False, gas_fields=["temperature","density"],gas_units=["K","g/cm**3"], gas_extrema=[[None,None],[None,None]], dark_fields=[('deposit', 'dark_density')], dark_units=["g/cm**3"], dark_extrema=[[None,None]], star_fields=[('deposit', 'stars_density')], star_units=["g/cm**3"], star_extrema=[[None,None]], filter=None, return_objects=False, return_fields=False, callbacks=[], width=width_thing,extra_width=None,depth=depth_thing, name="plot", format=".png", prefix="", normal_vector=[1.0,0.0,0.0], axis=[0,1,2], plot_images = True, weight_field = None, image_width = 1000, extra_image_width=None, suffix=None):

	"""
	This routine is designed to handle almost all of the visualisation routines that you can think of.

	this has a disgusting amount of kwargs, so lets break everything down

	arguments

	viz_type = slice, projection, off_axis_slice, off_axis_projection. Not available for all modules all these combinations. But it will know what you mean depending on the config that you load
	container = the data container or object of interest. This will accept either a YT raw snapshot or a pymses raw snapshot
	module = yt or pymses

	kwargs

	raw_snapshot = the raw snapshot of the simulation (i.e the whole box or a larger portion of the box). This may be used to make an additional cut for the datasource
	gas = do we want to plot gas data, default to True
	stars = do we want to plot stellar data, default to false
	dark = do we want to plot dark matter data, default to false

	gas_fields = a list of strings (or turples for some YT magic) of fields that we want to plot
	gas_units = the units of those fields

	dark_fields = see above
	dark_units = see above

	star_fields = see above
	star_units = see above

	filter = if there needs to be any additonal filters used... but really that should be handled before this, but this is more specific for YT
	return_objects = whether you wish to return a list of plot objects and images or not for further work
	callbacks = extra callbacks? this is YT specific
	
	width = YTArray containing the width of the image and its unit
	extra_width = YTArray if you want to do something a bit more fancy with the image dimensions

	depth = projection depth if applicable
	name = name of the plot if you want to rename it
	format = the output format
	
	normal_vector = the normal vector for the off axis plots
	axis = axis for on axis plots
	data_source = YT data source if applicable

	plot_images = if you actually want to plot something

	TODO use regex to rename dark_deposit with io_deposit if age does not exist
	"""

	# customise the width etc

#	if not [("all","particle_age")] in container:
#		print "WARNING, DM only run for particles. Viz may fail since the " + \
#			"age field does not exist in DM only runs"

	axis_unit = width.units
	if extra_width != None:
		width = ((width.v,width.units),(extra_width.v,width.units))
	else:
		width = ((width.v,width.units),(width.v,width.units))

	if extra_image_width != None:
		image_width = (image_width,extra_image_width)
	else:
		image_width = (image_width,image_width) 

	depth = (depth.v,depth.units)

	fields = []
	units = []
	extrema = []

	if gas_fields and gas == True:
		fields = fields + gas_fields
		units = units + gas_units
		extrema = extrema + gas_extrema	

	if star_fields and stars == True:
		fields = fields + star_fields
		units = units + star_units
		extrema = extrema + star_extrema

	if dark_fields and dark == True:
		fields = fields + dark_fields
		units = units + dark_units
		extrema = extrema + dark_extrema

	# correct length of the extrema
	# usually if this is the case, then we obviously don't care

	if len(fields) > len(extrema):
		diff = len(fields) - len(extrema)
		for i in range(0,diff):
			extrema = extrema + [[None, None]]

	for i in range(0, len(extrema)):
		if extrema[i][0] == None:
			extrema[i][0] = "min"
		if extrema[i][1] == None:
			extrema[i][1] = "max"


	print fields
	plots = {}
	frb = {}

	basis_vectors = ortho_find(normal_vector)
	north_vectors = ortho_find(normal_vector)

	print "basis_vectors", basis_vectors

	# check if basis vectors are in right units

	basis_vectors = YTArray(basis_vectors,"g*kpc**2/s")
	north_vectors = YTArray(north_vectors,"g*kpc**2/s")

	print basis_vectors
	# select the module
	if module == "yt":

		# select the type of plot
		if viz_type == "slice":
			# axis
			print "on axis slice plot"
			if 0 in axis:
				# x

				# to note, a slice plot is pretty irrelevent of the data source.
				# see e.g https://bpaste.net/show/f08dda811ae0
				# which will reproduce the same result
				print "plotting axis 0"		
				print width[0][0], "width"
				plot = yt.SlicePlot(container.ds,0,fields,center=container.center,width=width)
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), image_width[0])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				plots["0_plot"] = plot
				frb["0_frb"] = image

				if plot_images:
					plot_viz(plot,name,suffix)

			# y
			if 1 in axis:
				# x
				print "plotting axis 1"		
				plot = yt.SlicePlot(container.ds,1,fields,center=container.center,width=width)
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				plots["1_plot"] = plot
				frb["1_frb"] = image

				if plot_images:
					plot_viz(plot,name,suffix)
			# z

			if 2 in axis:
				# x
				print "plotting axis 2"
				plot = yt.SlicePlot(container.ds,2,fields,center=container.center,width=width)
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				plots["2_plot"] = plot
				frb["2_frb"] = image

				if plot_images:
					plot_viz(plot,name,suffix)

		if viz_type == "projection":

			# axis
			print "on axis projection plot"
			if 0 in axis:
				# x
				
				# in this instance.. the depth is along the axis that you are viewing "down".. so in a way is user defined
				data_source = gen_data_source(0,container,raw_snapshot,width,depth,axis_unit)
				print "plotting axis 0"		
				plot = yt.ProjectionPlot(container.ds,0,fields,center=container.center,width=width,weight_field=weight_field, data_source=data_source)
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				plots["0_plot"] = plot
				frb["0_frb"] = image

				if plot_images:
					plot_viz(plot,name,suffix)
			
				# y
			if 1 in axis:
				print "plotting axis 1"
				data_source = gen_data_source(1,container,raw_snapshot,width,depth,axis_unit)
				plot = yt.ProjectionPlot(container.ds,1,fields,center=container.center,width=width,weight_field=weight_field, data_source=data_source)
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				plots["1_plot"] = plot
				frb["1_frb"] = image

				if plot_images:
					plot_viz(plot,name,suffix)

				# z	
			if 2 in axis:
				print "plotting axis 2"
				data_source = gen_data_source(2,container,raw_snapshot,width,depth,axis_unit)
				plot = yt.ProjectionPlot(container.ds,2,fields,center=container.center,width=width,weight_field=weight_field, data_source=data_source)
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])
				
				plots["2_plot"] = plot
				frb["2_frb"] = image

				if plot_images:
					plot_viz(plot,name,suffix)


		if viz_type == "off_axis_slice":

				# off axis.. may as well do all 3
				# same rules RE the slice apply here too
			print "off axis slice plot"
			if 0 in axis:
				print "plotting axis 0"
				plot = yt.OffAxisSlicePlot(container.ds,basis_vectors[0],fields,center=container.center,width=width, north_vector=north_vectors[0])
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])
				
				plots["0_plot"] = plot
				frb["0_frb"] = image

				if plot_images:
					plot_viz(plot,name + "_axis_0",suffix)

			if 1 in axis:
				print "plotting axis 1"	
				plot = yt.OffAxisSlicePlot(container.ds,basis_vectors[1],fields,center=container.center,width=width, north_vector=north_vectors[0])
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])
				
				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				plots["1_plot"] = plot
				frb["1_frb"] = image

				if plot_images:
					plot_viz(plot,name + "_axis_1",suffix)


			if 2 in axis:
				print "plotting axis 2"
				plot = yt.OffAxisSlicePlot(container.ds,basis_vectors[2],fields,center=container.center,width=width, north_vector=north_vectors[0])
				image = plot.data_source.to_frb(yt.YTQuantity(width[0][0], axis_unit), [image_width[0], image_width[1]])
				
				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				plots["2_plot"] = plot
				frb["2_frb"] = image

				if plot_images:
					plot_viz(plot,name + "_axis_2",suffix)
		

		if viz_type == "off_axis_projection":

			print "off axis projection plot"


			if 0 in axis:	
				print "plotting axis 0"
				# since there is no way of doing an off axis box.. our data source is essentially a cube.. it ignore depth
				# oh, and you probably want this rather than your disk dataset

				#TODO make sure filtered datasts carry across... but it seems to be the case here.
				data_source = gen_data_source(0,container,raw_snapshot,width,width[0],axis_unit)
				plot = yt.OffAxisProjectionPlot(data_source.ds,basis_vectors[0],fields,center=container.center,width=width,depth=width[0],weight_field=weight_field, north_vector=north_vectors[0])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])
					plot.annotate_marker((4.0,6.92820), coord_system='plot',plot_args={'color':'black','s':400})

				image = plot.frb
				
				plots["0_plot"] = plot
				frb["0_frb"] = image

				if plot_images:
					plot_viz(plot,name + "_axis_0",suffix)

			if 1 in axis:
				print "plotting axis 1"
				data_source = gen_data_source(1,container,raw_snapshot,width,width[0],axis_unit)
				plot = yt.OffAxisProjectionPlot(data_source.ds,basis_vectors[1],fields,center=container.center,width=width,depth=width[0],weight_field=weight_field, north_vector=north_vectors[0])

				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])


				image = plot.frb
				
				plots["1_plot"] = plot
				frb["1_frb"] = image

				if plot_images:
					plot_viz(plot,name + "_axis_1",suffix)

			if 2 in axis:
				print "plotting axis 2"
				data_source = gen_data_source(2,container,raw_snapshot,width,width[0],axis_unit)
				plot = yt.OffAxisProjectionPlot(data_source.ds,basis_vectors[2],fields,center=container.center,width=width,depth=width[0],weight_field=weight_field, north_vector=north_vectors[0])


				# set the units
				plot.set_axes_unit("kpc")
				for i in range(0,len(fields)):
					plot.set_unit(fields[i],units[i])
					plot.set_zlim(fields[i],extrema[i][0],extrema[i][1])

				image = plot.frb
				
				plots["2_plot"] = plot
				frb["2_frb"] = image

				if plot_images:
					plot_viz(plot,name + "_axis_2",suffix)




	else:
		print "module not defined, please either use YT, pynbody or pymses"
		return

	if module == "pymses":
		from pymses.analysis.visualization import *
		from ramses_pp.modules.pymses import PymsesProjection

	if return_objects == True and return_fields == True:
		print plots, frb, fields
		return plots, frb, fields

	if return_objects == True:
		return plots, frb, None

	if return_fields == True:
		return None, None, fields
