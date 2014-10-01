# License:
#   Copyright (C) 2011 Thomas GUILLET, Damien CHAPON, Marc LABADENS. All Rights Reserved.
#
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.
import numpy as N
import tables as T
import re
import pymses
from tempfile import mkstemp
from os import listdir
from os.path import basename, dirname
from time import time
from pymses.utils import constants as C
from pymses.sources.ramses.output import RamsesOutput
from pymses.sources.ramses.octree import CameraOctreeDatasource
from pymses.sources.ramses.filename_utils import search_valid_outputs
from pymses.utils.regions import Box
from pymses.analysis import *
from pymses.analysis.visualization import *
from pymses.filters.cell_to_points import *
from pymses.filters.region_filter import *
from utils import get_appropriate_unit
from ..view.dialogs import *


class RamsesModel:#{{{
	
	tunit_dict = {'yr': C.year, 'Myr': C.Myr, 'Gyr': C.Gyr}

	def __init__(self):
		self.ramses_dir = None
		self.iout_list = None
		self.selected_iout_index = 0
		self.ro = None
		self.cmap = "jet"
		self.source = None
		self.reload_dset = True
		self.reload_dset_RGBA = True
		self.source_camera = None
		
	def reset(self):
		self.source = None
		self.reload_dset = True
		self.dset_loaded_RGBA = False
		self.source_camera = None
	
	def BuildSource(self, cam, fields_to_read=None, ngrid_max=None):
		if fields_to_read!=None:
			# remember fields to read and ngrid_max
			self.fields_to_read = fields_to_read
			self.ngrid_max = ngrid_max
			
		# CameraOctreeDatasource creation
		source = self.ro.amr_source(self.fields_to_read)
		esize = 0.5**(self.ro.info["levelmin"]+1)
		cod = CameraOctreeDatasource(cam, esize, source, \
				ngrid_max=self.ngrid_max, include_split_cells=True)
		self.source = cod.dset
		self.source_camera = cam
		self.reload_dset = True
		self.dset_loaded_RGBA = False
		# Clear map engine dictionary to force creation of new map engines
		# that take into account this source
		self.map_engine_dict.clear()

	def SearchValidOutputs(self, out_dir):#{{{
		"""
		Search for valid RAMSES outputs in the out_dir directory
		Records the ramses_dir and the sorted iout_list if valid outputs are found and return True
		Return False if none is found
		"""
		if out_dir is None:
			if self.ramses_dir is None:
				return False, 0
			else:
				out_dir = self.ramses_dir
		ilist = search_valid_outputs(out_dir)
		ok = (len(ilist) > 0)
		self.selected_iout_index = 0
		if not ok :
			# try to see if up directory is not already a good directory
			out_up_dir = dirname(out_dir)
			ilist = search_valid_outputs(out_up_dir)
			ok = (len(ilist) > 0)
			if ok:
				# Find output number iout
				outName = basename(out_dir)
				iout_regexp = re.compile("[0-9]{5}")
				iout_parse = iout_regexp.findall(outName)
				if len(iout_parse) > 0:
					iout = int(iout_parse[0])
					self.selected_iout_index = ilist.index(iout)
				self.iout_list = ilist
				self.ramses_dir = out_up_dir
		else:
			self.iout_list = ilist
			self.ramses_dir = out_dir
		return ok, self.selected_iout_index
	#}}}

	def GetvalidOutputsList(self):
		return self.iout_list

	def SetRamsesOutput(self, iout_index=None, iout=None):#{{{
		"""
		If iout is in self.iout_list :
		- Sets the chosen RamsesOutput according to the iout number
		- Init. the region size (boxlen)
		- Return True

		Return False if not.
		"""
		if iout_index is not None:
			if not iout_index in range(len(self.iout_list)):
				return False
			self.selected_iout_index = iout_index
			iout = self.iout_list[iout_index]
		else:
			self.selected_iout_index = self.iout_list.index(iout)
		self.ro = RamsesOutput(self.ramses_dir, iout)
		self.freeTheMemory()
		self.SetRegionSizeUnit(self.ro.info["unit_length"])
		return True
	#}}}

	def GetOutputTime(self):#{{{
		t = self.ro.info["unit_time"] * self.ro.info["time"]
		return get_appropriate_unit(t, RamsesModel.tunit_dict)
	#}}}
	
	def BuildPreviewCube(self):
		cam = self.get_los_camera(25)
		bb=cam.get_bounding_box()
		source = self.ro.amr_source(["rho"])
		from pymses.analysis import amr2cube
		t0 = time()
		res = cam.get_required_resolution()
		self.preview_cube = amr2cube(source, "rho", bb.min_coords, bb.max_coords, res)
		print "Preview amr2cube time = %.3fs"%(time() - t0), "size =",\
			self.preview_cube.shape, "octree level =", res
	
	def FindGalaxyAxis(self, gui):
		cam = self.get_los_camera(256)
		# get the rho map engine :
		# this code comes from model.region_finder.UpdateRegionFinderMaps()
		# TODO : regroup this duplicated code !
		me_key="density"
		if not self.map_engine_dict.has_key(me_key):
			self.map_engine_dict[me_key] = map_engine(self,
							me_key, self.ro, self.source)
		me = self.map_engine_dict[me_key]
		# get the map processor :
		if me.mp == None:
			me.mp = fft_projection.MapFFTProcessor(me.source,
				me.ramses_output.info, remember_data=True,
				cache_dset=self.amr_files_cache_dset, pre_flatten=True)
		me.mp.prepare_data(cam, me.field_list) # prepare data with the first camera
		from pymses.analysis import find_galaxy_axis
		cam.los_axis = find_galaxy_axis(me.mp.filtered_source, cam, nbSample=2000)
		self.LoadLosCamera(cam)
		gui.reset()
	
	def FindCenter(self, gui):
		cam = self.get_los_camera(256)
		# get the rho map engine :
		# this code comes from model.region_finder.UpdateRegionFinderMaps()
		# TODO : regroup this duplicated code !
		me_key="density"
		if not self.map_engine_dict.has_key(me_key):
			self.map_engine_dict[me_key] = map_engine(self,
							me_key, self.ro, self.source)
		me = self.map_engine_dict[me_key]
		# get the map processor :
		if me.mp == None:
			me.mp = fft_projection.MapFFTProcessor(me.source,
				me.ramses_output.info, remember_data=True,
				cache_dset=self.amr_files_cache_dset, pre_flatten=True)
		me.mp.prepare_data(cam, me.field_list) # prepare data with the first camera
		from pymses.analysis import find_center_of_mass
		cam.center = find_center_of_mass(me.mp.filtered_source, cam, nbSample=1E5)
		self.LoadLosCamera(cam)
		gui.reset()
	
	def FindLos(self, gui):
		cam = self.get_los_camera(256)
		# get the rho map engine :
		# this code comes from model.region_finder.UpdateRegionFinderMaps()
		# TODO : regroup this duplicated code !
		me_key="velocity"
		if not self.map_engine_dict.has_key(me_key):
			self.map_engine_dict[me_key] = map_engine(self,
							me_key, self.ro, self.source)
		me = self.map_engine_dict[me_key]
		# get the map processor :
		if me.mp == None:
			me.mp = fft_projection.MapFFTProcessor(me.source,
				me.ramses_output.info, remember_data=True,
				cache_dset=self.amr_files_cache_dset, pre_flatten=True)
		me.mp.prepare_data(cam, me.field_list) # prepare data with the first camera
		from pymses.analysis import find_los
		cam.los_axis = find_los(me.mp.filtered_source, cam, nbSample=1E5)
		self.LoadLosCamera(cam)
		gui.reset()
#}}}


class MapEngine(object):#{{{
	def __init__(self, model, field_list, ramses_output, use_stars_particle_source=False,
		     use_dm_particle_source=False, source=None):
		"""
			Turn the multiprocessing option to True to compute the region finder 
			maps in parallel only if there is a lot of memory availale on the machine
		"""
		self.model = model
		self.field_list = field_list
		self.ramses_output = ramses_output
		if use_stars_particle_source:
			self.source = self.ramses_output.particle_source(self.field_list, select_dark_matter=False)
		elif use_dm_particle_source:
			self.source = self.ramses_output.particle_source(self.field_list, select_stars=False)
		elif source == None and self.ramses_output != None:
			self.source = self.ramses_output.amr_source(self.field_list)
		else:
			self.source = source
		self.mp = None
		self.use_camera_lvlmax = True
		self.log_sensitive = True
		self.fraction = 0.1

	def MakeMaps(self, cam_dict, cmap, multiprocessing, FFTkernelSizeFactor, random_shift,
			cache_dset, stars_age_instensity_dimming=False):#{{{
		# Map computation
		self.map_h5files={}
		if self.mp == None:
			self.mp = fft_projection.MapFFTProcessor(self.source,
				self.ramses_output.info, remember_data=True,
				cache_dset=cache_dset, use_camera_lvlmax=self.use_camera_lvlmax, pre_flatten=True)
		self.mp.prepare_data(cam_dict.values()[0], self.field_list) # prepare data with the first camera
		map_name_list = []
		isStarsParticles = self.getMapName()=="stars_particles_Vlos"
		scal_func = self.getMapOperator(cam_dict.itervalues().next())
		if not self.log_sensitive:
			for c in cam_dict.values():
				c.log_sensitive = False
		if multiprocessing:
			try:
				from multiprocessing import Process, Pipe, cpu_count
				if cpu_count() == 1:
					raise(Exception) # don't use multiprocessing if there is only one cpu
				print 'Multiprocessing Activated'
				p={}
				parent_conn={}
				child_conn={}
				def process_cam(self, axis, cam, conn):
					# Map parameters
					mapname = "wx_AMRViewer_%s_%s"%(self.getMapName(),axis)
					fhandle, fname = mkstemp(prefix=mapname, suffix=".h5")
					mname = basename(fname)[:-3]
					tdir = dirname(fname)
					conn.send(fname)
					mmap = self.mp.process(scal_func, cam, surf_qty=self.IsSurfQuantity(),
						FFTkernelSizeFactor=FFTkernelSizeFactor, random_shift=random_shift,
						data_already_prepared=True,
						stars_age_instensity_dimming=(stars_age_instensity_dimming and isStarsParticles))
					mu = self.getMapUnit()
					au = self.ramses_output.info["unit_length"]
					self.map_h5files[axis] = save_map_HDF5(mmap, cam, unit=mu,
									scale_unit=au, hdf5_path=tdir, map_name=mname)
					conn.send(self.map_h5files[axis])
					conn.close()

				# Start process
				for axis, cam in cam_dict.iteritems():
					parent, child = Pipe()
					parent_conn[axis] = parent
					child_conn[axis] = child
					p[axis] = Process(target=process_cam, args=(self, axis, cam, child_conn[axis]))
					p[axis].start()
				# End Process
				for axis, cam in cam_dict.iteritems():
					map_name_list.append(parent_conn[axis].recv())
					self.map_h5files[axis] = parent_conn[axis].recv()
					p[axis].join()
			except Exception:
				print 'WARNING: multiprocessing unavailable'
				multiprocessing = False
		if not multiprocessing:
			for axis, cam in cam_dict.iteritems():
				# Map parameters
				mapname = "wx_AMRViewer_%s_%s"%(self.getMapName(),axis)
				fhandle, fname = mkstemp(prefix=mapname, suffix=".h5")
				mname = basename(fname)[:-3]
				tdir = dirname(fname)
				map_name_list.append(fname)
				mmap = self.mp.process(scal_func, cam, surf_qty=self.IsSurfQuantity(),
					FFTkernelSizeFactor=FFTkernelSizeFactor, random_shift=random_shift,
					data_already_prepared=True,
					stars_age_instensity_dimming=(stars_age_instensity_dimming and isStarsParticles))
				mu = self.getMapUnit()
				au = self.ramses_output.info["unit_length"]
				self.map_h5files[axis] = save_map_HDF5(mmap, cam, unit=mu,
						scale_unit=au, hdf5_path=tdir, map_name=mname)
		self.mp.free_filtered_source()# Free filtered data
		return map_name_list
	#}}}

	def MakeImage(self, map_name_list_or_img_dict, axis_keys, cmap, adaptive_gaussian_blur,\
			ramses_output=None, verbose = True):
		if isinstance(map_name_list_or_img_dict, list):
			def h5f_iter():
				for map_name in map_name_list_or_img_dict:
					yield map_name
			img_iter = save_HDF5_seq_to_img(h5f_iter, cmap=cmap,\
				adaptive_gaussian_blur=adaptive_gaussian_blur, fraction=self.fraction,\
				ramses_output=ramses_output, verbose=verbose,\
				log_sensitive=self.log_sensitive)
			imglist = [img for img in img_iter()]
			img_dict = dict(zip(axis_keys, imglist))
			return img_dict
		else:
			imglist=[]
			for axis in axis_keys:
				imglist.append(map_name_list_or_img_dict[axis])
			img_dict = dict(zip(axis_keys, imglist))
			return img_dict

	def GetMapArrays(self):#{{{
		"""
		Get the numpy.array of the maps contained in each HDF5 file of
		the self.h5files dict.
		The returned maps are maskred and clipped.
		"""
		ma = {}
		for axis, h5fname in self.map_h5files.iteritems():
			h5f = T.openFile(h5fname, 'r')
			cam = Camera.from_HDF5(h5f)
			map = h5f.getNode("/map/map").read()
			try :
				map_mask = h5f.getNode("/map/map_mask").read()
			except Exception:
				map_mask = numpy.ones(map.shape, "int8")
			h5f.close()
			nx, ny = map.shape
#			map = map.transpose()
#			map_mask = map_mask.transpose()
			if cam.log_sensitive and pymses.__revision__ >= '$Revision: 705$':
				map = apply_log_scale(map)
			vmin, vmax = get_map_range(map[(map_mask>0.)], cam.log_sensitive, None, self.fraction)
			map = N.clip(map[:,::-1], vmin, vmax) # map[:,::-1] <-> map.FLIP_LEFT_RIGHT
			if cam.log_sensitive:
				map = 10.**map
			ma[axis] = map
		return ma
	#}}}

	def getMapName(self):
		raise NotImplementedError()
	
	def getUnitDict(self):
		raise NotImplementedError()

	def getMapUnit(self):
		raise NotImplementedError()

	def getMapOperator(self, cam):
		raise NotImplementedError()

	def IsSurfQuantity(self):
		return False

	
#}}}

class RayTraceEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source, fl=[]):
		self.model = model
		self.ramses_output = ramses_output
		self.source = source
		self.rt = raytracing.RayTracer(ramses_output, fl)
		self.log_sensitive = False
		self.fraction = 0.1
	
	def MakeMaps(self, cam_dict, cmap, multiprocessing, FFTkernelSizeFactor,
			random_shift, cache_dset, stars_age_instensity_dimming=False):#{{{
		# Map computation
		self.map_h5files={}
		map_name_list = []
		axis, cam = cam_dict.keys()[0], cam_dict.values()[0]
		# Map parameters
		mapname = "wx_AMRViewer_%s_%s"%(self.getMapName(),axis)
		fhandle, fname = mkstemp(prefix=mapname, suffix=".h5")
		mname = basename(fname)[:-3]
		tdir = dirname(fname)
		map_name_list.append(fname)
		scal_func = self.getMapOperator(cam)
		mmap = self.rt.process(scal_func, cam, surf_qty=self.IsSurfQuantity())
		mu = self.getMapUnit()
		au = self.ramses_output.info["unit_length"]
		self.map_h5files[axis] = save_map_HDF5(mmap, cam, unit=mu,\
			scale_unit=au, hdf5_path=tdir, map_name=mname)
		return map_name_list
	#}}}

	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "levelmax"
	
	def getMapName(self):
		return "amr_max_ref_level"

	def getMapUnit(self):
		return C.none

	def getUnitDict(self):
		d = {"": C.none}
		return d
			
	def getMapOperator(self, cam):
		op = MaxLevelOperator()
		return op
#}}}

class RayTraceDensityEngine(RayTraceEngine, MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		RayTraceEngine.__init__(self, model, ramses_output, source, fl=["rho"])
		self.log_sensitive = True
	
	def MakeBothMaps(self, cam_dict, me_levelmax):#{{{
		# Map computation
		self.map_h5files={}
		map_name_list = []
		map_name_list_levelmax = []
		axis, cam = cam_dict.keys()[0], cam_dict.values()[0]
		# Map parameters
		mapname = "wx_AMRViewer_%s_%s"%(self.getMapName(),axis)
		mapname_levelmax = "wx_AMRViewer_%s_%s"%(me_levelmax.getMapName(),axis)
		fhandle, fname = mkstemp(prefix=mapname, suffix=".h5")
		fhandle, fname_levelmax = mkstemp(prefix=mapname_levelmax, suffix=".h5")
		mname = basename(fname)[:-3]
		mname_levelmax = basename(fname_levelmax)[:-3]
		tdir = dirname(fname)
		map_name_list.append(fname)
		map_name_list_levelmax.append(fname_levelmax)
		scal_func = self.getMapOperator(cam)
		
		rt = raytracing.OctreeRayTracer(self.source)
		mmap, mmap_levelmax = rt.process(scal_func, cam,
			surf_qty=self.IsSurfQuantity(), rgb=False, use_openCL=self.model.use_openCL,
			dataset_already_loaded=not self.model.reload_dset,
			reload_scalar_field=self.model.dset_loaded_RGBA)
		if self.model.use_openCL:
			self.model.reload_dset=False
			self.model.dset_loaded_RGBA = False
		mu = self.getMapUnit()
		au = self.ramses_output.info["unit_length"]
		self.map_h5files[axis] = save_map_HDF5(mmap, cam, unit=mu,\
			scale_unit=au, hdf5_path=tdir, map_name=mname)
		cam.log_sensitive = False
		mu = me_levelmax.getMapUnit()
		me_levelmax.map_h5files={}
		me_levelmax.map_h5files[axis] = save_map_HDF5(mmap_levelmax, cam, unit=mu,\
			scale_unit=au, hdf5_path=tdir, map_name=mname_levelmax)
		return map_name_list, map_name_list_levelmax, mmap_levelmax
	#}}}
	
	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "densityRT"
	
	def getMapName(self):
		return "gas_densityRT"

	def getMapUnit(self):
		return self.ramses_output.info["unit_density"]

	def getUnitDict(self):
		d = {"H/cc": C.H_cc}
		return d
			
	def getMapOperator(self, cam):
		#num_func = lambda dset: (dset["rho"]**2)
		#denom_func = lambda dset: (dset["rho"])
		#op = FractionOperator(num_func, denom_func)
		op = ScalarOperator(lambda dset: (dset["rho"]))
		op.max_alos=True
		return op

class MassWeightedDensityMapEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		MapEngine.__init__(self, model, ["rho"], ramses_output, source=source)

	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "density"

	def getMapName(self):
		return "gas_mw_density"

	def getMapUnit(self):
		return self.ramses_output.info["unit_density"]

	def getUnitDict(self):
		d = {"H/cc": C.H_cc}
		return d
		
	def getMapOperator(self, cam):
		num_func = lambda dset: (dset["rho"]**2 * dset.get_sizes()**3)
		denom_func = lambda dset: (dset["rho"] * dset.get_sizes()**3)
		return FractionOperator(num_func, denom_func)
#}}}

class SurfaceDensityMapEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		MapEngine.__init__(self, model, ["rho"], ramses_output, source=source)
	
	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "Sigma"
	
	def getMapName(self):
		return "gas_Sigma"

	def getUnitDict(self):
		d = {"Msun/pc^2": C.Msun/C.pc**2}
		return d
	
	def getMapUnit(self):
		u = self.ramses_output.info["unit_density"] * \
				self.ramses_output.info["unit_length"]
		return u
			
	def getMapOperator(self, cam):
		f = ScalarOperator(lambda dset: (dset["rho"] * dset.get_sizes()**3))
		return f
	
	def IsSurfQuantity(self):
		return True
#}}}

class MassWeightedTemperatureMapEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		MapEngine.__init__(self, model, ["P", "rho"], ramses_output, source=source)
	
	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "temperature"
	
	def getMapName(self):
		return "gas_mw_T"
	
	def getUnitDict(self):
		d = {"K": C.K}
		return d

	def getMapUnit(self):
		return self.ramses_output.info["unit_temperature"]
			
	def getMapOperator(self, cam):
		num_func = lambda dset: (dset["P"] * dset.get_sizes()**3)
		denom_func = lambda dset: (dset["rho"] * dset.get_sizes()**3)
		return FractionOperator(num_func, denom_func)
#}}}

# Line-of-sight velocity map engines {{{
class MassWeightedVelocityDispMapEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		MapEngine.__init__(self, model, ["vel", "rho"], ramses_output, source=source)
		self.log_sensitive = False
	
	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "velocity"

	def getMapName(self):
		return "gas_mw_V"
	
	def getUnitDict(self):
		d = {"km/s": C.km/C.s}
		return d
	
	def getMapUnit(self):
		return self.ramses_output.info["unit_velocity"]

	def getMapOperator(self, cam):
		num_func = lambda dset: (-N.sum(cam.los_axis[N.newaxis,:] * dset["vel"], axis=1) *\
				dset["rho"] * dset.get_sizes()**3)
		denom_func = lambda dset: (dset["rho"] * dset.get_sizes()**3)
		return FractionOperator(num_func, denom_func)

#}}}

class DmParticleMassMapEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		MapEngine.__init__(self, model, ["mass", "level"], ramses_output,\
			use_dm_particle_source=True, source=source)
		self.use_camera_lvlmax = False
	
	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "dm_particles"
	
	def getMapName(self):
		return "dm_particles_Vlos"
	
	def getUnitDict(self):
		d = {"Msun/pc^2": C.Msun/C.pc**2}
		return d
	
	def getMapUnit(self):
		u = self.ramses_output.info["unit_mass"] / \
				self.ramses_output.info["unit_length"]**2
		return u

	def getMapOperator(self, cam):
		return ScalarOperator(lambda dset: dset["mass"])
	
	def IsSurfQuantity(self):
		return True
#}}}

class StarsParticleMassMapEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		MapEngine.__init__(self, model, ["mass", "level", "epoch"], ramses_output,\
			use_stars_particle_source=True, source=source)
		self.use_camera_lvlmax = False
	
	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "stars_particles"
	
	def getMapName(self):
		return "stars_particles_Vlos"
	
	def getUnitDict(self):
		d = {"Msun/pc^2": C.Msun/C.pc**2}
		return d
	
	def getMapUnit(self):
		u = self.ramses_output.info["unit_mass"] / \
				self.ramses_output.info["unit_length"]**2
		return u

	def getMapOperator(self, cam):
		return ScalarOperator(lambda dset: dset["mass"])
	
	def IsSurfQuantity(self):
		return True
#}}}

class TransferRayTraceEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		self.model = model
		self.ramses_output = ramses_output
		self.source = source
		self.log_sensitive = False
		self.fraction = 0.1
	
	def MakeMaps(self, cam_dict, cmap, multiprocessing, FFTkernelSizeFactor,
		     random_shift, cache_dset, stars_age_instensity_dimming=False):#{{{
		# Map computation
		img_dict={}
		self.map_array = {}
		if self.source != None:
			rt = raytracing.OctreeRayTracer(self.source)
		else:
			rt = raytracing.OctreeRayTracer(self.ramses_output, ["rho"])
		for axis, cam in cam_dict.iteritems():
			cltf = ColorLinesTransferFunction( (-5.0, 2.0) )
			cltf.add_line(-2.0, 0.1)
			cltf.add_line(.0, 0.1)
			cltf.add_line(2., 0.1)
			cam.set_color_transfer_function(cltf)
			# We add 1e-8 to avoid NaN and -Inf log result
			# problems with approximative null values
			op = ScalarOperator(lambda dset: N.log10(dset["rho"]+1e-8))
			img_dict[axis] = rt.process(op, cam, use_openCL=self.model.use_openCL,
						dataset_already_loaded=not self.model.reload_dset,
						reload_scalar_field=not self.model.dset_loaded_RGBA).convert("RGBA")
			if self.model.use_openCL:
				self.model.reload_dset = False
				self.model.dset_loaded_RGBA = True
			self.map_array[axis] =  N.zeros(cam.get_map_size())
		return img_dict
	#}}}
	
	def GetMapArrays(self):#{{{
		return self.map_array

	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "transfer_function"
	
	def getMapName(self):
		return "transfer_function_map"

	def getMapUnit(self):
		return C.none

	def getUnitDict(self):
		d = {"": C.none}
		return d
#}}}

# Line-of-sight velocity map engines {{{
class CustomMapEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		dlg = FieldsToLoadFileDialog(self)
		dlg.ShowModal()
		s = dlg.fieldTextCtrl.GetValue()
		exec(s)
		MapEngine.__init__(self, model, field_to_load, ramses_output, source=source)
	
	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "custom"

	def getMapName(self):
		return "custom_map"
	
	def getMapUnit(self):
		return C.none

	def getUnitDict(self):
		return {"(?)": C.none}
	
	def getMapOperator(self, cam):
		dlg = customFileDialog(self)
		dlg.ShowModal()
		#self.log_sensitive = dlg.logCheckBox.GetValue()
		return dlg.op

#}}}

class hdf5ViewerEngine(MapEngine):#{{{
	def __init__(self, model, ramses_output, source):
		self.model = model
		self.ramses_output = ramses_output
		self.log_sensitive = False
		self.fraction = 0.1
	
	def GetMapArrays(self):#{{{
		return self.map_array

	@classmethod
	def is_map_engine_for(cls, map_type):
		return map_type == "hdf5"
	
	def getMapName(self):
		return "hdf5_map"
#}}}

def map_engine(model, map_type, ramses_output, source=None):#{{{
	for cls in MapEngine.__subclasses__():
		if cls.is_map_engine_for(map_type):
			return cls(model, ramses_output, source)
	raise ValueError("No MapEngine available for map type '%s'"%map_type)
#}}}
