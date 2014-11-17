#!/usr/bin/python
from ramses_pp import config
##########################
## Author, Ben Thompson
## just some simple routines to get AHF running it's built in analysis routines
## these get loaded into the Simulation Snapshot Object if the pynbody module is available
########################
import sys, os
import json, glob, uuid
import numpy as np
import re


class AHF():
	@staticmethod
	def run_ahf(self,gas=False):
		"""
		run AHF on all pynbody snapshots
		"""
		for i in range(1,self.num_snapshots()):
			shot = self.snapshot(i, module="pynbody")
			shot.tipsy(convert=True,gas=gas)
			shot.halos()
		return

	@staticmethod
	def run_ahf_merger(self):
		"""
		Generate the Halo merger tree from AHF outputs
		"""
		snapno = self.num_snapshots()
		ahf_files = []
		red = []
		mfiles = []
		for isnap in range(snapno,0,-1):
			tipsy_dir = str("%s/output_%05d/output_%05d_tipsy/" % (self.path(), isnap, isnap))
			if not os.path.isdir(tipsy_dir):
				print "AHF not run on this snapshot.. aborting"
				return
			if len(glob.glob('%s*.AHF_particles'%tipsy_dir)) == 0:
				print "AHF not run on this snapshot.. aborting"
				return
			ahf_files.append(glob.glob('%s*.AHF_particles'%tipsy_dir)[0])
			red.append( re.split(r'\.(?!\d)', glob.glob('%s*.AHF_particles'%tipsy_dir)[0])[2] )
			
		for m in range(snapno,1,-1):
			mfiles.append("%s/output_%05d/output_%05d_tipsy/output_%05d_fullbox.tipsy.%s.AHF" % (self.path(), m, m, m,str(red[snapno-m])))
		
		f = open((config.root_dir + "/applications/ahf_mtree_infile"), 'w')
		f.write(str(snapno))
		f.write("\n")
		for line in ahf_files:
			f.write(line)
			f.write("\n")
		for line in mfiles:
			f.write(line)
			f.write("\n")
		f.close()

		execommand = config.root_dir + "/applications/ahf-v1.0-084/bin/MergerTree < " + config.root_dir + "/applications/ahf_mtree_infile"
		os.system(execommand)
		print "cleaning up"
		execommand = "rm " + config.root_dir + "/applications/ahf_mtree_infile"
		os.system(execommand)

	@staticmethod
	def run_ahf_tracker(self,snaptime=None,halos=None,forward=True):
		"""
		take in a halo object from a snapshot of your choice  and track all or specific (halos_sel numpy array) individual halos backwards in time !!! WARNING, DO NOT RUN IN PARALLEL WITH ANOTHER AHF HALO TRACKER PROCESS"
		"""

		if snaptime==None:
			snaptime = self.num_snapshots()

		if snaptime > self.num_snapshots():
			print "too many snapshots, returning"
			return

		if halos == None:
			snap = self.snapshot(snaptime,module="pynbody")
			halo_arr = snap.halos()
			count = len(halo_arr)
			halos = np.arange(0,count)
		
		ahf_files = []
		mfiles = []
		mtree_files = np.zeros(snaptime-1, dtype=str)
		for isnap in range(snaptime,0,-1):
			tipsy_dir = str("%s/output_%05d/output_%05d_tipsy/" % (self.path(), isnap, isnap))
			ahf_files.append(glob.glob('%s*.AHF_halos'%tipsy_dir))
			if not os.path.isdir(tipsy_dir):
				print "AHF not run on this snapshot.. aborting"
				return
			if len(glob.glob('%s*.AHF_halos'%tipsy_dir)) == 0:
				print "AHF not run on this snapshot.. aborting"
				return


		for i in range (0, len(ahf_files)):
			word = ahf_files[i]
			word[0] = word[0][:-6]
			ahf_files[i] = word # trims the text "halos" from the main word... to give a suitable prefix

		# write the halo id file
		f = open((config.root_dir + "/applications/ahf_track_halos"),'w')
		for i in range(0,len(halos)):
			f.write(str(halos[i]) + "\n")
		f.close()

		f = open((config.root_dir + "/applications/ahf_track_infile"), 'w')

		## write the track infile
		for i in range(0,snaptime):
			f.write(str(ahf_files[i][0]) + "\n")
		f.close()

		# write the redshift file

		print snaptime
		print self.info_snap(snaptime)['z']
		red = self.avail_redshifts(zmin=self.info_snap(snaptime)['z']) # redshifts based on RAMSES simulation redshifts here
		print red
		print len(red)
		if red[(snaptime-1)] < 0.0000000000: # sometimes, the final snapshot can be less than 0
			red[(snaptime-1)] = 0

		z = open((config.root_dir + "/applications/zfile"),'w')
		for i in range(snaptime-1,-1,-1): # technically since we're going back in time.. this actually works index wise
			z.write(str(red[i]) + "\n")
		z.close()

		execommand =  config.root_dir + "/applications/ahf-v1.0-084/bin/ahfHaloHistory " + config.root_dir + "/applications/ahf_track_halos " + config.root_dir + "/applications/ahf_track_infile " + config.root_dir + "/applications/zfile"
		os.system(execommand)
		
		# finally storing the halo tracker data from this timestep into the timestep directory
		print "cleaning up"
		execommand = "rm " + config.root_dir + "/applications/ahf_track_halos" 
		os.system(execommand)
		execommand = "rm " + config.root_dir + "/applications/ahf_track_infile"
		os.system(execommand)
		execommand = "rm " + config.root_dir + "/applications/zfile"
		os.system(execommand)
		
		for h in range(0,len(halos)):
			path = ("%s/output_%05d/output_%05d_tipsy/ahf_halo_track/" %  (self.path(), snaptime, snaptime))
			if not os.path.isdir(path): os.mkdir(path)
			fname = ("%s/halo_%07d.dat" % (os.getcwd(), halos[h]))
			execommand = "mv " + fname + " " + path
			os.system(execommand)
