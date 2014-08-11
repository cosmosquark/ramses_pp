#!/usr/bin/python
from ramses_pp import config
##########################
## Author, Ben Thompson
## just some simple routines to get AHF running it's built in analysis routines
########################
import sys, os
import json, glob, uuid



class AHF():
	@staticmethod
	def run_ahf(self,gas=True):
		for i in range(1,self.num_snapshots()):
			shot = self.snapshot(i)
			shot.tipsy(gas)
			shot.halos()
		return

	@staticmethod
	def run_ahf_merger(self):
		snapno = self.num_snapshots()
		ahf_files = []
		red = []
		mtree_files = np.zeros(snapno-1, dtype=str)
		for isnap in range(snapno,0,-1):
			tipsy_dir = str("%s/output_%05d/output_%05d_tipsy/" % (self.path(), isnap, isnap))
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

		execommand = "." + config.root_dir + "/applications/ahf-v1.0-084/bin/MergerTree < " + config.root_dir + "applications/ahf_mtree_infile"
		os.system(execommand)

	
	def run_ahf_tracker(self)
		snapno = self.num_snapshots()
		ahf_files = []
		mfiles = []
		mtree_files = np.zeros(snapno-1, dtype=str)
		for isnap in range(snapno,0,-1):
			tipsy_dir = str("%s/output_%05d/output_%05d_tipsy/" % (sim.path(), isnap, isnap))
			ahf_files.append(glob.glob('%s*.AHF_halos'%tipsy_dir))

		for i in range (0, len(ahf_files)):
			word = ahf_files[i]
			word[0] = word[0][:-6]
			ahf_files[i] = word # trims the text "halos" from the main word... to give a suitable prefix

		f = open((config.root_dir + "applications/ahf_track_infile"), 'w')

		## write the track infile
		for i in range(0,snapno):
			f.write(str(ahf_files[i][0]) + "\n")
		f.close()

		# write the redshift file

		red = sim.avail_redshifts() # redshifts based on RAMSES simulation redshifts here
		if red[(sim.num_snapshots()-1)] < 0.00:
			red[(sim.num_snapshots()-1)] = 0

		z = open((config.root_dir + "applications/zfile"),'w')
		for i in range(snapno-1,-1,-1):
			z.write(str(red[i]) + "\n")
		z.close

		execommand = "." + config.root_dir + "/applications/ahf-v1.0-084/bin/ahfHaloHistory " + config.root_dir + "/applications/ahf_track_infile " + config.root_dir + "/applications/zfile"
		os.system(execommand)
