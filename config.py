'''
Module config object to store app wide settings
'''

#Enable/disable modules
quick_import = True

yt_enabled = True
pymses_enabled = True
pynbody_enabled = True

default_module='pynbody'

ramses_f90_dir='/gpfs/home/bthompson1/ramses_pp/ramses_pp/ramses_f90'
json_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp/data'
simulation_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp/simulations'
applications_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp/applications'
root_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp'
vide_catalogue_root = '/starpulse/pool/cosmic/bthompson/voids' # the directory where you run VIDE

rockstar_base=""
rockstar_autogrp=False

override = True
verbose = True

#Methods

def list_modules():
	print 'Enabled modules:'
	print 'yt: ', yt_enabled
	print 'pymses: ', pymses_enabled
	print 'pynbody: ', pynbody_enabled
