'''
Module config object to store app wide settings
'''


# TODO - Rewrite this entirely (nice to have a system like pynbody where you can
#	do config_parser.getboolean('Rockstar', "some_prop"))

#Enable/disable modules
yt_enabled = True
pymses_enabled = True
pynbody_enabled = True
quick_import = True

# core settings
default_module='yt'

ramses_f90_dir='/gpfs/home/bthompson1/ramses_pp/ramses_pp/ramses_f90'
json_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp/data'
simulation_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp/simulations'
applications_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp/applications'
root_dir = '/gpfs/home/bthompson1/ramses_pp/ramses_pp'

# void settings
vide_catalogue_root = '/starpulse/pool/cosmic/bthompson/voids' # the directory where you run VIDE
void_finder = "vide"


rockstar_base=""
rockstar_autogrp=False
rockstar_load_ctree=False
default_finder="AHF"

override = True
verbose = True

#Methods

def list_modules():
	print 'Enabled modules:'
	print 'yt: ', yt_enabled
	print 'pymses: ', pymses_enabled
	print 'pynbody: ', pynbody_enabled
