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

default_module = 'yt'

ramses_f90_dir='/home/d/ds/ds381/Code/ramses-rt/trunk/ramses/utils/f90/'
json_dir = '/home/d/ds/ds381/.local/lib/python2.6/site-packages/ramses_pp/data'
override = True
verbose = True
rockstar_base='rockstar/full_cat/rockstar_halos/'
rockstar_autogrp=False

#Methods

def list_modules():
	print 'Enabled modules:'
	print 'yt: ', yt_enabled
	print 'pymses: ', pymses_enabled
	print 'pynbody: ', pynbody_enabled