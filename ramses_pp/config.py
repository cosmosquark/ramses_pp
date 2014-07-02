'''
Module config object to store app wide settings
'''

#Enable/disable modules
quick_import = True

yt_enabled = False
pymses_enabled = True
pynbody_enabled = True

default_module='pynbody'

ramses_f90_dir='../ramses_f90'
json_dir = '../data'
override = True
verbose = True

#Methods

def list_modules():
	print 'Enabled modules:'
	print 'yt: ', yt_enabled
	print 'pymses: ', pymses_enabled
	print 'pynbody: ', pynbody_enabled
