'''
Module config object to store app wide settings
'''

#Enable/disable modules
quick_import = True

yt_enabled = False
pymses_enabled = True
pynbody_enabled = True

ramses_f90_dir='/home/d/ds/ds381/Code/ramses-rt/trunk/ramses/utils/f90/'
json_dir = '/home/d/ds/ds381/.local/lib/python2.6/site-packages/ramses_pp/data'
override = True
verbose = True

#Methods

def list_modules():
	print 'Enabled modules:'
	print 'yt: ', yt_enabled
	print 'pymses: ', pymses_enabled
	print 'pynbody: ', pynbody_enabled
