'''
Based on Snapshot.py from the Hamu project https://github.com/samgeen/Hamu

@author dsullivan
'''

import abc

#Abstract projection class

class Profile():
	'''
	Abstract class for projecting a single ramses snapshot. Pymses/Pynbody specific classes to be based on this class.
	'''
	__metaclass__ = abc.ABCMeta

	def __init__(self, snapshot, sim_type):
		'''
		Constructor
		'''
		self._snapshot = snapshot
		self._type = sim_type