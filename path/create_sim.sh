#!/usr/bin/env python

import sys, os

if (len(sys.argv) != 3):
	print 'Usage: python %s <name> <full_path>'%sys.argv[0]
	sys.exit(1)

from modules import Simulation

name = sys.argv[1]
#The fully qualified path
path = sys.argv[2]

if path == '.':
	path = os.getcwd()

sim = Simulation.create(name, path)
print 'Created simulation: %s'%sim.name()
print sim
