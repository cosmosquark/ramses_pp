#!/usr/bin/env python

import sys

if (len(sys.argv) != 2):
	print 'Usage: python %s <name> '%sys.argv[0]
	sys.exit(1)

from ramses_pp.modules import Simulation

name = sys.argv[1]
sim = Simulation.load(name)
rmdir = False

print "Delete the data directory for the simulation '%s'?"%name
print "Yes/No:"
s = raw_input('--> ')

if (s.lower() == 'yes'):
	print "Are you sure? (Yes/No):"
	s = raw_input('--> ')
	rmdir
	
print rmdir

sim.delete(rmdir)

print 'Deleted simulation: %s'%name
