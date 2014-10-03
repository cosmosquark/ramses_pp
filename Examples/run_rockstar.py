from mpi4py import MPI
import yt, sys, os
yt.enable_parallelism()
from yt.config import ytcfg
from yt.analysis_modules.halo_finding.rockstar.api import \
	RockstarHaloFinder
from yt.data_objects.particle_filters import \
	particle_filter
from yt.data_objects.time_series import DatasetSeries
import numpy as np
from ramses_pp.modules import Simulation

#ytcfg.set('yt', '__global_parallel_size', str(MPI.COMM_WORLD.Get_size()))
#ncpu = 32
ncpu = 1
#ncpu = MPI.COMM_WORLD.Get_size()
ytcfg.set('yt', '__global_parallel_size', str(ncpu))

# create a particle filter to remove star particles
@yt.particle_filter("dark_matter", requires=[('particle_mass')])
def dark_matter(pfilter, data):
	if ('all', 'particle_age') in data.ds.field_list:
		return data[("all", "particle_age")] == 0.0
	else:
		arr = np.zeros(len(data['particle_mass']))
		return arr == 0.0

def setup_ds(ds):
	#Return only dark matter particles, and assert that the filter holds
	print 'Here'
	assert(ds.add_particle_filter("dark_matter"))

#from yt.data_objects.particle_filters import add_particle_filter
#add_particle_filter("dark_matter", function=dark_matter, filtered_type='all', requires=["page"])

#Load the simulation
simname = sys.argv[1]
sim = Simulation.load(simname)

#cd to the simulation directory
print 'Changing directory to: %s'%sim.path()
os.chdir(sim.path())
#Check if a rockstar/ directory exists. If it doesn't, make it
if not os.path.exists('%s/rockstar/'%sim.path()):
	os.mkdir('%s/rockstar/'%sim.path())

#cd to rockstar/
os.chdir('rockstar/')
print 'In dir: %s'%os.getcwd()
print 'Starting rockstar...'

outputs = np.arange(1, sim.num_snapshots()+1)
dirs = []
#Add the datasets
for ioutput in outputs:
	#ds = yt.load('../output_%05d/info_%05d.txt'%(ioutput, ioutput))
	ds = sim.snapshot(ioutput, module='yt').raw_snapshot()
	#assert(ds.add_particle_filter("dark_matter"))
	dirs.append(ds)

#es = yt.load('../output_*/info_*.txt')
es = DatasetSeries(dirs, setup_function=setup_ds)
#es = DatasetSeries(dirs)

readers = int(ncpu/4.)
#Reserve one cpu for the server
writers = ncpu - readers - 1
print 'Running rockstar with %i writers and %i readers'%(writers, readers)
rh = RockstarHaloFinder(es, num_readers=readers, num_writers=writers,
			particle_type="dark_matter")
rh.run()
