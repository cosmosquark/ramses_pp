#from mpi4py import MPI
import yt
yt.enable_parallelism()
from yt.config import ytcfg
from yt.analysis_modules.halo_finding.rockstar.api import \
	RockstarHaloFinder
from yt.data_objects.particle_filters import \
	particle_filter
from yt.data_objects.time_series import DatasetSeries
import numpy as np

#ytcfg.set('yt', '__global_parallel_size', str(MPI.COMM_WORLD.Get_size()))
ncpu = 32
ytcfg.set('yt', '__global_parallel_size', str(ncpu))

# create a particle filter to remove star particles
@yt.particle_filter("dark_matter", requires=["particle_age"])
def _dm_filter(pfilter, data):
	return data["particle_age"] == 0.0

def setup_ds(ds):
	ds.add_particle_filter("dark_matter")

#from yt.data_objects.particle_filters import add_particle_filter
#add_particle_filter("dark_matter", function=DarkMatter, filtered_type='all', requires=["particle_age"])

start = 30
end = 215
outputs = np.arange(start, end+1)
dirs = []
for ioutput in outputs:
	dirs.append('../output_%05d/info_%05d.txt'%(ioutput, ioutput))

#es = yt.load('../output_*/info_*.txt')
es = DatasetSeries(dirs, setup_function=setup_ds)

writers = 20
readers = 11
rh = RockstarHaloFinder(es, num_readers=readers, num_writers=writers,
			particle_type="dark_matter")
rh.run()

#for ds in es:
#	dd = ds.all_data()
#	print dd.particles.source['particle_age']
