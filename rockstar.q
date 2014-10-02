#!/bin/sh
#$ -N yt_rk
#$ -M ds381@sussex.ac.uk
#$ -m bea
#$ -cwd
#$ -pe openmpi_mixed_4 32
#$ -q mps.q@@mps_amd
#$ -S /bin/bash
# source modules environment:
module add sge
module add gcc/4.8.1
module add openmpi/gcc/64/1.7.3

source /mnt/lustre/scratch/phys/ds381/yt_3.0/yt-x86_64/bin/activate
echo "Python set to:"
which python

mpirun -np ${NSLOTS} --mca btl ^openib python ~/.local/lib/python2.6/site-packages/ramses_pp/run_rockstar.py 25Mpc_128_nbody --parallel
