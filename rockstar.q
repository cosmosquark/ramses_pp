#!/bin/sh
#$ -N aton_256
#$ -M ds381@sussex.ac.uk
#$ -m bea
#$ -cwd
#$ -pe openmpi 16
#$ -q mps_gpu.q
#$ -S /bin/bash
# source modules environment:
module add sge
module add gcc/4.8.1
module add cuda/5.5
module add openmpi/gcc/64/1.7.3

mpirun -np 16 /home/ds381/Code/ramses_cudaton/trunk/ramses/bin/ramses3d $1 > run.log

