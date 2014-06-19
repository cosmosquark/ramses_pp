#!/bin/sh
#$ -N r2ascii
#$ -m bea
#$ -cwd
#$ -pe openmpi 16
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge
module add gcc/4.8.1
module add openmpi/gcc/64/1.7.3

mpirun -np $NSLOTS ./ramses_part2ascii $1 $2  > run.log

