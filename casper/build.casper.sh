#!/bin/bash -l

#PBS -N MURbui 
#PBS -A 
#PBS -q casper
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB:ngpus=1
#PBS -l gpu_type=v100
#PBS -l walltime=00:10:00
#PBS -e build.err 
#PBS -o build.out 

#cd <project_dir>/MURaM_main

module purge
module load nvhpc/22.2
module load openmpi/4.1.1
module load cuda/11.4.0

make clean
make

 
