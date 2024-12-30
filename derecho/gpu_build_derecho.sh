#!/bin/bash -l

#PBS -N MURbui 
#PBS -A UCSU0085 
#PBS -q main@gusched01
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB:ngpus=1
##PBS -l gpu_type=a100
#PBS -l walltime=00:10:00
#PBS -e build.err 
#PBS -o build.out 

module purge
module load ncarenv
module load nvhpc/24.3
module load cuda
module load craype
module load cray-mpich

module load ncarcompilers
module list

cp Make_defs.derecho_gpu Make_defs

make clean
make

 
