#!/bin/bash -l

#SBATCH -A m4093_g
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 0:10:00
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 128
#SBATCH --gpus-per-task=1
#SBATCH -e MURaMBuild.err
#SBATCH -o MURaMBuild.out

module load gpu
module load PrgEnv-nvidia
module load cray-fftw

export CRAY_ACCEL_TARGET=nvidia80

make clean
make

 
