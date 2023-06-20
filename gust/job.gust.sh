#!/bin/bash -l

#PBS -N MurJob
#PBS -A UCSU0085 
#PBS -q main
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB:ngpus=1
#PBS -l gpu_type=a100
#PBS -l walltime=00:03:00
#PBS -e run.err 
#PBS -o run.out 

#cd $cwd
export TMPDIR=/glade/gust/scratch/$USER/temp
mkdir -p $TMPDIR

module purge
module load ncarenv/22.12
module load nvhpc/22.11
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module list

#stack
export LIBS_HOME=
export FFTW3_HOME=$LIBS_HOME/fftw
export HEFFTE_HOME=$LIBS_HOME/heffte
export PATH=$FFTW3_HOME/bin:$HEFFTE_HOME/bin:$PATH
export LD_LIBRARY_PATH=$FFTW3_HOME/lib:$HEFFTE_HOME/lib:$LD_LIBRARY_PATH
export NVHPC_CUDA_HOME=/glade/u/apps/gust/22.12/spack/opt/spack/cuda/11.7.1/gcc/7.5.0/

export MURaM_FFTW_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1
export CRAY_ACCEL_TARGET=nvidia80
export CUDA_VISIBLE_DEVICES=0,1,2,3


#export PCAST_COMPARE=rel=7,patchall,summary,file=pcast.dat
#export PCAST_COMPARE=create,file="/glade/gust/scratch/cmille73/muram/nvhpc_2211/MURaM_main/TEST/Test_3D/golden.dat"
#export PCAST_COMPARE=compare,rel=7,patchall,summary,file="/glade/gust/scratch/cmille73/muram/nvhpc_2211/MURaM_main/TEST/Test_3D/golden.dat"

### Source the test environment
echo "job starts"
echo $cwd
echo "PATH: "
echo $PATH
echo "LD_LIBRARY_PATH: "
echo $LD_LIBRARY_PATH

./clean

ulimit -s unlimited
nvidia-smi

mpiexec -n 1 -ppn 1 ./mhd3d.x > MURaM.out
#cuda-memcheck ./mhd3d.x > MURaM.out
