#!/bin/bash -l

#PBS -N MurJob
#PBS -A P22100000 
#PBS -q main
#PBS -l select=1:ncpus=64:mpiprocs=4:ngpus=4
#PBS -l walltime=00:05:00

#cd $cwd
export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

module purge

module load ncarenv
module load nvhpc/24.3
module load cuda
module load cray-mpich
module load craype
module load ncarcompilers
module load cray-libsci

module list

#stack
#export LIBS_HOME=/glade/work/rempel/MURaM_derecho_asd/FFTW_LIBS/
#export FFTW3_HOME=$LIBS_HOME/fftw_gpu_24.3
#export HEFFTE_HOME=$LIBS_HOME/heffte_gpu_24.3
#export PATH=$FFTW3_HOME/bin:$HEFFTE_HOME/bin:$PATH
#export LD_LIBRARY_PATH=$FFTW3_HOME/lib:$HEFFTE_HOME/lib:$LD_LIBRARY_PATH
#export NVHPC_CUDA_HOME=/glade/u/apps/common/23.04/spack/opt/spack/cuda/11.7.1/

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

#mpiexec -ppn 4 /glade/u/apps/derecho/23.06/opt/utils/bin/get_local_rank ./mhd3d.x > muram_log.$PBS_JOBID
mpiexec -ppn 4  ./mhd3d.x > muram_log.$PBS_JOBID
#cuda-memcheck ./mhd3d.x > MURaM.out
