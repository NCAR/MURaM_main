#!/bin/bash -l

##PBS -S /glade/u/apps/dav/opt/nvidia-mps/mps_bash
#PBS -N MUR111 
#PBS -A 
#PBS -q casper
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB:ngpus=1
#PBS -l gpu_type=v100
#PBS -l walltime=00:30:00
#PBS -e MURaM.err 
#PBS -o MURaM.out 
### Set TMPDIR as recommended

#cd $cwd
export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

module purge
module load nvhpc/22.2
module load openmpi/4.1.1
module load cuda/11.4.0
module list


#fftw
export FFTW3_HOME=/glade/work/cmille73/nvhpc222/fftw_339
export HEFFTE_HOME=/glade/work/cmille73/nvhpc222/heffte
export PATH=$FFTW3_HOME/bin:$HEFFTE_HOME/bin:$PATH
export LD_LIBRARY_PATH=$FFTW3_HOME/lib:$HEFFTE_HOME/lib:$LD_LIBRARY_PATH

#export FFTW_PLAN_WITH_NTHREADS=1
export MURaM_FFTW_THREADS=2

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

export CUDA_LAUNCH_BLOCKING=0
export UCX_TLS=rc,sm,cuda_copy,cuda_ipc,gdr_copy
export OMPI_MCA_pml=ucx
export OMPI_MCA_btl=self,vader,tcp,smcuda #openib
export UCX_RNDV_SCHEME=get_zcopy
export UCX_RNDV_THRESH=0
export UCX_MAX_RNDV_RAILS=1
export UCX_MEMTYPE_CACHE=n

mpirun -n 1 ./mhd3d.x > MURaM.out

