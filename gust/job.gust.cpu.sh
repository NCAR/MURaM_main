#!/bin/bash -l

#PBS -N MurJob
#PBS -A UCSU0085 
#PBS -q main
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB
##PBS -l gpu_type=a100
#PBS -l walltime=01:30:00
#PBS -e run.err 
#PBS -o run.out 

#cd $cwd
export TMPDIR=/glade/gust/scratch/$USER/temp
mkdir -p $TMPDIR

module purge
module load ncarenv/22.12
module load nvhpc/22.11
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module list

#stack
export LIBS_HOME=/glade/gust/scratch/cmille73/nvhpc_2211
export FFTW3_HOME=$LIBS_HOME/fftw
export HEFFTE_HOME=$LIBS_HOME/heffte
export PATH=$FFTW3_HOME/bin:$HEFFTE_HOME/bin:$PATH
export LD_LIBRARY_PATH=$FFTW3_HOME/lib:$HEFFTE_HOME/lib:$LD_LIBRARY_PATH

export MURaM_FFTW_THREADS=1

### Source the test environment
echo "job starts"
echo $cwd
echo "PATH: "
echo $PATH
echo "LD_LIBRARY_PATH: "
echo $LD_LIBRARY_PATH

./clean

ulimit -s unlimited

mpiexec -n 1 -ppn 1 ./mhd3d.x > MURaM.out
