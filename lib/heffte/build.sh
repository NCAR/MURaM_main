#!/bin/bash -l

# Set up your environment ##########
SCRIPT_DIR=$(pwd)
INSTALL_DIR=$SCRIPT_DIR

FFTW=ON
FFTW_HOME=/glade/scratch/efwright/MURaM_main/lib/fftw3

CUDA=ON
CUDA_HOME=/glade/u/apps/dav/opt/cuda/11.4.0

CXX=pgc++
MPICXX=mpic++


module purge

# Intel CPU
module load cmake/3.18.2
#module load intel
#module load mpt

# OpenACC GPU
module load nvhpc/21.9
module load openmpi/4.1.1
module load cuda/11.4
####################################
git clone https://bitbucket.org/icl/heffte.git
cd heffte
mkdir -p build
cd build

#rm -v !("build.sh")

cmake \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG" \
-D BUILD_SHARED_LIBS=OFF \
-D CMAKE_INSTALL_PREFIX=$INSTALL_DIR \
-D Heffte_ENABLE_FFTW=$FFTW \
-D FFTW_ROOT=$FFTW_HOME \
-D Heffte_ENABLE_CUDA=$CUDA \
-D CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
-D CMAKE_CXX_COMPILER=$CXX \
-D MPI_CXX_COMPILER=$MPICXX \
../

make
make install

cd ../..

