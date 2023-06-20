#!/bin/bash -l

#PBS -N MURstk 
#PBS -A UCSU0085
#PBS -q main@gusched01
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB
#PBS -l walltime=00:50:00
#PBS -e stackbuild.err 
#PBS -o stackbuild_heffte.out 

module purge
module load ncarenv/22.12
module load nvhpc/22.11
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cmake
module list

# By default, script will install fftw and heffte libraries in current location
# Modify the INSTALL_DIR path to install in a different location
export INSTALL_DIR=$(pwd)
export FFTW_INSTALL_DIR=$INSTALL_DIR/fftw
export HEFFTE_INSTALL_DIR=$INSTALL_DIR/heffte
export MPICH_GPU_SUPPORT_ENABLED=1

# download FFTW library
cd $INSTALL_DIR
#wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvzf fftw-3.3.10.tar.gz
cd fftw-3.3.10

# Build single and double precision versions 
./configure CC=mpicc CFLAGS="-O3 -fPIC" --prefix=$FFTW_INSTALL_DIR --enable-threads --enable-mpi

make
make install
make clean

./configure CC=mpicc CFLAGS="-O3 -fPIC" --prefix=$FFTW_INSTALL_DIR --enable-threads --enable-mpi --enable-float
make
make install

# download heffte library
#wget https://bitbucket.org/icl/heffte/get/5caa90cf028e.zip
#mv heffte.zip heffte_build.zip
#unzip heffte_build.zip

cd icl-heffte-5caa90cf028e
mkdir build
cd build

export CRAY_ACCEL_TARGET=nvidia80

#export MPI_CXX=$MPI_ROOT/bin/mpic++


cmake \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG" \
-D BUILD_SHARED_LIBS=OFF \
-D CMAKE_INSTALL_PREFIX=$HEFFTE_INSTALL_DIR \
-D Heffte_ENABLE_FFTW=ON \
-D FFTW_ROOT=$FFTW_INSTALL_DIR \
-D Heffte_ENABLE_CUDA=ON \
-D CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
-D CMAKE_CXX_COMPILER=$NCAR_ROOT_NCARCOMPILERS/bin/nvc++ \
-D MPI_CXX_COMPILER=$NCAR_ROOT_NCARCOMPILERS/bin/mpi/mpic++ \
$INSTALL_DIR/icl-heffte-5caa90cf028e

make
make install

cd $INSTALL_DIR
rm -rf icl-heffte-5caa90cf028e
rm heffte_build.zip



