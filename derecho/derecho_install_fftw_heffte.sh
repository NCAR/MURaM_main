#!/bin/bash -l

#PBS -N MURstk 
#PBS -A P22100000 
#PBS -q main
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB
#PBS -l walltime=00:50:00
#PBS -e stackbuild.err 
#PBS -o stackbuild_heffte.out 

module purge
module load ncarenv
module load nvhpc/24.3
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cmake
module list

# By default, script will install fftw and heffte libraries in current location
# Modify the INSTALL_DIR path to install in a different location
export INSTALL_DIR=$(pwd)
export FFTW_INSTALL_DIR=$INSTALL_DIR/fftw_24.3
export HEFFTE_INSTALL_DIR_CPU=$INSTALL_DIR/heffte_cpu_24.3
export HEFFTE_INSTALL_DIR_GPU=$INSTALL_DIR/heffte_gpu_24.3
export MPICH_GPU_SUPPORT_ENABLED=0

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

cd $INSTALL_DIR
rm -f fftw-3.3.10

# download heffte library
#wget https://bitbucket.org/icl/heffte/get/5caa90cf028e.zip
#mv heffte.zip heffte_build.zip
#unzip heffte_build.zip

unzip 5caa90cf028e.zip
cd icl-heffte-5caa90cf028e
mkdir build
cd build

cmake \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG" \
-D BUILD_SHARED_LIBS=OFF \
-D CMAKE_INSTALL_PREFIX=$HEFFTE_INSTALL_DIR_CPU \
-D Heffte_ENABLE_FFTW=ON \
-D FFTW_ROOT=$FFTW_INSTALL_DIR \
-D Heffte_ENABLE_CUDA=OFF \
-D CMAKE_CXX_COMPILER=nvc++ \
-D MPI_CXX_COMPILER=mpic++ \
$INSTALL_DIR/icl-heffte-5caa90cf028e

make
make install

cd $INSTALL_DIR
rm -rf icl-heffte-5caa90cf028e


unzip 5caa90cf028e.zip
cd icl-heffte-5caa90cf028e
mkdir build
cd build

export CRAY_ACCEL_TARGET=nvidia80
export MPICH_GPU_SUPPORT_ENABLED=1

cmake \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG" \
-D BUILD_SHARED_LIBS=OFF \
-D CMAKE_INSTALL_PREFIX=$HEFFTE_INSTALL_DIR_GPU \
-D Heffte_ENABLE_FFTW=ON \
-D FFTW_ROOT=$FFTW_INSTALL_DIR \
-D Heffte_ENABLE_CUDA=ON \
-D CMAKE_CXX_COMPILER=nvc++ \
-D MPI_CXX_COMPILER=mpic++ \
$INSTALL_DIR/icl-heffte-5caa90cf028e

make
make install

cd $INSTALL_DIR
rm -rf icl-heffte-5caa90cf028e

