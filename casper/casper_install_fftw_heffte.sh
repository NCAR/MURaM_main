#!/bin/bash -l

#PBS -N MURstk 
#PBS -A NTDD0004
#PBS -q casper
##PBS -q preview
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB:ngpus=1
#PBS -l gpu_type=v100
#PBS -l walltime=00:50:00
#PBS -e stackbuild.err 
#PBS -o stackbuild.out 


module purge
module load ncarenv/1.3
module load nvhpc/22.2
module load ncarcompilers/0.5.0
module load openmpi/4.1.1
module load cuda/11.4.0

module list


# By default, script will install fftw and heffte libraries in current location
# Modify the INSTALL_DIR path to install in a different location
export INSTALL_DIR=$(pwd)
export FFTW_DIR=$INSTALL_DIR/fftw
export HEFFTE_DIR=$INSTALL_DIR/heffte



# download FFTW library
cd $INSTALL_DIR
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvzf fftw-3.3.10.tar.gz
cd fftw-3.3.10

# Build single and double precision versions 

./configure CC=pgcc CFLAGS="-O3"  --prefix=$FFTW_DIR --enable-threads --enable-mpi

make
make install
make clean

./configure CC=pgcc CFLAGS="-O3"  --prefix=$FFTW_DIR --enable-threads --enable-mpi --enable-float
make
make install

# clean up and remove fftw build directory

cd ..
rm -rf fftw-3.3.10

# download heFFTe library

wget https://bitbucket.org/icl/heffte/get/5caa90cf028e.zip
mv 5caa90cf028e.zip heffte_build.zip
unzip heffte_build.zip

cd icl-heffte-5caa90cf028e
mkdir build
cd build

module load cmake

export MPI_CXX=$MPI_ROOT/bin/mpic++


cmake \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG" \
-D BUILD_SHARED_LIBS=OFF \
-D CMAKE_INSTALL_PREFIX=$HEFFTE_DIR \
-D Heffte_ENABLE_FFTW=ON \
-D FFTW_ROOT=$FFTW_DIR \
-D Heffte_ENABLE_CUDA=ON \
-D CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
-D CMAKE_CXX_COMPILER=/glade/u/apps/opt/nvhpc/22.2/Linux_x86_64/22.2/compilers/bin/nvc++ \
-D MPI_CXX_COMPILER=$MPI_CXX \
$INSTALL_DIR/icl-heffte-5caa90cf028e

make
make install

cd $INSTALL_DIR
rm -rf icl-heffte-5caa90cf028e
rm heffte_build.zip 
