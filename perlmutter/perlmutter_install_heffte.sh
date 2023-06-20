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


#module purge
#module load ncarenv/1.3
#module load nvhpc/22.2
#module load ncarcompilers/0.5.0
#module load openmpi/4.1.1
#module load cuda/11.4.0
module purge
module load PrgEnv-nvidia
module load craype-x86-milan
module load cray-fftw/3.3.8.13
module load cudatoolkit
module list

export CRAY_ACCEL_TARGET=nvidia80

# By default, script will install fftw and heffte libraries in current location
# Modify the INSTALL_DIR path to install in a different location
export INSTALL_DIR=$(pwd)
export FFTW_DIR=/opt/cray/pe/fftw/3.3.8.13/x86_milan/
export HEFFTE_DIR=$INSTALL_DIR/heffte



# download FFTW library
cd $INSTALL_DIR

# download heFFTe library

#wget https://bitbucket.org/icl/heffte/get/5caa90cf028e.zip
#mv 5caa90cf028e.zip heffte_build.zip
#unzip heffte_build.zip

cd icl-heffte-5caa90cf028e
mkdir build
cd build

module load cmake

#export MPI_CXX=$MPI_ROOT/bin/mpic++


export CUDA_cufft_LIBRARY=/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/math_libs/lib64/libcufft.so

cmake \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG" \
-D BUILD_SHARED_LIBS=OFF \
-D CMAKE_INSTALL_PREFIX=$HEFFTE_DIR \
-D Heffte_ENABLE_FFTW=ON \
-D FFTW_ROOT=$FFTW_DIR \
-D Heffte_ENABLE_CUDA=ON \
-D CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
-D CMAKE_CXX_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/bin/nvc++ \
-D MPI_CXX_COMPILER=CC \
-D CUDA_cufft_LIBRARY=/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/math_libs/lib64/libcufft.so \
$INSTALL_DIR/icl-heffte-5caa90cf028e

make
make install

#cd $INSTALL_DIR
#rm -rf icl-heffte-5caa90cf028e
#rm heffte_build.zip 
