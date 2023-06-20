#!/bin/bash -l

# Setup Environment ###############

SCRIPT_DIR=$(pwd)
OUTPUT_DIR=$SCRIPT_DIR
C_COMPILER=pgcc

module purge

# Intel
#module load intel
#module load mpt

# NVHPC
module load nvhpc/21.9
module load openmpi/4.1.1
###################################

wget -nc ftp://fftw.org/pub/fftw/fftw-3.3.10.tar.gz
tar -x --skip-old-files -f fftw-3.3.10.tar.gz
cd fftw-3.3.10

./configure CC=$C_COMPILER CFLAGS="-O3" --prefix=$OUTPUT_DIR --enable-mpi --enable-threads --enable-float
make
make install
make clean
./configure CC=$C_COMPILER CFLAGS="-O3" --prefix=$OUTPUT_DIR --enable-mpi --enable-threads
make
make install

cd ..

