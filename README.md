#######################
Libraries:
The code requires a parallel version of fftw3! Please install fftw3 (www.fftw.org/) with --enable-mpi

#######################
Compilation:
Customize the file "Make_defs". In general the following varibales need to be set:

MURaM_HOME_DIR = MURaM home directory
FFTW3_HOME = location of fftw3 library

OPT = optimization flages

DBG = debugging flags

CC     = C compiler
LD     = C++ compiler
CCC    = C++ compiler

compile the code by typing "make" in the MURaM_HOME_DIR
the executable will be found in MURaM_HOME_DIR/mhd3d.x

#######################
Test Problems:

Two test problems are provided that are small enough to run on a single core. To run these problems copy
MURaM_HOME_DIR/src/mhd3d.x into either of the following directories:
TEST/Test_2D
TEST/Test_3D

Change into the test-directory and type "./clean" to reset the test-problem. Run the code with system specific
command, i.e. mpirun -np 2 mhd3d.x

Note: The processors layout needs to be set in parameters.dat -> procs = 2 1 1

#######################
Verify for correctnes:
In order to benchmark against a reference run use the Jupyter notebook TEST/verify_run.ipynb

###################

WARNING - 2D likely does not function!!!!!!!

##################
Known Issues:

NVHPC does not work on CPUs with optimization higher than -O1, it does work if compiled with -acc=multicore 
