#######################
How to reference the code:
Please reference the following papers when using this code:

Original code:
Vögler, A., Shelyag, S., Schüssler, M., Cattaneo, F., Emonet, T., Linde, T., "Simulations of magneto-convection in the solar photosphere.  Equations, methods, and results of the MURaM code", Astronomy and Astrophysics, v.429, p.335-351 (2005)

Updated numerical diffusivities:
Rempel, M., "Numerical Simulations of Quiet Sun Magnetism: On the Contribution from a Small-scale Dynamo", The Astrophysical Journal, Volume 789, Issue 2, article id. 132, 22 pp. (2014).

Extended version with coronal physics (the version in this repository is based on this):
Rempel, M., "Extension of the MURaM Radiative MHD Code for Coronal Simulations", The Astrophysical Journal, Volume 834, Issue 1, article id. 10, 23 pp. (2017).

Initial GPU Refactoring:
Wright, E., Przybylski, D., Rempel, M., Miller, C., Suresh, S., Su, S., Loft, R., Chandrasekaran, S., "Refactoring the MPS/University of Chicago Radiative MHD(MURaM) Model for GPU/CPU Performance Portability Using OpenACC Directives", eprint arXiv:2107.08145


#######################
Libraries:
The code requires a parallel version of fftw3 and heffte! 
See the install scripts in lib/fftw and lib/heffte for further details, 

#######################
Compilation:
Customize the file "Make_defs". In general the following varibales need to be set:

MURaM_HOME_DIR = MURaM home directory
FFTW3_HOME = location of fftw3 library
HEFFTE_HOME = location of heffte library

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
