# MURaM Automated Debugging (PCAST)

PGI PCAST features are found in include/muramacc.H. These features are disabled for automatically for other compilers. You do not need to run GPU code to take advantage of PCAST; you can compare two CPU codes with PCAST as long as both are using the PGI compiler. To enable these debug features you must add the following line to the Make_defs file:

DGB = -DDEBUG -DPGICOMPARE

PCAST works by first doing a "golden run" to generate reference results. For this, it's recommended to use the master branch. Build the master branch with the above flags. It is also recommended to set a small number of iterations in parameter.dat file, as PCAST will produce larger files the most iterations that are run. Then before you run the code, set the following environment variable:

$ export PGI_COMPARE=create,file=filename.dat

Give it a recognizable filename. If doing MPI runs, you differentiate the filenames by rank. For example in OpenMPI:

$ export PGI_COMPARE=create,file=%q{OMPI_COMM_WORLD_RANK}.dat

After that has finished running, you should see the generated file(s). If you are testing the GPU code, switch to the most recent GPU branch, and rebuild, including the debug flags again in Make_defs. Make sure that you are running the same test case, with the same number of MPI ranks and core layout. Now we want to use the reference file we just generated. There are a lot of options to do this, but this is a typical example:

$ export PGI_COMPARE=summary,rel=1,report=1,patch,filename.dat

The code will run again, and any errors that exceed the given threshold will be printed out with their actual value and their experimental value.

[This is the PGI website for more info, and a list of all of the PGI_COMPARE options.](https://www.pgroup.com/resources/pcast.htm)