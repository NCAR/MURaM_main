/* 
 *   Copyright (C) 1998-2000 University of Chicago. 
 *   See COPYRIGHT notice in top-level directory.
 */

#include <mpi.h>
#include "run.H"
#include "grid.H"
#include <stdlib.h>
#include <stdio.h>
#include <fftw3-mpi.h>
#include "physics.H"
#include "solver.H"
#include "init.H"
#include "icbc.H"
#include "io.H"
#include "comm_split.H"
#include "rt/rt.h"
#include <iostream>

using std::cout;
using std::endl;

double start_time,rst_time;

int main(int argc, char** argv) {

  RunData Run;
  GridData Grid;
  PhysicsData Physics;
  
  RTS * rts = 0;
  
  double clock;

#ifdef MURAM_FFTW
  int threads_ok;
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED, &provided);
  threads_ok = provided >= MPI_THREAD_FUNNELED;
  if (threads_ok)
    threads_ok = fftw_init_threads();

  fftw_mpi_init();

  int nthreads = atoi(getenv("MURaM_FFTW_THREADS"));  
  if (nthreads <=0)
    nthreads = 1;
  if (threads_ok)
    fftw_plan_with_nthreads(nthreads);
#else
  MPI_Init(&argc,&argv);
#endif

  start_time=MPI_Wtime();
 
  if( Initialize(Run,Grid,Physics,rts)==restart ) {
    clock=MPI_Wtime();
    RestoreSolution(Run,Grid,Physics); 
    rst_time = MPI_Wtime()-clock;
    if (Run.rank == 0){
      cout << "Restored in " << rst_time << " seconds" << endl;
    }
  } else {
    if (Run.rank == 0) cout << Run.backfile << " not found. Aborting ... " << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  ComputeSolution(Run,Grid,Physics,rts);

  comm_split_finalize();
  IO_Finalize();
  int ierr = MPI_Finalize();

  if(rts) delete rts;
  if (Run.rank==0) fprintf(stderr,"deleted rts\n");
       
  return 0;
}
