#include <mpi.h>
#include <iostream>
#include <string.h>
#include "run.H"

using namespace std;

RunData::RunData() {
  rank = 0;
  zrank = 0;
  autostep = 0;
  iteration = 0;
  globiter = 0;
  maxiter = 1000000;
  anlfreq = maxiter;
  resfreq = maxiter;
  backfreq = maxiter;
  outcad = 0.0;
  slicefreq = maxiter;
  dt = 0.;
  time = 0.;
  Tmax = 1.0;
  CFL = 0.8;
  CFL_tvd = 0.95;
  maxWtime = 1.0e7;
  verbose = 0;
 
  // Defaults save temp and pres, set rest to 0
  eos_output[0] = 1;
  eos_output[1] = 1;
  for (int ee=2;ee<13;ee++) eos_output[ee] = 0;
  for (int tt=0;tt<11;tt++) diag_output[tt] = 0;

  
  // Flag to turn on/off diagnostic saves (Tvar-related)
  diagnostics = 0;

  // Only run diagnostics on the last step of the RK. to be saved first step of the next iteration.

  need_diagnostics = 0;

  HAVG = 0;
  DEM = 0;
  RT_HAVG = 0;

  strcpy(resfile,"result");
  strcpy(anlfile,"anls.log");
  strcpy(backfile,"backup.dat");

  strcpy(eos_name,"eos.dat");
  strcpy(kap_name,"kappa.dat");

  strcpy(comment,"");
  strcpy(path_3D,"");
  strcpy(path_2D,"");  
}

void RunData::Synchronize(real delt) {
  MPI_Allreduce(&delt,&dt,1,REALTYPE,MPI_MIN,MPI_COMM_WORLD);
  dt *= CFL;
}

void RunData::Show() const {
  if(rank ==0) {
    cout << " ------------ Run Parameter Settings -------------" << endl;
    cout << "rank        = " << rank << endl
	 << "autostep    = " << autostep << endl
	 << "comment     = " << comment << endl
	 << "anlfile     = " << anlfile << endl
	 << "resfile     = " << resfile << endl
	 << "backfile    = " << backfile << endl
	 << "eos_name    = " << eos_name << endl
	 << "kap_name    = " << kap_name << endl
	 << "iteration   = " << iteration << endl
	 << "globiter    = " << globiter << endl
	 << "maxiter     = " << maxiter << endl
	 << "anlfreq     = " << anlfreq << endl
	 << "resfreq     = " << resfreq << endl
	 << "backfreq    = " << backfreq << endl
         << "outcad      = " << outcad << endl
	 << "slicefreq   = " << slicefreq << endl
	 << "dt          = " << dt << endl
	 << "time        = " << time << endl
	 << "Tmax        = " << Tmax << endl
	 << "CFL         = " << CFL << endl
	 << "CFL_tvd     = " << CFL_tvd << endl
	 << "maxWtime    = " << maxWtime << endl
	 << "verbose     = " << verbose << endl
	 << "path_3D     = " << path_3D << endl
	 << "path_2D     = " << path_2D << endl
	 << "DEM         = " << DEM << endl
	 << "HAVG        = " << HAVG << endl
	 << "RT_HAVG     = " << RT_HAVG << endl
	 << "diagnostics = " << diagnostics << endl
         << "eos_output  = " << eos_output[0] << " " << eos_output[1] << " " << eos_output[2]
	 << " " << eos_output[3] << " " << eos_output[4] << " " << eos_output[5]
	 << " " << eos_output[6] << " " << eos_output[7] << " " << eos_output[8]
	 << " " << eos_output[9] << " " << eos_output[10] << " " << eos_output[11]
	 << " " << eos_output[12] << " " << eos_output[13] << endl
	 << "diagnostics = " << diagnostics << endl
         << "diag_output = " << diag_output[0] << " " << diag_output[1] << " " << diag_output[2]
	 << " " << diag_output[3] << " " << diag_output[4] << " " << diag_output[5]
	 << " " << diag_output[6] << " " << diag_output[7] << " " << diag_output[8]
	 << " " << diag_output[9] << " " << diag_output[10] << endl;
  }
}
