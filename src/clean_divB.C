#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "exchange.H"
#include "muramacc.H"
#include <iostream>
#include <algorithm>

using std::min;
using std::max;
using std::cout;
using std::endl;

/*
  Add contributions from divB heating to Qres for diagnostic reasons.
  Removed ef, since it does not work with the SR version.
*/


#define GLOOP(gs,i,j,d1,d2) \
  for((j)=(gs)[(d2)][0];(j)<=(gs)[(d2)][1];(j)++) \
  for((i)=(gs)[(d1)][0];(i)<=(gs)[(d1)][1];(i)++)

#define XZ_LOOP(G,i,k) \
  for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
  for((i)=(G).lbeg[0];(i)<=(G).lend[0];(i)++)

#define YZ_LOOP(G,j,k) \
  for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
  for((j)=(G).lbeg[1];(j)<=(G).lend[1];(j)++)

void Clean_divB(const RunData&  Run, GridData& Grid, const PhysicsData&  Physics) {

  double ttime, atime, etime, time;
  ttime=MPI_Wtime();
  etime = 0.0;
  atime = 0.0;

  static int divB_ini_flag = 1;

  const int it_max      = (int) Physics.divB[i_divB_itmax];
  const double divB_err = Physics.divB[i_divB_err];

  static double divB_coeff = 1;
  static double divB_damp  = 0.65;
  /********************************************************************/
  static int next_iter = 3;

  register int i,j,k,iter;

  const int nvar = Physics.NVAR;
  const int i_beg = Grid.lbeg[0];
  const int i_end = Grid.lend[0];
  const int kbeg = Grid.lbeg[2];
  const int kend = Grid.lend[2];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];

  int blcksz;

  const int need_diagnostics = Run.need_diagnostics;

  int next[3];

  double dxmin,dxmax,err_unit,err,conv;
  double err_loc[2],err_glo[2];

  static double divB_fac;

  double invdt=1.0/Run.dt;
    
  blcksz=1;
  err_unit=0.0;
  for(int v=0;v<3;v++){
    next[v] = Grid.stride[v];
    blcksz *= (Grid.lsize[v]+2*Grid.ghosts[v]);
    err_unit += 1.0/Grid.dx[v];
  }

  err_unit = 0.5/err_unit*3.5449;
 
#pragma acc enter data copyin(next[:3])
 
  dxmin = min(Grid.dx[0],Grid.dx[1]);
  dxmax = max(Grid.dx[0],Grid.dx[1]);
  if (Grid.NDIM == 1){
    dxmin = Grid.dx[0];
    dxmax = Grid.dx[0];
  }
  if (Grid.NDIM == 3){
    dxmin = min(dxmin,Grid.dx[2]);
    dxmax = max(dxmin,Grid.dx[2]);
  }
  
  double* divB = Grid.divB;
  double* phi  = Grid.phi;
  double* U    = (double*) Grid.U;
  const int bufsize = Grid.bufsize;
  const int NDIM = Grid.NDIM;

  if (divB_ini_flag == 1){    
    if (Run.rank == 0) { 
      cout << "divB_ini: divB_coeff = "  << divB_coeff  << endl;
      cout << "divB_ini: divB_damp  = "  << divB_damp   << endl; 
      cout << "divB_ini: it_max     = "  << it_max      << endl;
      cout << "divB_ini: divB_err   = "  << divB_err    << endl;
    }
    divB_ini_flag =0;
  }

  divB_fac = divB_coeff*min(dxmin,0.4*dxmax)*dxmin;

  /********************************************************************/
  
  time=MPI_Wtime();
  exchange_B_acc(Grid);    
  etime+=MPI_Wtime()-time;

  err_loc[0] = 0.0;
  err = 0.0;

#pragma acc parallel loop gang collapse(2) \
 present(Grid[:1],U[:bufsize*8], divB[:bufsize],next[:3],Grid.w1[:3],Grid.w2[:3]) \
 reduction(max:err)
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++) {
    int off0 = j*next[1]+k*next[2];
    
#pragma acc loop vector
    for(i=i_beg;i<=i_end;i++)
      divB[off0+i]  = 0.0;
    
#pragma acc loop seq
    for(int v=0;v<NDIM;v++){
      double ww1 = Grid.w1[v];
      double ww2 = Grid.w2[v];
      int next_d = next[v];
      int offBv = off0*nvar+5+v;
      
      int offp1 = offBv + nvar*next_d;
      int offp2 = offBv + nvar*2*next_d;
      int offn1 = offBv - nvar*next_d;
      int offn2 = offBv - nvar*2*next_d;
      
      #pragma ivdep
#pragma acc loop vector
      for(i=i_beg;i<=i_end;i++){
	divB[off0+i] +=
	  (ww1*(U[offp1+i*nvar]-U[offn1+i*nvar])+
	   ww2*(U[offp2+i*nvar]-U[offn2+i*nvar]));
	
      }
    }

    //err = 0.0;  
    #pragma ivdep 
#pragma acc loop vector reduction(max:err)
    for(i=i_beg;i<=i_end;i++)
      err = max(err,fabs(divB[off0+i]));

    //err_loc[0] = max(err_loc[0],err);
    
  } // end 1st kern

  err_loc[0] = err;
  
  PGI_COMPARE(divB, double, Grid.bufsize, "divB", "clean_divB.C",
              "Clean_divB", 1)
  
  for(iter=0;iter<next_iter;iter++){

#pragma acc parallel loop  gang collapse(2) \
 present(divB[:bufsize], phi[:bufsize])
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
      int off0 = j*next[1]+k*next[2];
      
      #pragma ivdep
#pragma acc loop vector
      for(i=i_beg;i<=i_end;i++)
	phi[off0+i] += divB_fac*divB[off0+i] - divB_damp*phi[off0+i];
    }

    PGI_COMPARE(phi, double, Grid.bufsize, "phi", "clean_divB.C",
                "Clean_divB", 2)

    time=MPI_Wtime();    
    exchange_single_acc(Grid,phi);
    etime+=MPI_Wtime()-time;

    PGI_COMPARE(phi, double, Grid.bufsize, "phi", "clean_divB.C",
                "Clean_divB", 3)
    
    /* vertical boundaries divB anti-symmetric */ 
    // This can probably be done differently now that vertical is x

    if( Grid.is_gbeg[0] ) {
      i = Grid.lbeg[0];
#pragma acc parallel loop collapse(2) \
 present(Grid[:1], phi[:bufsize])
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
	phi[Grid.node(i-1,j,k)] = -phi[Grid.node(i,j,k)]; 
	phi[Grid.node(i-2,j,k)] = -phi[Grid.node(i+1,j,k)];      
      }

      PGI_COMPARE(phi, double, Grid.bufsize, "phi", "clean_divB.C",
                  "Clean_divB", 4)

    }
   
 
    if( Grid.is_gend[0] ) {
      i = Grid.lend[0];
#pragma acc parallel loop collapse(2) \
 present(Grid[:1], phi[:bufsize])
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
	phi[Grid.node(i+1,j,k)] = -phi[Grid.node(i,j,k)]; 
	phi[Grid.node(i+2,j,k)] = -phi[Grid.node(i-1,j,k)];   
      }
      PGI_COMPARE(phi, double, Grid.bufsize, "phi", "clean_divB.C",
                  "Clean_divB", 5)
    }

#pragma acc parallel loop collapse(2) \
 present(Grid[:1], U[:bufsize*8], phi[:bufsize], divB[:bufsize], Grid.w1[:3],Grid.w2[:3], \
  Grid.Qres[:bufsize])
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
      int off0 = j*next[1]+k*next[2];

#pragma acc loop seq
      for(int v=0;v<NDIM;v++){
      double ww1 = Grid.w1[v];
      double ww2 = Grid.w2[v];
      int next_d = next[v];
	  int offBv = off0*nvar+5+v;
      
	  int offp1 = off0 + next_d;
	  int offp2 = off0 + 2*next_d;
	  int offn1 = off0 - next_d;
	  int offn2 = off0 - 2*next_d;
	
        #pragma ivdep
#pragma acc loop vector
	for(i=i_beg;i<=i_end;i++){
	  
	  U[offBv+i*nvar] +=
	    (ww1*(phi[offp1+i]-phi[offn1+i])+
	     ww2*(phi[offp2+i]-phi[offn2+i]));

	}
      }

      #pragma ivdep
#pragma acc loop vector
      for(i=i_beg;i<=i_end;i++)
	U[(off0+i)*nvar+4] += phi[off0+i]*divB[off0+i];

      if(need_diagnostics){
        #pragma ivdep
#pragma acc loop vector
	for(i=i_beg;i<=i_end;i++){	  
	  Grid.Qres[off0+i] += phi[off0+i]*divB[off0+i]*invdt;	  	
	}
      }
            
    }

    ACC_COMPARE(U, real, bufsize*8)
    PGI_COMPARE(U, double, Grid.bufsize*8, "U", "clean_divB.C",
                "Clean_divB", 6)
    if(need_diagnostics) {
      PGI_COMPARE(Grid.Qres, double, Grid.bufsize, "Qres",
                  "clean_divB.C", "Clean_divB", 7)
    }

    time=MPI_Wtime();
    exchange_B_acc(Grid);    
    etime+=MPI_Wtime()-time;

    PGI_COMPARE(U, double, Grid.bufsize*8, "U", "clean_divB.C",
                "Clean_divB", 8)
    
    err = 0.0;
#pragma acc parallel loop gang collapse(2) \
 present(Grid[:1],U[:bufsize*8], divB[:bufsize],Grid.w1[:3],Grid.w2[:3]) \
 reduction(max:err)
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
      int off0 = j*next[1]+k*next[2];
    
#pragma acc loop vector
      for(i=i_beg;i<=i_end;i++)
	divB[off0+i]  = 0.0;
    
#pragma acc loop seq
      for(int v=0;v<NDIM;v++){
      double ww1 = Grid.w1[v];
      double ww2 = Grid.w2[v];
      int next_d = next[v];
      
      int offBv = off0*nvar+5+v;
      
	  int offp1 = offBv + nvar*next_d;
      int offp2 = offBv + nvar*2*next_d;
      int offn1 = offBv - nvar*next_d;
      int offn2 = offBv - nvar*2*next_d;

        #pragma ivdep
#pragma acc loop vector
	for(i=i_beg;i<=i_end;i++){
	  divB[off0+i] +=
	    (ww1*(U[offp1+i*nvar]-U[offn1+i*nvar])+
	     ww2*(U[offp2+i*nvar]-U[offn2+i*nvar]));
	  
	}
      }

      #pragma ivdep
#pragma acc loop vector reduction(max:err)
      for(i=i_beg;i<=i_end;i++)
	err = max(err,fabs(divB[off0+i]));
    }

    err_loc[1] = err;

    PGI_COMPARE(divB, double, Grid.bufsize, "divB", "clean_divB.C",
                "Clean_divB", 9)

  }
#pragma acc exit data delete(next[:3])
  time=MPI_Wtime();      
  MPI_Allreduce(err_loc,err_glo,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  atime+=MPI_Wtime()-time;
  
  err_glo[0] *= err_unit;
  err_glo[1] *= err_unit;

  err = err_glo[1]/err_glo[0];

  conv = pow(err_glo[0]/err_glo[1],1./float(iter));

  if(err > divB_err)
    next_iter += 1;
  else
    next_iter -= 1;

  next_iter = min(next_iter,it_max);
  next_iter = max(1,next_iter);
  
  if ( (Run.rank == 0) && (Run.verbose  > 0) ){
    cout << "divB_max [" << iter << "] " << err_glo[0] << " -> " << err_glo[1] << " [G]   ("  << conv << ")"  << endl;
  }

  ttime = MPI_Wtime() - ttime;

  if (Run.rank == 0){
    if(Run.verbose  > 2){
      cout << "divB time : " << ttime << ' ' << etime << ' ' << atime << ' '
           << (etime+atime)/ttime << endl;
    }
  }

}
