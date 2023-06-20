#include <mpi.h>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <algorithm>
#include "grid.H"
#include "physics.H"
#include "precision.h"
#include "eos.H"
#include "run.H"
#include "muramacc.H"
#include <stdio.h>
#include <iostream>

using std::cout;
using std::endl;

using std::min;
using std::max;

double eos_gamma = 1.65;
double eos_mu = 0.62;
double eos_mu_n = 1.3;
double mh = 1.67262158e-24;
double kb = 1.380622e-16; 

inline int invalid_eos(double a, double b){
  int status = 0;

  if( (a != a) or (a <= 0.) or (a == std::numeric_limits<double>::infinity()) )
    status = 1;

  if( (b != b) or (b <= 0.) or (b == std::numeric_limits<double>::infinity()) )
    status = 1;

  return status;
}

void ConsToPrim(GridData& Grid, const PhysicsData& Physics, const RunData& Run) {

  //double time,s_time;
  //static double t_time = 0.0 ,c_time = 0.0 , r_time = 0.0;
  //static int call_count = 0;
  //s_time = MPI_Wtime();
  
  int i,j,k,node,i_d,i_e;
  int kbeg, kend, jbeg, jend, ibeg, iend;
  const int bufsize = Grid.bufsize;

  double ep,lr,dn,vv,xx;

  const double eps_max=exp(xeps[N_eps-1])-eps_off;

  int sz = Grid.lsize[0]+2*Grid.ghosts[0];
  const int v_nvar = Grid.v_nvar;
  /* ideal gas law for out of bound values */
  const double Rgas=8.314e7;

  const double c_temp= (eos_gamma-1.)*eos_mu/Rgas;
  const double c_pres= (eos_gamma-1.);  

  cState* U = Grid.U;

  double var[5][sz], 
  eps[sz], 
  lgr[sz], pres[sz], temp[sz],nel[sz],amb[sz],rhoi[sz],coeff0,coeff1,coeff2,coeff3;

  kbeg = Grid.lbeg[2]-Grid.ghosts[2];
  kend = Grid.lend[2]+Grid.ghosts[2];
  jbeg = Grid.lbeg[1]-Grid.ghosts[1];
  jend = Grid.lend[1]+Grid.ghosts[1];
  ibeg = Grid.lbeg[0]-Grid.ghosts[0];
  iend = Grid.lend[0]+Grid.ghosts[0];

#pragma acc enter data copyin(Grid[:1], U[:bufsize], Grid.pres[:bufsize],           \
  Grid.temp[:bufsize], Grid.ne[:bufsize],                      \
  Grid.rhoi[:bufsize], Grid.amb[:bufsize],                     \
  xeps[:N_eps], xlr[:N_lr], p_eostab[:N_eps][:N_lr],           \
  T_eostab[:N_eps][:N_lr], s_eostab[:N_eps][:N_lr],            \
  ne_eostab[:N_eps][:N_lr], rhoi_eostab[:N_eps][:N_lr],        \
  amb_eostab[:N_eps][:N_lr])

//cout << "Before loop" << endl;
#pragma acc parallel loop collapse(2) gang                     \
 default(present)                                   \
 private(var[:5][:sz],         \
   eps[:sz],                           \
  lgr[:sz], pres[:sz], temp[:sz], nel[:sz],                    \
  amb[:sz], rhoi[:sz])
  for(k=kbeg; k<=kend; k++){
    for(j=jbeg; j<=jend; j++){
      //time = MPI_Wtime();
#pragma acc loop vector private(node)
      for(i=0;i<sz;i++){
        node = Grid.node(i,j,k);
        var[0][i] = U[node].d;
        var[1][i] = U[node].M.x;
        var[2][i] = U[node].M.y;
        var[3][i] = U[node].M.z;
        var[4][i] = U[node].e;
        //var[5][i] = U[node].B.x;
        //var[6][i] = U[node].B.y;
        //var[7][i] = U[node].B.z;
      }
      //r_time += MPI_Wtime()-time;

      //time = MPI_Wtime();
#pragma acc loop vector private(dn, vv) 
      for(i=0;i<sz;i++){
        dn        = 1.0/var[0][i];
        var[1][i] = var[1][i]*dn;
        var[2][i] = var[2][i]*dn;
        var[3][i] = var[3][i]*dn;

        vv  = (var[1][i]*var[1][i]+
        var[2][i]*var[2][i]+
	var[3][i]*var[3][i]);
        var[4][i] = var[4][i] - 0.5*vv*var[0][i];
  
        eps[i] = var[4][i]*dn +eps_off;

        lgr[i] = log(var[0][i]);
      }

#pragma acc loop vector private(ep, lr, i_d, i_e,coeff0,coeff1,coeff2,coeff3)
      for(i=0;i<sz;i++){

        ep = log(eps[i]);
        lr = lgr[i];

	i_d = (int) ( (lr-xlr[0])*inv_del_lr );
        i_d = max(0,i_d);
        i_d = min(N_lr-2,i_d);
        
        i_e = (int) ((ep-xeps[0])*inv_del_eps );
        i_e = max(0,i_e);
        i_e = min(N_eps-2,i_e);
  
        coeff0 = (ep-xeps[i_e])   * (lr-xlr[i_d]);
        coeff1 = (ep-xeps[i_e])   * (xlr[i_d+1]-lr);
        coeff2 = (xeps[i_e+1]-ep) * (lr-xlr[i_d]);
        coeff3 = (xeps[i_e+1]-ep) * (xlr[i_d+1]-lr);
          
        double pp  = coeff0 * (double) p_eostab[i_e+1][i_d+1];
        pp += coeff1 * (double) p_eostab[i_e+1][i_d];
        pp += coeff2 * (double) p_eostab[i_e][i_d+1];
        pp += coeff3 * (double) p_eostab[i_e][i_d];
        
        double tt  = coeff0 * (double) T_eostab[i_e+1][i_d+1];
        tt += coeff1 * (double) T_eostab[i_e+1][i_d];
        tt += coeff2 * (double) T_eostab[i_e][i_d+1];
        tt += coeff3 * (double) T_eostab[i_e][i_d];

        double nn  = coeff0 * (double) ne_eostab[i_e+1][i_d+1];
        nn += coeff1 * (double) ne_eostab[i_e+1][i_d];
        nn += coeff2 * (double) ne_eostab[i_e][i_d+1];
        nn += coeff3 * (double) ne_eostab[i_e][i_d];
       
        double ii = coeff0 * (double) rhoi_eostab[i_e+1][i_d+1]; 
        ii += coeff1 * (double) rhoi_eostab[i_e+1][i_d]; 
        ii += coeff2 * (double) rhoi_eostab[i_e][i_d+1]; 
        ii += coeff3 * (double) rhoi_eostab[i_e][i_d];

        double aa  = coeff0 * (double) amb_eostab[i_e+1][i_d+1];
        aa += coeff1 * (double) amb_eostab[i_e+1][i_d];
        aa += coeff2 * (double) amb_eostab[i_e][i_d+1];
        aa += coeff3 * (double) amb_eostab[i_e][i_d];

        pres[i] = exp(pp);
        temp[i] = exp(tt);
        nel[i] = exp(nn);
        rhoi[i] = exp(ii);
        amb[i] = exp(aa);
      }
      
#pragma acc loop vector 
      for(i=0;i<sz;i++){
    
      // Out of upper table bounds
      if( (lgr[i] < xlr[0]) ) {
        pres[i]=c_pres*var[4][i];
        temp[i]=c_temp*eps[i];
        // Use PV = nrt for ne
        nel[i]= pres[i]/(temp[i]*kb)-var[0][i]/(eos_mu_n*mh);
      }
  
      if (eps[i] > 0.8*eps_max) {
        xx = max(0.0,5.0*(eps_max-eps[i])/eps_max);
        pres[i]=pres[i]*xx+(1.0-xx)*c_pres*var[4][i];
        temp[i]=temp[i]*xx+(1.0-xx)*c_temp*eps[i];
        // Use PV = nrt for ne
        nel[i]= pres[i]/(temp[i]*kb)-var[0][i]/(eos_mu_n*mh);
      }
  
      if(invalid_eos(pres[i],temp[i])){
        pres[i]=c_pres*var[4][i];
        temp[i]=c_temp*eps[i];
      }

      // if nel is inf use pv=nrt
      if(nel[i] == std::numeric_limits<double>::infinity()) nel[i]= pres[i]/(temp[i]*kb)-var[0][i]/(eos_mu_n*mh);

      // If electron number goes too small, or NaNs, or is still inf set to a minimum
      if ((nel[i] < 1.0) || (nel[i]!=nel[i]) || (nel[i] == std::numeric_limits<double>::infinity())) nel[i] = 1.0;

      } // end loop

      //time = MPI_Wtime();
#pragma acc loop vector private(node) 
      for(i=ibeg;i<=iend;i++){
        node = Grid.node(i,j,k);
        U[node].M.x=var[1][i];
        U[node].M.y=var[2][i];
        U[node].M.z=var[3][i];
        U[node].e  =var[4][i];
        
        Grid.pres[node] = pres[i];
        Grid.temp[node] = temp[i];
        Grid.ne[node] = nel[i];
        Grid.rhoi[node] = rhoi[i];
        Grid.amb[node] = amb[i];

      }
      //r_time += MPI_Wtime()-time;

    }
  }

  PGI_COMPARE(Grid.U, double, Grid.bufsize*8, "U", "cons_to_prim.C", "ConsToPrim", 1)
  PGI_COMPARE(Grid.pres, double, Grid.bufsize, "pres", "cons_to_prim.C", "ConsToPrim", 2)
  PGI_COMPARE(Grid.temp, double, Grid.bufsize, "temp", "cons_to_prim.C", "ConsToPrim", 3)
  PGI_COMPARE(Grid.ne, double, Grid.bufsize, "ne", "cons_to_prim.C", "ConsToPrim", 4)
  PGI_COMPARE(Grid.rhoi, double, Grid.bufsize, "rhoi", "cons_to_prim.C", "ConsToPrim", 5)
  PGI_COMPARE(Grid.amb, double, Grid.bufsize, "amb", "cons_to_prim.C", "ConsToPrim", 6)

}

