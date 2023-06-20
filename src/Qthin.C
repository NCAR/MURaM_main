#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "run.H"
#include "grid.H"
#include "comm_split.H"
#include <iostream>
#include <algorithm>

#define YZ_LOOP(G,j,k) \
for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
for((j)=(G).lbeg[1];(j)<=(G).lend[1];(j)++)

using std::min;
using std::max;
using std::cout;
using std::endl;

// N are the table sizes
// the axes temperature (tax), column depth (cmassax), hydrogen column depth (tauhax)
// And the tables,
// Neutral fractions (_nf), excitaiton (_exc) and Escape probabilities (_esc)
// for Hydrogen, Calcium and Magnesium
int N_T;
double *tax;
double *C_thin;

double dT_Tab;

// Abundances relative to hydrgoen from Asplund 2009
const double mp   = 1.67262158E-24;

// Electron and hydrogen number density from simple H/He mix.
const double XX = 0.7;
const double nh = XX/mp; 

/*****************************************************************************/
void Get_Radloss(const RunData&  Run, GridData& Grid,const PhysicsData& Physics){

  static int radloss_ini_flag = 1;
  
  register int i,j,k,off,node,n,n1,n2;

  const int i_beg    = Grid.lbeg[0];
  const int i_end    = Grid.lend[0];
/*
  int next[3];
  for(i=0;i<3;i++)
    next[i] = Grid.stride[i];
*/

  //int next0 = Grid.stride[0];
  int next1 = Grid.stride[1];
  int next2 = Grid.stride[2];

  double t1,t2,r1,r2,ff,pt,pr,ts,rr,qq,qloss,t_a,t_b,tmin,tmax;
    
  // Chianti assumes Q=-n_e*n_H*L(T)
  const double X_H    = 0.7;
  const double ne_ref = 1e9;
  const double ne_par = sqrt(0.5*(1.0+X_H)*X_H)*6e23/ne_ref;

  const double inv_pmax = 1.0/Physics.rt[i_rt_pre_cut];
  
  static int ntab=0;

  static double* T_tab;
  static double* Q_tab;
  
  static double lgT0,del0;
 
  if(radloss_ini_flag){
    std::ifstream fptr("Radloss_Chianti.dat",std::ios::in);
    if(fptr){
      fptr.precision(16);
      fptr >> ntab;
      T_tab = new double[ntab];
      Q_tab = new double[ntab];
      for(i=0;i<ntab;i++)
	fptr >> T_tab[i] >> Q_tab[i];
      fptr.close();
    } else {
      cout << "Radloss Table not found .... abort" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    lgT0=log(T_tab[0]);
    del0=1.0/(log(T_tab[1])-log(T_tab[0]));

    for(i=0;i<ntab;i++){
      T_tab[i]  = log(T_tab[i]);
      Q_tab[i] *= ne_ref*ne_ref;
     }

    if(Run.rank == 0){
      cout << "Radloss: use log(T), log(rho)" << endl;
      cout << "Radloss: use X_H = " << X_H << endl;
      cout << "Radloss: pressure_cutoff = " << 1./inv_pmax << endl;
    }
   
    radloss_ini_flag = 0;

#pragma acc enter data copyin(T_tab[:ntab])
#pragma acc enter data copyin(Q_tab[:ntab])
  }

  const int ibeg = Grid.lbeg[0];
  const int iend = Grid.lend[0];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];
  const int kbeg = Grid.lbeg[2];
  const int kend = Grid.lend[2];
  const int bufsize = Grid.bufsize;
  double qthin[2][i_end+2];

#pragma acc parallel loop collapse(2) gang \
 present(Grid[:1], Grid.Qthin[:bufsize], Grid.temp[:bufsize], \
         Grid.pres[:bufsize], Grid.U[:bufsize], \
         T_tab[:ntab], Q_tab[:ntab]) \
 private(off, qthin[:2][:i_end+2])
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++) {
//#pragma acc cache(qthin[:2][:i_end+2])
    off = j*next1+k*next2;

#pragma acc loop vector
    for(i=i_beg-1;i<=i_end+1;i++) {
      qthin[0][i] = 0.0;
      qthin[1][i] = 0.0;
      //Grid.Qthin[off+i]=0.0;
    }

    //node = off+i_beg-1;
    //t1 = log(Grid.temp[node]);
    //r1 = log(Grid.U[node].d*ne_par);
    
#pragma acc loop vector private(node, t1, t2, r1, r2, n1, n2, \
 t_a, t_b, ff, ts, pr, pt, rr, qq, qloss)
    for(i=i_beg;i<=i_end+1;i++){
      node = off+i;
     
      t1 = log(Grid.temp[node-1]);
      r1 = log(Grid.U[node-1].d*ne_par);

      t2 = log(Grid.temp[node]);
      r2 = log(Grid.U[node].d*ne_par);

      tmin = min(t1,t2);
      tmax = max(t1,t2);

      n1 = (int) ( (tmin-lgT0)*del0 );
      n2 = (int) ( (tmax-lgT0)*del0 );

      n1 = max(0,n1);
      n2 = min(ntab-2,n2);
      
#pragma acc loop seq
      for(n=n1;n<=n2;n++){
	
	t_a=max(tmin,T_tab[n]);
	t_b=min(tmax,T_tab[n+1]);
	
	if(t_b-t_a > 1e-6 ){
	  ff = (t_b-t_a)/(tmax-tmin);
	  ts = 0.5*(t_a+t_b);
	  pr = (t2-ts)/(t2-t1);
	  pt = (T_tab[n+1]-ts)/(T_tab[n+1]-T_tab[n]);
	} else if (t_b-t_a >= 0.) {
	  ff = 1.0;
	  ts = 0.5*(t_a+t_b);
	  pr = 0.5;
	  pt = (T_tab[n+1]-ts)/(T_tab[n+1]-T_tab[n]);
	} else {
	  ff=0.0;
	  pr=0.5;
	  pt=0.5;
	}
	
	rr = exp(pr*r1+(1.0-pr)*r2);
	qq = pt*Q_tab[n]+(1.0-pt)*Q_tab[n+1];
	
	qloss = -rr*rr*qq*ff;

	//Grid.Qthin[node-1] += qloss*pr;
	//Grid.Qthin[node]   += qloss*(1.0-pr);
	qthin[0][i-1] += qloss*pr;
        qthin[1][i] += qloss*(1.0-pr);
      }
      //t1 = t2;
      //r1 = r2;
    }

#pragma acc loop vector
    for(i=i_beg-1; i<=i_end+1; i++)
      Grid.Qthin[off+i] = qthin[0][i] + qthin[1][i];

    // remove radiative loss in high pressure regions 
#pragma acc loop vector
    for(i=i_beg;i<=i_end;i++)
      Grid.Qthin[off+i] *= max(0.0,1.-pow(Grid.pres[off+i]*inv_pmax,2));
  }

}
