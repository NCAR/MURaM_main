#include <mpi.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "muramacc.H"
#include <iostream>

using std::cout;
using std::endl;

double inv_va2max,v_lim=0.0,e_lim=0.0;;

double lfac(double va2) {
  double x2,x4,s;
  x2=va2*inv_va2max;
  x4=x2*x2;
  s=1.0+x2;
  return s/(s+x4);
}

void Adjust_Valf_Max(const RunData& Run,const GridData& Grid,
		     const PhysicsData& Physics) {
  
  static int ini_flag = 1;

  const double max_fill  = Physics.params[i_param_max_fill];
  const double cs_fac    = 2.0;
  const double vv_fac    = 3.0;
  const double eps_decay = 0.99;
  
  register int i,j,k,node;
  
  double l_max[3], g_max[3], l_sum[2], g_sum[2], va_max, vv, ee;

  static double va_max_old;
  
  const int ibeg = Grid.lbeg[0];
  const int iend = Grid.lend[0];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];
  const int kbeg = Grid.lbeg[2];
  const int kend = Grid.lend[2];
  const int bufsize = Grid.bufsize;

  double l_max0 = 0.0;
  double l_max1 = 0.0;
  double l_max2 = 0.0;
  double l_sum0 = 0.0;
  double l_sum1 = 0.0;

//#pragma acc update device(Grid.U[:bufsize])
//#pragma acc update device(Grid.pres[:bufsize])

#pragma acc parallel loop collapse(3)                              \
 private(node, vv, ee)                                             \
 reduction(max:l_max0) reduction(max:l_max1) reduction(max:l_max2) \
 reduction(+:l_sum0) reduction(+:l_sum1)                           \
 present(Grid[:1], Grid.U[:bufsize], Grid.pres[:bufsize])
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++)
  for(i=ibeg; i<=iend; i++) {
    node = Grid.node(i,j,k);
    vv   = Grid.U[node].M.abs();
    ee   = Grid.U[node].e/Grid.U[node].d;
    
    l_max0 = fmax(l_max0,vv);
    l_max1 = fmax(l_max1,5/3.*Grid.pres[node]/Grid.U[node].d);
    l_max2 = fmax(l_max2,ee);
    
    if(vv >= v_lim*0.95){
      l_sum0 += 1.0;
    }

     if(ee >= e_lim*0.95){
      l_sum1 += 1.0;
    }
  }

  l_max[0] = l_max0;
  l_max[1] = l_max1;
  l_max[2] = l_max2;
  l_sum[0] = l_sum0;
  l_sum[1] = l_sum1;

  PGI_COMPARE(l_max, double, 3, "l_max", "Adjust_Valf_Max.C", "Adjust_Valf_Max", 1)
  PGI_COMPARE(l_sum, double, 2, "l_sum", "Adjust_Valf_Max.C", "Adjust_Valf_Max", 2)
  
  MPI_Allreduce(l_max,g_max,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  g_max[1] = sqrt(g_max[1]);
  
  MPI_Allreduce(l_sum,g_sum,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(Physics.params[i_param_va_adjust] > 0.0){
    // first call set v_lim to v_max, after that adjust    
    if(ini_flag){
      v_lim = fmin(g_max[0],Physics.tchk[i_tchk_vmax]);
      e_lim = fmin(g_max[2],Physics.tchk[i_tchk_eps_max]);

      va_max_old = 0.0;
      
      if(Run.rank == 0)
	cout << "Dynamically limit vmax & emax! max_fill = " << max_fill << endl;
      
      ini_flag = 0;
    } else {
      if(g_sum[0] >= 4*max_fill)
	v_lim *= 1.04;
      else if(g_sum[0] >= max_fill)
	v_lim *= 1.01;
      else if(g_sum[0] < 0.5*max_fill)
	v_lim *= 0.99;

      if(g_sum[1] >= 4*max_fill)
	e_lim *= 1.04;
      else if(g_sum[1] >= max_fill)
	e_lim *= 1.01;
      else if(g_sum[1] < 0.5*max_fill)
	e_lim *= 0.99;
    }
    
    va_max = fmax(Physics.params[i_param_va_max],
	Physics.params[i_param_va_adjust]*vv_fac*g_max[0]);
    va_max = fmax(va_max,Physics.params[i_param_va_adjust]*cs_fac*g_max[1]);

    va_max = fmax(va_max,va_max_old*eps_decay);
    
    v_lim = fmax(v_lim,Physics.params[i_param_va_max]/vv_fac);
    v_lim = fmin(v_lim,Physics.tchk[i_tchk_vmax]);

    e_lim = fmax(e_lim,pow(Physics.params[i_param_va_max]/cs_fac,2));
    e_lim = fmin(e_lim,Physics.tchk[i_tchk_eps_max]);
  }  else {
    va_max = Physics.params[i_param_va_max];
    v_lim  = Physics.tchk[i_tchk_vmax];
    e_lim  = Physics.tchk[i_tchk_eps_max];
  }

  inv_va2max = 1.0/pow(va_max,2);

  va_max_old = va_max;
  
  if( (Run.rank == 0) && (Run.verbose >0) )
    std::cout << "Adjust_va_max [" << Run.globiter<<"]   " << g_max[0] << "   " << g_max[1] << "   " << va_max << "   |   "
	 << v_lim << " (" << g_sum[0]<< ")   " << e_lim << " (" << g_sum[1] << ") " << std::endl;
  
}
