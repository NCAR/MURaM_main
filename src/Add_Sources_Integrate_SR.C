#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "src_int_tck.H"
#include "limit_va.H"
#include "exchange.H"
#include <stdio.h>
#include "ACCH.h"

#define XZ_LOOP(G,i,k) \
  for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
  for((i)=(G).lbeg[0];(i)<=(G).lend[0];(i)++)

#define YZ_LOOP(G,j,k) \
  for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
  for((j)=(G).lbeg[1];(j)<=(G).lend[1];(j)++)

#define YZ_FULL_LOOP(G,j,k) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lend[1]+(G).ghosts[1];(j)++)

#define GLOOP(gs,i,j,d1,d2) \
  for((j)=(gs)[(d2)][0];(j)<=(gs)[(d2)][1];(j)++) \
  for((i)=(gs)[(d1)][0];(i)<=(gs)[(d1)][1];(i)++)

using std::cout;
using std::endl;
using std::min;
using std::max;

void get_damping(const RunData&, const GridData&, const PhysicsData&, double*);
/*****************************************************************************/
void Source_Integrate_Tcheck(const RunData& Run, GridData& Grid, 
			     const PhysicsData& Physics, const int stage) {
  
  int i,j,k,node,off0,v;
  static int ini_flag=1;

  const double wdt = Run.dt/double(maxstage+1-stage);
  const int damp_on = Physics.dmp[i_dmp_switch];

  const int vsize = Grid.vsize;
  const int bufsize = Grid.bufsize;
  const int i_beg = Grid.lbeg[0];
  const int i_end = Grid.lend[0];
  const int kbeg  = Grid.lbeg[2];
  const int kend  = Grid.lend[2];
  const int jbeg  = Grid.lbeg[1];
  const int jend  = Grid.lend[1];
  
  int next1,next2;
  next1 = Grid.stride[1];
  next2 = Grid.stride[2];
   
  const double eps_min = Physics.tchk[i_tchk_eps_min];
  const double rho_min = Physics.tchk[i_tchk_rho_min];

  const bool ambipolar = (Physics.params[i_param_ambipolar] > 0.0);
  const bool spitzer   = (Physics.params[i_param_spitzer] > 0.0);

  double vmax,vv,ekin,eps_max;
  double ee_low, ee_up,ee,efaci,e_dmpi, Qrad[vsize];

  if(v_lim > 0.0)
    vmax = v_lim;
  else
    vmax = Physics.tchk[i_tchk_vmax];

  if(e_lim > 0.0)
    eps_max = e_lim;
  else
    eps_max = Physics.tchk[i_tchk_eps_max];

  const double g = Physics.g[0];

  double dn,fsr_x,fsr_y,fsr_z,va2,bb,bf,s,x2,x4;
  double rh,vx,vy,vz,bx,by,bz,rd;

  const double invdt  = 1.0/Run.dt; 
  const double invRdt = Physics.rt[i_rt_cfl]/Run.dt;
  
  const double ambvel_max = Physics.params[i_param_ambvel_max];
  const int need_diagnostics = Run.need_diagnostics;
  const int ext_cor = Physics.rt_ext[i_ext_cor];
//#pragma acc enter data copyin(need_diagnostics)
  double boris0,boris1,boris2,boris3;

  int sz = Grid.lsize[0]+2*Grid.ghosts[0];
  static double* dmp;

  if(ini_flag)
  {
    dmp = (double*) ACCH::Malloc(sz*sizeof(double));
    memset(dmp,0,sz*sizeof(double));
#pragma acc update device(dmp[:sz])
    ini_flag=0;
  }
  //static double* dmp =  (double*) ACCH::Malloc(sz*sizeof(double)); //malloc(sz*sizeof(double));
//#pragma acc enter data create(dmp[:sz])

  if((stage == 1)&&(damp_on)) {
    get_damping(Run,Grid,Physics,dmp);
//#pragma acc update device(dmp[:sz])
  }
  
#pragma acc data present(Grid[:1], Grid.U[:bufsize], Grid.U0[:bufsize],          \
  Grid.Res[:bufsize], Grid.Qtot[:bufsize],                       \
  Grid.Qthin[:bufsize],                                          \
  Grid.tvar4[:bufsize], Grid.tvar5[:bufsize],                    \
  Grid.tvar6[:bufsize], Grid.tvar7[:bufsize],                    \
  Grid.sflx[:bufsize], Grid.sflx0[:bufsize],                     \
  Grid.Rflx[:bufsize], Grid.v_amb[:bufsize],                     \
  Grid.v0_amb[:bufsize], Grid.R_amb[:bufsize],                   \
  dmp[:sz])
{

if(need_diagnostics){
#pragma acc parallel loop collapse(2) gang                       \
 present(Grid[:1], Grid.U[:bufsize], Grid.U0[:bufsize],          \
  Grid.Res[:bufsize], Grid.Qtot[:bufsize],                       \
  Grid.Qthin[:bufsize],                                          \
  Grid.tvar4[:bufsize], Grid.tvar5[:bufsize],                    \
  Grid.tvar6[:bufsize], Grid.tvar7[:bufsize],                    \
  Grid.sflx[:bufsize], Grid.sflx0[:bufsize],                     \
  Grid.Rflx[:bufsize], Grid.v_amb[:bufsize],                     \
  Grid.v0_amb[:bufsize], Grid.R_amb[:bufsize],                   \
  dmp[:sz])                                                      \
 private(off0,      \
   Qrad[:vsize])
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++) {
    off0 = j*next1+k*next2;

    if(ext_cor==1){
      #pragma ivdep
#pragma acc loop vector private(node)
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Qrad[i] = Grid.Qtot[node]+Grid.Qthin[node];
      }
    } else {
      #pragma ivdep
#pragma acc loop vector private(node) 
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Qrad[i] = Grid.Qtot[node];
      }
    }

    #pragma ivdep
#pragma acc loop vector private(node, e_dmpi, efaci, dn) 
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      dn      = invRdt*Grid.U[node].e;
      efaci = dn/fmax(dn,fabs(Qrad[i]));

      // gravity + damping + RT
      dn       = g*Grid.U[node].d-dmp[i];
      e_dmpi = dn*Grid.U[node].M.x;
      Grid.Res[node].M.x += dn;
      Grid.Res[node].e   += e_dmpi + Qrad[i]*efaci;
        // diagnostics
        Grid.tvar4[node] = efaci;
        Grid.tvar5[node] = e_dmpi;
     } 

    // Boris correction
    #pragma ivdep
#pragma acc loop vector private(node, rh, rd, vx, vy, vz, \
 bx, by, bz, fsr_x, fsr_y, fsr_z, bb, va2, x2, x4, bf, dn,boris0,boris1,boris2,boris3)
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      rh = Grid.U[node].d;
      rd = Grid.Res[node].d;
      vx = Grid.U[node].M.x;
      vy = Grid.U[node].M.y;
      vz = Grid.U[node].M.z;
      bx = Grid.U[node].B.x;
      by = Grid.U[node].B.y;
      bz = Grid.U[node].B.z;

      fsr_x = Grid.Res[node].M.x;
      fsr_y = Grid.Res[node].M.y;
      fsr_z = Grid.Res[node].M.z;

      bb  = bx*bx+by*by+bz*bz;
      va2 = bb/rh;

      // approximation of 1-1/sqrt(1+(va/c)^4)
      x2 = va2*inv_va2max;
      x4 = x2*x2;
      bf = x4/(1+x2+x4);

      // convert -div(rho*v:v) to -rho*(v*grad)v by adding v*div(rho*v)
      fsr_x -= vx*rd;
      fsr_y -= vy*rd;
      fsr_z -= vz*rd;

      dn = (fsr_x*bx+fsr_y*by+fsr_z*bz)/fmax(1e-100,bb);

      boris0 = -bf*(fsr_x-dn*bx);
      boris1 = -bf*(fsr_y-dn*by);
      boris2 = -bf*(fsr_z-dn*bz);

      boris3 = boris0*vx+boris1*vy+boris2*vz;
      
      Grid.Res[node].M.x += boris0;
      Grid.Res[node].M.y += boris1;
      Grid.Res[node].M.z += boris2;
      Grid.Res[node].e   += boris3;
      Grid.tvar6[node] =  boris3;
    }

    // time integration - MHD
    #pragma ivdep
#pragma acc loop vector private(node) 
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;
      Grid.U[node]   = Grid.U0[node] + wdt*Grid.Res[node];
      Grid.Res[node].d = 0;
      Grid.Res[node].M.x = 0;
      Grid.Res[node].M.y = 0;
      Grid.Res[node].M.z = 0;
      Grid.Res[node].e = 0;
      Grid.Res[node].B.x = 0;
      Grid.Res[node].B.y = 0;
      Grid.Res[node].B.z = 0;
    }

    // time integration - hyperbolic conduction   
    if(spitzer){
      #pragma ivdep
#pragma acc loop vector private(node)
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Grid.sflx[node] = Grid.sflx0[node]+wdt*Grid.Rflx[node];
        Grid.Rflx[node] = 0;
      }
    }

    if(ambipolar){
      //#pragma ivdep
#pragma acc loop vector private(node,vv,s) 
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Grid.v_amb[node] = Grid.v0_amb[node]+wdt*Grid.R_amb[node];
        vv = Grid.v_amb[node].abs();
        s  = ambvel_max/fmax(ambvel_max,vv);
        Grid.v_amb[node] *= s;
      }
#ifdef _OPENACC
     #pragma acc loop vector private(node) 
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Grid.R_amb[node].x = 0;
        Grid.R_amb[node].y = 0;
        Grid.R_amb[node].z = 0;
      }
#else
      memset(&Grid.R_amb[off0+i_beg],0.0,(i_end-i_beg+1)*sizeof(Vector));
#endif
    }

    // Tcheck - catch extreme values before they cause problems
    #pragma ivdep
#pragma acc loop vector private(node, ee_low,ee_up,dn, vv, s, ekin,ee) 
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      dn = fmax(rho_min,Grid.U[node].d);

      vv = Grid.U[node].M.abs()/dn;
      s  = vmax/fmax(vmax,vv);

      ekin = 0.5*dn*vv*vv;
      ee = Grid.U[node].e;
      ee_low = eps_min*dn+ekin*s*s;
      ee_up  = eps_max*dn+ekin*s*s;
      
      Grid.U[node].e  = fmin(fmax(ee,ee_low),ee_up);


      Grid.U[node].d  = dn;
      Grid.U[node].M *= s;
      Grid.tvar7[node] += (Grid.U[node].e-ee)*invdt;
      //Grid.U[node].e  = fmin(fmax(ee[i],ee_low[i]),ee_up[i]);
    }

  }// end loop
 } else {  //does not need diagnostics
	   //
#pragma acc parallel loop collapse(2) gang                       \
 present(Grid[:1], Grid.U[:bufsize], Grid.U0[:bufsize],          \
  Grid.Res[:bufsize], Grid.Qtot[:bufsize],                       \
  Grid.Qthin[:bufsize],                                          \
  Grid.tvar4[:bufsize], Grid.tvar5[:bufsize],                    \
  Grid.tvar6[:bufsize], Grid.tvar7[:bufsize],                    \
  Grid.sflx[:bufsize], Grid.sflx0[:bufsize],                     \
  Grid.Rflx[:bufsize], Grid.v_amb[:bufsize],                     \
  Grid.v0_amb[:bufsize], Grid.R_amb[:bufsize],                   \
  dmp[:sz])                                                      \
 private(off0,      \
  Qrad[:vsize])
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++) {
    off0 = j*next1+k*next2;

    if(ext_cor==1){
      #pragma ivdep
#pragma acc loop vector private(node)
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Qrad[i] = Grid.Qtot[node]+Grid.Qthin[node];
      }
    } else {
      #pragma ivdep
#pragma acc loop vector private(node) 
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Qrad[i] = Grid.Qtot[node];
      }
    }

    #pragma ivdep
#pragma acc loop vector private(node, e_dmpi,efaci, dn) 
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      dn      = invRdt*Grid.U[node].e;
      efaci = dn/fmax(dn,fabs(Qrad[i]));

      // gravity + damping + RT
      dn       = g*Grid.U[node].d-dmp[i];
      e_dmpi = dn*Grid.U[node].M.x;
      Grid.Res[node].M.x += dn;
      Grid.Res[node].e   += e_dmpi + Qrad[i]*efaci;
    }


    // Boris correction
    #pragma ivdep
#pragma acc loop vector private(node, rh, rd, vx, vy, vz, \
 bx, by, bz, fsr_x, fsr_y, fsr_z, bb, va2, x2, x4, bf, dn,boris0,boris1,boris2,boris3)
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      rh = Grid.U[node].d;
      rd = Grid.Res[node].d;
      vx = Grid.U[node].M.x;
      vy = Grid.U[node].M.y;
      vz = Grid.U[node].M.z;
      bx = Grid.U[node].B.x;
      by = Grid.U[node].B.y;
      bz = Grid.U[node].B.z;

      fsr_x = Grid.Res[node].M.x;
      fsr_y = Grid.Res[node].M.y;
      fsr_z = Grid.Res[node].M.z;

      bb  = bx*bx+by*by+bz*bz;
      va2 = bb/rh;

      // approximation of 1-1/sqrt(1+(va/c)^4)
      x2 = va2*inv_va2max;
      x4 = x2*x2;
      bf = x4/(1+x2+x4);

      // convert -div(rho*v:v) to -rho*(v*grad)v by adding v*div(rho*v)
      fsr_x -= vx*rd;
      fsr_y -= vy*rd;
      fsr_z -= vz*rd;

      dn = (fsr_x*bx+fsr_y*by+fsr_z*bz)/fmax(1e-100,bb);

      boris0 = -bf*(fsr_x-dn*bx);
      boris1 = -bf*(fsr_y-dn*by);
      boris2 = -bf*(fsr_z-dn*bz);

      boris3 = boris0*vx+boris1*vy+boris2*vz;
      
      Grid.Res[node].M.x += boris0;
      Grid.Res[node].M.y += boris1;
      Grid.Res[node].M.z += boris2;
      Grid.Res[node].e   += boris3;
    }

    // time integration - MHD
    #pragma ivdep
#pragma acc loop vector private(node) 
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;
      Grid.U[node]   = Grid.U0[node] + wdt*Grid.Res[node];
      Grid.Res[node].d = 0;
      Grid.Res[node].M.x = 0;
      Grid.Res[node].M.y = 0;
      Grid.Res[node].M.z = 0;
      Grid.Res[node].e = 0;
      Grid.Res[node].B.x = 0;
      Grid.Res[node].B.y = 0;
      Grid.Res[node].B.z = 0;
    }

    // time integration - hyperbolic conduction   
    if(spitzer){
      #pragma ivdep
#pragma acc loop vector private(node)
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Grid.sflx[node] = Grid.sflx0[node]+wdt*Grid.Rflx[node];
        Grid.Rflx[node] = 0;
      }
    }

    if(ambipolar){
      //#pragma ivdep
#pragma acc loop vector private(node,vv,s) 
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Grid.v_amb[node] = Grid.v0_amb[node]+wdt*Grid.R_amb[node];
        vv = Grid.v_amb[node].abs();
        s  = ambvel_max/fmax(ambvel_max,vv);
        Grid.v_amb[node] *= s;
      }
#ifdef _OPENACC
     #pragma acc loop vector private(node) 
      for(i=i_beg;i<=i_end;i++){
        node = off0+i;
        Grid.R_amb[node].x = 0;
        Grid.R_amb[node].y = 0;
        Grid.R_amb[node].z = 0;
      }
#else
      memset(&Grid.R_amb[off0+i_beg],0.0,(i_end-i_beg+1)*sizeof(Vector));
#endif
     }

    // Tcheck - catch extreme values before they cause problems
    #pragma ivdep
#pragma acc loop vector private(node, ee,ee_low,ee_up,dn, vv, s, ekin) 
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      dn = fmax(rho_min,Grid.U[node].d);

      vv = Grid.U[node].M.abs()/dn;
      s  = vmax/fmax(vmax,vv);

      ekin = 0.5*dn*vv*vv;
      ee = Grid.U[node].e;
      ee_low = eps_min*dn+ekin*s*s;
      ee_up  = eps_max*dn+ekin*s*s;
      
      Grid.U[node].e  = fmin(fmax(ee,ee_low),ee_up);

      Grid.U[node].d  = dn;
      Grid.U[node].M *= s;
      //Grid.U[node].e  = fmin(fmax(ee[i],ee_low[i]),ee_up[i]);
    }

  } // end loop

} //end if
} //end data

}

/*****************************************************************************/
void SaveCons(GridData& Grid,const PhysicsData& Physics) {

  const int i_beg    = Grid.lbeg[0];
  const int i_end    = Grid.lend[0];
  const int kbeg = Grid.lbeg[2];
  const int kend = Grid.lend[2];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];

  int off0,i,j,k;
  int bufsz = i_end-i_beg+1;
  const int bufsize = Grid.bufsize;

  int next1,next2;
  next1 = Grid.stride[1];
  next2 = Grid.stride[2];

#pragma acc parallel loop collapse(2) gang              \
 present(Grid[:1], Grid.U[:bufsize], Grid.U0[:bufsize]) \
 private(off0) async
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++) {
    off0 = j*next1+k*next2;
#pragma acc loop vector
    for(i=i_beg; i<=i_end; i++) {
      Grid.U0[off0+i] = Grid.U[off0+i];
    }
  }

  if(Physics.params[i_param_spitzer] > 0.0) {
#pragma acc parallel loop collapse(2) gang                    \
 present(Grid[:1], Grid.sflx[:bufsize], Grid.sflx0[:bufsize]) \
 private(off0) async
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
      off0 = j*next1+k*next2;
#pragma acc loop vector
      for(i=i_beg; i<=i_end; i++) {
        Grid.sflx0[off0+i] = Grid.sflx[off0+i];
      }
    }
  }

  if(Physics.params[i_param_ambipolar] > 0.0) {
#pragma acc parallel loop collapse(2) gang                      \
 present(Grid[:1], Grid.v_amb[:bufsize], Grid.v0_amb[:bufsize]) \
 private(off0) async
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
      off0 = j*next1+k*next2;
#pragma acc loop vector
      for(i=i_beg; i<=i_end; i++) {
        Grid.v0_amb[off0+i] = Grid.v_amb[off0+i];
      }
    }
  }
}
 
/*****************************************************************************/
void TCheck(const RunData&  Run, GridData& Grid, const PhysicsData& Physics) {

  static int ini_flag = 1;

  // called after Integrate(); U contains conservative Vars.

  register int i,j,k,v,off0,node;
 
  const int vsize = Grid.vsize;

  const double eps_min = Physics.tchk[i_tchk_eps_min];
  const double rho_min = Physics.tchk[i_tchk_rho_min];
  
  double dn,vmax,vv,ekin,s,eps_max;
  double ee_low, ee_up,eei;

  const int i_beg = Grid.lbeg[0];
  const int i_end = Grid.lend[0];
  const int kbeg = Grid.lbeg[2];
  const int kend = Grid.lend[2];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];

  const int need_diagnostics = Run.need_diagnostics;
  const int bufsize = Grid.bufsize;

  double invdt;
  if(Run.dt > 0.0)
    invdt = 1.0/Run.dt;
  else
    invdt = 0.0;

  int next1,next2;
  next1 = Grid.stride[1];
  next2 = Grid.stride[2];

  if(v_lim > 0.0)
    vmax = v_lim;
  else
    vmax = Physics.tchk[i_tchk_vmax];

  if(e_lim > 0.0)
    eps_max = e_lim;
  else
    eps_max = Physics.tchk[i_tchk_eps_max];
  
  if(ini_flag){
    if (Run.rank == 0) {
      cout << "tcheck: rho_min   = "  << rho_min   << endl;
      cout << "tcheck: eps_min   = "  << eps_min   << endl;
      cout << "tcheck: eps_max   = "  << eps_max   << endl;    
      cout << "tcheck: vmax      = "  << vmax      << endl;
    }
  
    ini_flag = 0;
  }

#pragma acc parallel loop collapse(2) gang             \
 present(Grid, Grid.U[:bufsize], Grid.tvar7[:bufsize]) \
 private(off0)
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++) {
    off0 = j*next1+k*next2;

    #pragma ivdep
#pragma acc loop vector private(node,eei, ee_low,ee_up,dn, vv, s, ekin)
    for(i=i_beg;i<=i_end;i++){
      node = off0+i;

      dn = fmax(rho_min,Grid.U[node].d);

      vv = Grid.U[node].M.abs()/dn;
      s  = vmax/fmax(vmax,vv);

      ekin = 0.5*dn*vv*vv;

      ee_low = eps_min*dn+ekin*s*s;
      ee_up  = eps_max*dn+ekin*s*s;
      eei = Grid.U[node].e;
      Grid.U[node].e  = fmin(fmax(eei,ee_low),ee_up);
      //ee[i] = Grid.U[node].e;

      Grid.U[node].d  = dn;
      Grid.U[node].M *= s;
      if(need_diagnostics){
        Grid.tvar7[node] += (Grid.U[node].e-eei)*invdt;
       }
      //Grid.U[node].e  = fmin(fmax(ee[i],ee_low[i]),ee_up[i]);
    }

  }
}

/*****************************************************************************/
void get_damping(const RunData&  Run, const GridData& Grid,
		 const PhysicsData& Physics, double* my_mean) {

  static int ini_flag = 1;

  int i,j,k,node;

  double hh,hmax,hphot,flxmax,vmax;
  double sbuf[4],rbuf[4];

  const int hsize = Grid.gsize[1]*Grid.gsize[2];
  const int vsize = Grid.lsize[0]+2*Grid.ghosts[0];

  //static double *lbuf    = (double*) ACCH::Malloc(vsize*sizeof(double)); //malloc(vsize*sizeof(double));
  //static double *rh_mean = (double*) ACCH::Malloc(vsize*sizeof(double)); //malloc(vsize*sizeof(double));
  static double *lbuf;
  static double *rh_mean;
  static double *rlbuf;
  static double rho_max;

  const double tau_ref = 1.0/Physics.dmp[i_dmp_tau_ref];
  const double vel_ref = Physics.dmp[i_dmp_vel_ref];
  const double tau_min = 1.0/Physics.dmp[i_dmp_tau_min];

  const int kbeg = Grid.lbeg[2];
  const int kend = Grid.lend[2];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];
  const int ibeg = Grid.lbeg[0];
  const int iend = Grid.lend[0];

  //-------------------------------------------------------------------

  if (ini_flag == 1){
    lbuf = (double*) ACCH::Malloc(vsize*sizeof(double));
    rlbuf = (double*) ACCH::Malloc(vsize*sizeof(double));
    rh_mean = (double*) ACCH::Malloc(vsize*sizeof(double));
    hh = 0.0;
#pragma acc parallel loop present(lbuf[:vsize],rlbuf[:vsize])
    for(i=0;i<vsize;i++) {
        lbuf[i] = 0.0;
        rlbuf[i] = 0.0;
     }

#pragma acc parallel loop collapse(3) \
 reduction(max:hh) \
 present(Grid[:1], Grid.U[:Grid.bufsize], lbuf[:vsize])
    for(k=kbeg; k<=kend; k++) {
      for(j=jbeg; j<=jend; j++) {
        for(i=ibeg; i<=iend; i++) {
          node     = Grid.node(i,j,k);
          hh       = fmax(hh,Grid.U[node].d);
#pragma acc atomic update
          lbuf[i] += Grid.U[node].d/hsize;  
        }
      }
    }

    MPI_Allreduce(&hh,&rho_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#pragma acc host_data use_device(lbuf, rh_mean)
{
    MPI_Allreduce(lbuf,rh_mean,vsize,MPI_DOUBLE,MPI_SUM,YZ_COMM);
}

    if ( Run.rank == 0 ) {
      cout << "damp: tau_ref = " << tau_ref         << endl;
      cout << "damp: vel_ref = " << vel_ref         << endl;
      cout << "damp: tau_min = " << tau_min         << endl;
      cout << "damp: mflx    = " << vel_ref*rho_max << endl;
    }

    ini_flag = 0;
  }

  //-----------------------------------------------------------------

#pragma acc parallel loop present(lbuf[:vsize])
  for(j=0;j<vsize;j++) lbuf[j] = 0.0;
#pragma acc parallel loop collapse(3) \
 present(Grid[:1], Grid.U[:Grid.bufsize], lbuf[:vsize])
  for(k=kbeg; k<=kend; k++) {
    for(j=jbeg; j<=jend; j++) {
      for(i=ibeg; i<=iend; i++) {
        node = Grid.node(i,j,k);
#pragma acc atomic update
        lbuf[i] += Grid.U[node].M.x*Grid.U[node].d/hsize;  
      }
    }
  }

#pragma acc host_data use_device(lbuf, rlbuf)
{  
  MPI_Allreduce(lbuf,rlbuf,vsize,MPI_DOUBLE,MPI_SUM,YZ_COMM);
}

  hmax  = 0.0; hphot = 0.0; flxmax  = 0.0; vmax  = 0.0;
    
  //for(j=Grid.lbeg[0];j <=Grid.lend[0];j++){
 
#pragma acc parallel loop present(rlbuf[:vsize],my_mean[:vsize])
  for(j=ibeg; j<=iend;j++){
     my_mean[j] = rlbuf[j];
  }
#pragma acc parallel loop \
 reduction(max:hmax) reduction(max:hphot) reduction(max:flxmax) reduction(max:vmax) \
 present(Grid[:1],my_mean[:vsize])
  for(j=ibeg; j<=iend; j++) {

    hh = tau_ref*pow(my_mean[j]/(vel_ref*rho_max),4);
    hh = fmin(tau_min,hh);

    hmax=fmax(hmax,hh);
    if(rh_mean[j] < 1.e-06) hphot=fmax(hphot,hh);
    flxmax=fmax(flxmax,fabs(my_mean[j]));
    vmax=fmax(vmax,fabs(my_mean[j]/rh_mean[j]));

    my_mean[j] *= hh;
  }

  if(Run.verbose >1){
    sbuf[0] = flxmax;
    sbuf[1] = vmax;
    sbuf[2] = hmax;
    sbuf[3] = hphot; 
    MPI_Reduce(sbuf,rbuf,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (Run.rank == 0){
      cout << "damp: flx_max  = " << rbuf[0] << endl;
      cout << "damp: vel_max  = " << rbuf[1] << endl;
      cout << "damp: tau_max  = " << rbuf[2] << endl;
      cout << "damp: tau_phot = " << rbuf[3] << endl;
    }
  }

}

