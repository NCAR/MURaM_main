#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "run.H"
#include "grid.H"
#include "eos_precision.h"

using std::cout;
using std::endl;

int N_eps,N_lr,N_lp3,N_s3;
#pragma acc declare create(N_eps, N_lr, N_lp3, N_s3)
double del_eps,del_lr,del_lp3,del_s3;
double inv_del_eps,inv_del_lr,inv_del_lp3,inv_del_s3;
#pragma acc declare create(inv_del_eps, inv_del_lr)
double eps_off, ss_off;
#pragma acc declare create(ss_off, eps_off, inv_del_lp3, inv_del_s3)

double *xeps;
#pragma acc declare create(xeps)
double *xlr;
#pragma acc declare create(xlr)
eos_real **p_eostab;
#pragma acc declare create(p_eostab)
eos_real **T_eostab;
eos_real **s_eostab;
#pragma acc declare create(s_eostab)
eos_real **ne_eostab; 
eos_real **rhoi_eostab;
eos_real **amb_eostab;

double *xlp3;
#pragma acc declare create(xlp3)
double *xs3;
#pragma acc declare create(xs3)
eos_real **d3_eostab;
#pragma acc declare create(d3_eostab)
eos_real **eps3_eostab;
#pragma acc declare create(eps3_eostab)

eos_real** array_2d_contiguous(int nx,int ny){
  eos_real **array = new eos_real*[nx];

  array[0] = new eos_real[nx*ny];
  
  for(int i=0;i<nx;i++)
    array[i] = array[0]+i*ny;

  return array;
}  

void eos_init(GridData& Grid, RunData& Run, const PhysicsData& Physics) { 

  int i,j,gsize,ind1,ind2,ind3,ind4,ind5,ind6;
  int buffsize;    
  eostab_real* buff;
  FILE* pFile;
  eostab_real header[14];
  double eps0,eps1,lr0,lr1,lp0,lp1,s0,s1;

  if (Run.rank == 0)
      cout << "EOS: Use ln(eps), ln(rho), ln(p), ln(T), s, ln(rhoi),ln(amb),ln(ne) EOS interpolation" << endl;

  gsize = (Grid.lsize[0]+2.*Grid.ghosts[0])*(Grid.lsize[1]+2.*Grid.ghosts[1])*(Grid.lsize[2]+2.*Grid.ghosts[2]);

  // read in tables size
  pFile = fopen(Run.eos_name,"r");

  if(pFile == NULL){
      if (Run.rank == 0) cout << "eos_init: EOS file - " << Run.eos_name << " not found" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
  } else {
     if (Run.rank == 0) cout << "eos_init: Opened file - " << Run.eos_name << endl;
  } 
  fread(header,sizeof(eostab_real),14,pFile);

  N_eps = (int) header[0];
#pragma acc update device(N_eps)
  N_lr =  (int) header[1];
#pragma acc update device(N_lr)
  N_lp3 = (int) header[2];
#pragma acc update device(N_lp3)
  N_s3 =  (int) header[3];
#pragma acc update device(N_s3)

  eps0 = (double) header[4];
  eps1 = (double) header[5];
  lr0 =  (double) header[6];
  lr1 =  (double) header[7];
  lp0 =  (double) header[8];
  lp1 =  (double) header[9];
  s0 =   (double) header[10];
  s1 =   (double) header[11];
  eps_off = (double) header[12];
#pragma acc update device(eps_off)
  ss_off =  (double) header[13];
#pragma acc update device(ss_off)
  //Dynamically allocate EOS tables
  xeps = new double[N_eps];
  xlr  = new double[N_lr];

  p_eostab    = array_2d_contiguous(N_eps,N_lr);
  T_eostab    = array_2d_contiguous(N_eps,N_lr);
  s_eostab    = array_2d_contiguous(N_eps,N_lr);
  ne_eostab   = array_2d_contiguous(N_eps,N_lr);
  rhoi_eostab = array_2d_contiguous(N_eps,N_lr);
  amb_eostab  = array_2d_contiguous(N_eps,N_lr);

  //Dynamically allocate inverse EOS tables
  xlp3 = new double[N_lp3];
  xs3  = new double[N_s3];
  
  d3_eostab   = array_2d_contiguous(N_lp3,N_s3);
  eps3_eostab = array_2d_contiguous(N_lp3,N_s3);

  //Define epsilon and log rho axes values.
  del_eps = (eps1-eps0)/(double (N_eps-1));
  del_lr  = (lr1-lr0)/(double (N_lr-1));

  inv_del_eps = 1.0/del_eps;
#pragma acc update device(inv_del_eps)
  inv_del_lr = 1.0/del_lr;
#pragma acc update device(inv_del_lr)

  for (i=0; i<N_eps; i++) {
    xeps[i] = eps0 + del_eps * ((double) i);
  }
  for (i=0; i<N_lr;  i++) {
    xlr[i] =  lr0 + del_lr * ((double) i);
  }

  //Define logp and entropy axes values.
  del_lp3 = (lp1-lp0)/(double (N_lp3-1));
  del_s3  = (s1-s0)/(double (N_s3-1));

  inv_del_lp3 = 1.0/del_lp3;
#pragma acc update device(inv_del_lp3)
  inv_del_s3  = 1.0/del_s3;
#pragma acc update device(inv_del_s3)

  for (i=0; i<N_lp3; i++) {
    xlp3[i] = lp0+del_lp3 * ((double) i);
  }
#pragma acc enter data copyin(xlp3[:N_lp3])
  for (i=0; i<N_s3;  i++) {
    xs3[i]  = s0+del_s3 * ((double) i);
  }
#pragma acc enter data copyin(xs3[:N_s3])

  // Temporary buffer
  buffsize = 6*N_eps*N_lr;
  buff = new eostab_real[buffsize];

  fread(buff,sizeof(eostab_real),buffsize,pFile);
  
  for (i=0; i< N_eps; i++){ 
    for (j=0; j< N_lr; j++){
      ind1=j*N_eps+i;
      ind2=ind1+N_eps*N_lr;
      ind3=ind2+N_eps*N_lr;
      ind4=ind3+N_eps*N_lr;
      ind5=ind4+N_eps*N_lr;
      ind6=ind5+N_eps*N_lr;
      // Merged 'higher' tables, all tables log except for entropy.
      p_eostab[i][j]    = (eos_real) ((double) buff[ind1])/(del_eps*del_lr);
      T_eostab[i][j]    = (eos_real) ((double) buff[ind2])/(del_eps*del_lr);
      s_eostab[i][j]    = (eos_real) ((double) buff[ind3])/(del_eps*del_lr);
      ne_eostab[i][j]   = (eos_real) ((double) buff[ind4])/(del_eps*del_lr);
      rhoi_eostab[i][j] = (eos_real) ((double) buff[ind5])/(del_eps*del_lr); 
      amb_eostab[i][j]  = (eos_real) ((double) buff[ind6])/(del_eps*del_lr); 
    }
  }

  delete[] buff;

  //===========Reverse tables: eps=eps(p,s) and rho=rho(p,s) ================
  // Recyle temporary buffer
  
  buffsize=2*N_lp3*N_s3;
  buff=new eostab_real[buffsize];

  fread(buff,sizeof(eostab_real),buffsize,pFile);
  
  fclose(pFile);

  for (i=0; i< N_lp3; i++){ 
    for (j=0; j< N_s3; j++){
      ind1=j*N_lp3+i;
      ind2=ind1+N_lp3*N_s3;
      eps3_eostab[i][j] = ((double) buff[ind1])/(del_lp3*del_s3);
      d3_eostab[i][j]   = ((double) buff[ind2])/(del_lp3*del_s3);
    }
  }
#pragma acc enter data copyin(d3_eostab[:N_lp3][:N_s3])
#pragma acc enter data copyin(eps3_eostab[:N_lp3][:N_s3]) 
  delete[] buff;

#pragma acc enter data copyin(xeps[:N_eps])
#pragma acc enter data copyin(xlr[:N_lr])
#pragma acc enter data copyin(p_eostab[:N_eps][:N_lr])
#pragma acc enter data copyin(T_eostab[:N_eps][:N_lr])
#pragma acc enter data copyin(s_eostab[:N_eps][:N_lr])
#pragma acc enter data copyin(ne_eostab[:N_eps][:N_lr])
#pragma acc enter data copyin(rhoi_eostab[:N_eps][:N_lr])
#pragma acc enter data copyin(amb_eostab[:N_eps][:N_lr])

}

/***********************************************************************/
#pragma acc routine seq
double d3_interp(double pp, double ss){
  int i,j;
  double logp, dd,ss1;

  ss1 = ss + ss_off;
  logp = log(pp);

  i = (int) ( (logp-xlp3[0])*inv_del_lp3    );
  if (i < 0)       i=0;
  if (i > N_lp3-2) i=N_lp3-2;

  j = (int) ( (ss1-xs3[0])*inv_del_s3 );
  if (j < 0)      j=0;
  if (j > N_s3-2) j=N_s3-2;

  dd =  (logp-xlp3[i])   * (ss1-xs3[j])   * d3_eostab[i+1][j+1]
       +(logp-xlp3[i])   * (xs3[j+1]-ss1) * d3_eostab[i+1][j]
       +(xlp3[i+1]-logp) * (ss1-xs3[j])   * d3_eostab[i][j+1]
       +(xlp3[i+1]-logp) * (xs3[j+1]-ss1) * d3_eostab[i][j]; 

  dd = exp(dd);

  return dd;  
}

#pragma acc routine seq
double eps3_interp(double pp, double ss){
  int i,j;
  double logp, ee,ss1;

  ss1 = ss + ss_off;
  logp = log(pp);

  i = (int) ( (logp-xlp3[0])*inv_del_lp3    );
  if (i < 0)       i=0;
  if (i > N_lp3-2) i=N_lp3-2;

  j = (int) ( (ss1-xs3[0])*inv_del_s3 );
  if (j < 0)      j=0;
  if (j > N_s3-2) j=N_s3-2;

  ee =  (logp-xlp3[i])   * (ss1-xs3[j])   * eps3_eostab[i+1][j+1]
       +(logp-xlp3[i])   * (xs3[j+1]-ss1) * eps3_eostab[i+1][j]
       +(xlp3[i+1]-logp) * (ss1-xs3[j])   * eps3_eostab[i][j+1]
       +(xlp3[i+1]-logp) * (xs3[j+1]-ss1) * eps3_eostab[i][j]; 

    ee = exp(ee);

  return ee-eps_off;
}

#pragma acc routine seq
double s_interp(double ee, double dd) {

  int i,j;
  double logr, ss,ee1;

   ee1 = log(ee+eps_off);
   logr = log(dd);

  i = (int) ((ee1-xeps[0])*inv_del_eps );
  if (i < 0)         i=0;
  if (i > N_eps - 2) i=N_eps-2;

  j = (int) ( (logr-xlr[0])*inv_del_lr );
  if (j < 0)        j=0;
  if (j > N_lr - 2) j=N_lr - 2;  
 
  ss =  (logr-xlr[j])   * (ee1-xeps[i])   * s_eostab[i+1][j+1]
       +(logr-xlr[j])   * (xeps[i+1]-ee1) * s_eostab[i][j+1]
       +(xlr[j+1]-logr) * (ee1-xeps[i])   * s_eostab[i+1][j]
       +(xlr[j+1]-logr) * (xeps[i+1]-ee1) * s_eostab[i][j]; 
 
  return ss = ss-ss_off;
}

double T_interp(double ee, double dd){

  int i,j;
  double logr, logT,ee1;

  ee1 = log(ee+eps_off);
  logr = log(dd);

  i = (int) ((ee1-xeps[0])*inv_del_eps );
  if (i < 0)         i=0;
  if (i > N_eps - 2) i=N_eps-2;
  
  j = (int) ( (logr-xlr[0])*inv_del_lr );
  if (j < 0)        j=0;
  if (j > N_lr - 2) j=N_lr - 2;  
 
  logT =  (logr-xlr[j])   * (ee1-xeps[i])   * T_eostab[i+1][j+1]
       +(logr-xlr[j])   * (xeps[i+1]-ee1) * T_eostab[i][j+1]
       +(xlr[j+1]-logr) * (ee1-xeps[i])   * T_eostab[i+1][j]
       +(xlr[j+1]-logr) * (xeps[i+1]-ee1) * T_eostab[i][j]; 

  logT = exp(logT);
  
  return logT;
}

#pragma acc routine seq
double p_interp(double ee, double dd) {

  int i,j;
  double logr, logp,ee1;

  ee1 = log(ee+eps_off);
  logr = log(dd);

  i = (int) ((ee1-xeps[0])*inv_del_eps );
  if (i < 0)         i=0;
  if (i > N_eps - 2) i=N_eps-2;
  
  j = (int) ( (logr-xlr[0])*inv_del_lr );
  if (j < 0)        j=0;
  if (j > N_lr - 2) j=N_lr - 2;  
 
  logp =  (logr-xlr[j])   * (ee1-xeps[i])   * p_eostab[i+1][j+1]
       +(logr-xlr[j])   * (xeps[i+1]-ee1) * p_eostab[i][j+1]
       +(xlr[j+1]-logr) * (ee1-xeps[i])   * p_eostab[i+1][j]
       +(xlr[j+1]-logr) * (xeps[i+1]-ee1) * p_eostab[i][j]; 

  logp = exp(logp);
  
  return logp;
}

double ne_interp(double ee, double dd) {

  int i,j;
  double logr, nel,ee1;

  ee1 = log(ee+eps_off);
  logr = log(dd);

  i = (int) ((ee1-xeps[0])*inv_del_eps );
  if (i < 0)         i=0;
  if (i > N_eps - 2) i=N_eps-2;
  
  j = (int) ( (logr-xlr[0])*inv_del_lr );
  if (j < 0)        j=0;
  if (j > N_lr - 2) j=N_lr - 2;  
 
  nel =  (logr-xlr[j])   * (ee1-xeps[i])   * ne_eostab[i+1][j+1]
       +(logr-xlr[j])   * (xeps[i+1]-ee1) * ne_eostab[i][j+1]
       +(xlr[j+1]-logr) * (ee1-xeps[i])   * ne_eostab[i+1][j]
       +(xlr[j+1]-logr) * (xeps[i+1]-ee1) * ne_eostab[i][j]; 

  nel = exp(nel);
  
  return nel;
}
double rhoi_interp(double ee, double dd) {

  double ee1 = log(ee+eps_off);
  double logr = log(dd);

  int i = (int) ((ee1-xeps[0])*inv_del_eps );
  if (i < 0)         i=0;
  if (i > N_eps - 2) i=N_eps-2;
  
  int j = (int) ( (logr-xlr[0])*inv_del_lr );
  if (j < 0)        j=0;
  if (j > N_lr - 2) j=N_lr - 2;  
 
  double rhoil =  (logr-xlr[j])   * (ee1-xeps[i])   * rhoi_eostab[i+1][j+1]
       +(logr-xlr[j])   * (xeps[i+1]-ee1) * rhoi_eostab[i][j+1]
       +(xlr[j+1]-logr) * (ee1-xeps[i])   * rhoi_eostab[i+1][j]
       +(xlr[j+1]-logr) * (xeps[i+1]-ee1) * rhoi_eostab[i][j]; 

  rhoil = exp(rhoil);
  
  return rhoil;
}
double amb_interp(double ee, double dd) {

  int i,j;
  double logr, ambl,ee1;

  ee1 = log(ee+eps_off);
  logr = log(dd);

  i = (int) ((ee1-xeps[0])*inv_del_eps );
  if (i < 0)         i=0;
  if (i > N_eps - 2) i=N_eps-2;
  
  j = (int) ( (logr-xlr[0])*inv_del_lr );
  if (j < 0)        j=0;
  if (j > N_lr - 2) j=N_lr - 2;  
 
  ambl =  (logr-xlr[j])   * (ee1-xeps[i])   * amb_eostab[i+1][j+1]
       +(logr-xlr[j])   * (xeps[i+1]-ee1) * amb_eostab[i][j+1]
       +(xlr[j+1]-logr) * (ee1-xeps[i])   * amb_eostab[i+1][j]
       +(xlr[j+1]-logr) * (xeps[i+1]-ee1) * amb_eostab[i][j]; 

  ambl = exp(ambl);
  
  return ambl;
}
