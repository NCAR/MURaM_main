#include <cmath>
#include <string.h>
#include <fstream>
#include <errno.h>
#include "mem.h"
#include "exchange.H"
#include "rt.h"
#include "comm_split.H"
#include "ACCH.h"
#include <algorithm>

using std::min;
using std::max;
using std::cout;
using std::endl;

namespace
{
  int bound1, bound2, bound3;
  int b1_i, b2_i, b3_i;
  int step1, step2, step3;
  int str1, str2, str3;
  int inustr1, inustr2, inustr3; // coeff stride
}

void RTS::IntegrateSetup(
  int yi_i, int xi_i, int zi_i, int ystep, int xstep, int zstep
)
{
  bool ix, iy, iz;
  ix = (ixstep[0] == 1 && ixstep[1] == 1    &&
        ixstep[2] == 1 && ixstep[3] == 1)   ||
       (ixstep[0] == -1 && ixstep[1] == -1  &&
        ixstep[2] == -1 && ixstep[3] == -1);
  iy = (iystep[0] == 1 && iystep[1] == 1    &&
        iystep[2] == 1 && iystep[3] == 1)   ||
       (iystep[0] == -1 && iystep[1] == -1  &&
        iystep[2] == -1 && iystep[3] == -1);
 iz = (izstep[0] == 1 && izstep[1] == 1    &&
       izstep[2] == 1 && izstep[3] == 1)   ||
      (izstep[0] == -1 && izstep[1] == -1  &&
       izstep[2] == -1 && izstep[3] == -1);

          if(ix) {
            bound1=nx-1; bound2=ny-1; bound3=nz-1;
            b1_i=xi_i; b2_i=yi_i; b3_i=zi_i;
            step1=xstep; step2=ystep; step3=zstep;
            str1=nz; str2=nx*nz; str3=1;
            inustr1=(nz-1); inustr2=(nx-1)*(nz-1); inustr3=1;
          } else if(iy) {
            bound1=ny-1; bound2=nx-1; bound3=nz-1;
            b1_i=yi_i; b2_i=xi_i; b3_i=zi_i;
            step1=ystep; step2=xstep; step3=zstep;
            str1=nx*nz; str2=nz; str3=1;
            inustr1=(nx-1)*(nz-1); inustr2=(nz-1); inustr3=1;
          } else if(iz) {
            bound1=nz-1; bound2=ny-1; bound3=nx-1;
            b1_i=zi_i; b2_i=yi_i; b3_i=xi_i;
            step1=zstep; step2=ystep; step3=xstep;
            str1=1; str2=nx*nz; str3=nz;
            inustr1=1; inustr2=(nx-1)*(nz-1); inustr3=(nz-1);
          } else {
            cout << "Error in RTS::Driver: No valid dependency.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
}

RTS *rt_new(GridData &Grid,RunData &Run,PhysicsData &Physics)
{
  int rttype;
  if(Run.rank==0){ // check solver type only
    rttype = Physics.rt[i_rt_type];
  }

  MPI_Bcast(&rttype,1,MPI_INT,0,MPI_COMM_WORLD);

  // Switch for chromospheric extension, three options here.
  // reads in parameters file
  // RT_DEFAULT - 0 - LTE rt, no scattering, optically thin and NLTE line losses can be incorporated.
  // SCATTER - 1 - rt with multigroup scattering as in Skartlein (2000) and Hayek (2010)
  // if nothing else return error

  switch(rttype){
      case(RT_DEFAULT): return new RTS(Grid,Run,Physics);
/*      case(RT_SCATTER): return new RTS_SCATTER(Grid,Run,Physics); No scattering for now */
      default: return 0;
  }
}


RTS::~RTS(void)
{
  for(int j = 0;j < cart_sizes[2]; j++) delete [] comm_col[j];
  delete [] comm_col;

///*
  ACCH::Free(I_o, ny*nx*sizeof(double)); // RT Grid
  ACCH::Free(tr_switch, nx*ny*nz*sizeof(int)); // switch for the lowest point at which we have reached the transition region
  ACCH::Free(lgTe, nx*ny*nz*sizeof(double)); // temperature
  ACCH::Free(lgPe, nx*ny*nz*sizeof(double)); // temperature
  ACCH::Free(T_ind, nx*ny*nz*sizeof(int)); // temperature
  ACCH::Free(P_ind, nx*ny*nz*sizeof(int)); // temperature
  ACCH::Free(rho, nx*ny*nz*sizeof(double));  // RT grid

  ACCH::Free(Tau, nx*ny*nz*sizeof(double));  // RT grid
  ACCH::Free(Qt, ny-yo*nx-xo*nz-zo*sizeof(double)); // MHD grid

  ACCH::Free(Jt, nx*ny*nz*sizeof(double));
  ACCH::Free(St, nx*ny*nz*sizeof(double)); // RT grid

// frequency dependent quantities
  ACCH::Free(B, nx*ny*nz*sizeof(double));
  ACCH::Free(kap, nx*ny*nz*sizeof(double));
  
  if (rttype ==1){
    ACCH::Free(sig, nx*ny*nz*sizeof(double));
    ACCH::Free(abn, nx*ny*nz*sizeof(double));
  }
  
  ACCH::Free(J_band, nx*ny*nz*sizeof(double));
  ACCH::Free(I_n, nx*ny*nz*sizeof(double));
  ACCH::Free(I_n1, nx*ny*nz*sizeof(double));

  ACCH::Free(Fx, nx*ny*nz*sizeof(double));
  ACCH::Free(Fy, nx*ny*nz*sizeof(double));
  ACCH::Free(Fz, nx*ny*nz*sizeof(double));

  ACCH::Free(coeff, nx*ny*nz*2*sizeof(double));
  ACCH::Free(coeff1, nx*ny*nz*sizeof(double));
  ACCH::Free(coeff2, nx*ny*nz*sizeof(double));

  ACCH::Free(Qtemp, ny*nx*nz*2);
  ACCH::Free2D<double>(sbuf, ny, nx);
  ACCH::Free2D<double>(rbuf, ny, nx);
  ACCH::Free3D<double>(Col_out, Nbands, col_nz, col_nvar);
   
  del_i5dim(numits,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);
  
  if (NDIM>1){
    ACCH::Free(x_sbuf, Nbands*2*2*2*NMU*ny*nz*sizeof(double));
    ACCH::Free(x_rbuf, Nbands*2*2*2*NMU*ny*nz*sizeof(double));
    ACCH::Free(x_oldbuf, Nbands*2*2*2*NMU*ny*nz*sizeof(double));
  }

  if (NDIM==3){
    ACCH::Free(y_sbuf, Nbands*2*2*2*NMU*nx*nz*sizeof(double));
    ACCH::Free(y_rbuf, Nbands*2*2*2*NMU*nx*nz*sizeof(double));
    ACCH::Free(y_oldbuf, Nbands*2*2*2*NMU*nx*nz*sizeof(double));
  }

  ACCH::Free(z_sbuf, Nbands*2*2*2*NMU*nx*ny*sizeof(double));
  ACCH::Free(z_rbuf, Nbands*2*2*2*NMU*nx*ny*sizeof(double));
  ACCH::Free(z_oldbuf, Nbands*2*2*2*NMU*nx*ny*sizeof(double));

  ACCH::Free(tab_T, NT*sizeof(double));
  ACCH::Free(tab_p, Np*sizeof(double));
  ACCH::Free(invT_tab, NT*sizeof(double));
  ACCH::Free(invP_tab, Np*sizeof(double));
  
  if(fullodf) {
    ACCH::Free(nu_tab, Nlam*sizeof(float));
  }

  if (N5000){
    ACCH::Free(B_5000_tab, NT*sizeof(float));
    ACCH::Free(kap_5000_tab, NT*Np*sizeof(float));
  }

  ACCH::Free(kap_tab, Nbands*NT*Np*sizeof(float));
  ACCH::Free(B_tab, Nbands*NT*sizeof(float));

  ACCH::Free(I_band, nx*ny);

  ACCH::Free(tI_n, nx*ny*nz*sizeof(double));
  ACCH::Free(trho, nx*ny*nz*sizeof(double));
  ACCH::Free(tkap, nx*ny*nz*sizeof(double));
  ACCH::Free(tB, nx*ny*nz*sizeof(double));



  ACCH::Delete(this, sizeof(RTS));
  //fprintf(stdout,"deleting RTS mem");
}
 
double RTS::tau(int z,int x,int y){
  int ind = y*nx*nz + x*nz + z;
  int y_o = yo*nx*nz;
  int x_o = xo*nz;
  int z_o = zo;
  double Tau_local = Tau[ind] +
                     Tau[ind - y_o] +
                     Tau[ind - x_o] +
                     Tau[ind - z_o] +
                     Tau[ind - y_o - x_o] +
                     Tau[ind - y_o - z_o] +
                     Tau[ind - x_o - z_o] +
                     Tau[ind - y_o - x_o - z_o];
  return Tau_local * 0.125;
}

// Unused
double RTS::Qtot(int z,int x,int y)
{  
  //return Qt[y-yl-yo][x-xl-xo][z-zl-zo];
  return Qt[(((y-yl-yo)*(nx-xo)+(x-xl-xo))*(nz-zo)+(z-zl-zo))];
}


double RTS::Jtot(int z,int x,int y){
  int ind = y*nx*nz + x*nz + z;
  int y_o = yo*nx*nz;
  int x_o = xo*nz;
  int z_o = zo;
  double J_local = Jt[ind] +
                   Jt[ind - y_o] +
                   Jt[ind - x_o] +
                   Jt[ind - z_o] +
                   Jt[ind - y_o - x_o] +
                   Jt[ind - y_o - z_o] +
                   Jt[ind - x_o - z_o] +
                   Jt[ind - y_o - x_o - z_o];
  return J_local * 0.125;
}

double RTS::Stot(int z,int x,int y){
  double S_local = St[y*nx*nz+x*nz+z] +
                   St[(y-yo)*nx*nz+x*nz+z] +
                   St[y*nx*nz+(x-xo)*nz+z] +
                   St[y*nx*nz+x*nz+(z-zo)] +
                   St[(y-yo)*nx*nz+(x-xo)*nz+z] +
                   St[(y-yo)*nx*nz+x*nz+(z-zo)] +
                   St[y*nx*nz+(x-xo)*nz+(z-zo)] +
                   St[(y-yo)*nx*nz+(x-xo)*nz+(z-zo)];
  S_local *= 0.125;

  return S_local;
}

void RTS::UpdateIout()
{
  ACCH::UpdateCPU(I_o, nx*ny*sizeof(double));
}

double RTS::Iout(int x,int y)
{
  double io = I_o[(y-yl)*nx+(x-xl)] + I_o[(y-yo-yl)*nx+(x-xl)] + I_o[(y-yl)*nx+(x-xo-xl)] + I_o[(y-yo-yl)*nx+(x-xo-xl)];
  io *= 0.25;
  return io;
}

RTS::RTS(GridData&Grid,RunData &Run,PhysicsData &Physics){
 
  NDIM = Grid.NDIM;
  rttype = Physics.rt[i_rt_type];
  eps_const = Physics.rt[i_rt_epsilon];

  verbose=Run.verbose;

  call_count=0;

  xl=Grid.lbeg[1];
  xh=Grid.lend[1];
  yl=Grid.lbeg[2];
  yh=Grid.lend[2];
  zl=Grid.lbeg[0];
  zh=Grid.lend[0];
//
  xl-=(xo=(Grid.lbeg[1]-Grid.lend[1])?1:0);
  yl-=(yo=(Grid.lbeg[2]-Grid.lend[2])?1:0);
  zl-=(zo=1);

//
  nx=xh-xl+1;
  ny=yh-yl+1;
  nz=zh-zl+1;
  myrank=Run.rank;

  for(int v=0;v<3;v++)
    next[v] = Grid.stride[v];

  int ndim=3,cart_periods[3];  

  MPI_Cart_get(cart_comm,ndim,cart_sizes,cart_periods,lrank);
  MPI_Cart_coords(cart_comm, Run.rank, ndim, lrank);
 
  for(int nleft,nright,nd=0;nd<ndim;nd++){
    MPI_Cart_shift(cart_comm,nd,1,&nleft,&nright);
    leftr[nd] = nleft;
    rightr[nd] = nright;
    
    isgbeg[nd]=Grid.is_gbeg[nd];
    isgend[nd]=Grid.is_gend[nd];
  }

  int ***colranks=i3dim(0,cart_sizes[2]-1,0,cart_sizes[1]-1,0,cart_sizes[0]-1);
  for(int j=0; j<cart_sizes[2];j++)
    for(int k=0; k<cart_sizes[1];k++)
      for(int i=0; i<cart_sizes[0];i++){
        int coords[3]={i,k,j};
        MPI_Cart_rank(cart_comm,coords,&(colranks[j][k][i]));
      }

  char carlson_fp[16];

  if (NMU == 3)
    strcpy(carlson_fp,"./carlson3.dat");
  else if (NMU == 6)
    strcpy(carlson_fp,"./carlson6.dat");
  else if (NMU == 10)
    strcpy(carlson_fp,"./carlson10.dat");
  else {
    fprintf(stdout," NMU %i not supported, Carlson 3,6,10 currently included.\n ",NMU);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  std::ifstream fptr(carlson_fp,std::ios::in);

  if(fptr){
    fptr.precision(16);
    for (int i=0; i < NMU; i++)
      fptr >> xmu[0][i] >> xmu[1][i] >> xmu[2][i] >> wmu[i];
    fptr.close();
  } else {
    fprintf(stdout,"carlson%i.dat not found. Aborting \n",NMU);
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int* ranks=new int [cart_sizes[0]];
  MPI_Group MPI_GROUP_WORLD;

  MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);
  MPI_Group ** grp_col = new MPI_Group * [cart_sizes[2]];
  comm_col = new MPI_Comm * [cart_sizes[2]];

  for(int j=0; j<cart_sizes[2];j++){
    grp_col[j]=new MPI_Group [cart_sizes[1]];
    comm_col[j] = new MPI_Comm [cart_sizes[1]];
  }

  for(int j=0; j< cart_sizes[2]; j++){
    for(int k=0; k< cart_sizes[1]; k++){
      for(int i=0; i< cart_sizes[0]; i++)
        ranks[i]=colranks[j][k][cart_sizes[0]-1-i];
      MPI_Group_incl(MPI_GROUP_WORLD,cart_sizes[0],ranks,&(grp_col[j][k]));
      MPI_Comm_create(MPI_COMM_WORLD,grp_col[j][k],&(comm_col[j][k]));
    }
  }


  // load kappa bins

  ACCH::Copyin(this, sizeof(RTS));

  load_bins(Run.kap_name);

// output intensity
  I_o = (double*) ACCH::Malloc(nx*ny*sizeof(double));        // RT grid

  Fr_mean = (double*) ACCH::Malloc(Nbands*sizeof(double));
  gFr_mean = (double*) ACCH::Malloc(Nbands*sizeof(double));
  memset(Fr_mean,0,Nbands*sizeof(double));
  memset(gFr_mean,0,Nbands*sizeof(double));

  lgTe = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // temperature
  lgPe = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // pressure
  T_ind = (int*) ACCH::Malloc(nx*ny*nz*sizeof(int)); // pressure
  P_ind = (int*) ACCH::Malloc(nx*ny*nz*sizeof(int)); // pressure
  rho = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // RT grid
  tr_switch = (int*) ACCH::Malloc(nx*ny*nz*sizeof(int));

  Tau = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));  // RT grid

  Qt = (double*) ACCH::Malloc((ny-yo)*(nx-xo)*(nz-zo)*sizeof(double)); // MHD grid
  Jt = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // MHD grid
  St = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double)); // MHD grid
  sbuf = (double**) ACCH::Malloc2D<double>(ny, nx);
  rbuf = (double**) ACCH::Malloc2D<double>(ny, nx);
  Qtemp = (double*) ACCH::Malloc(ny*nx*nz*2*sizeof(double));
// frequency dependent quantities
  B = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  kap = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  I_n = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  I_n1 = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  if (rttype==0){
    sig=kap;
    abn=kap;
  } else {
    sig = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
    abn = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  }

  J_band = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  Fx = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  Fy = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  Fz = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  coeff = (double*) ACCH::Malloc(nx*ny*nz*2*sizeof(double));
  coeff1 = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  coeff2 = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  /* Column outputs
   * 0: J_col - Angle averaged intensity
   * 1: S_col - Total Source function (=B in LTE)
   * 2: kap_col - total opacity
   * 3: abs_col - absorption opacity (= kap in LTE)
   * 4: sig_col - scattering opacity (= sig in LTE)
   * 5: B_col - Planck function (=S in LTE)
   * 6: tau_col - band dependant tau
   * 7: Qj_col - Q from radiative energy imbalance
   * 8: Qf_col - Q from flux divergence
   */

  col_offz=Grid.beg[0]-Grid.gbeg[0];
  col_nz = Grid.lsize[0];
  col_nzt = Grid.gsize[0];
  col_nvar = 9;

  Col_out = ACCH::Malloc3D<double>(Nbands, col_nz, col_nvar);

  memset(Col_out[0][0],0,col_nz*Nbands*col_nvar*sizeof(double));

  numits = i5dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);

  if (NDIM==3){
    y_sbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*nx*nz*sizeof(double));
    y_rbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*nx*nz*sizeof(double));
    y_oldbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*nx*nz*sizeof(double));
  }

  if (NDIM>1){
    x_sbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*ny*nz*sizeof(double));
    x_rbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*ny*nz*sizeof(double));
    x_oldbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*ny*nz*sizeof(double));
  }
  z_sbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*ny*nx*sizeof(double));
  z_rbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*ny*nx*sizeof(double));
  z_oldbuf = (double*) ACCH::Malloc(Nbands*2*2*2*NMU*ny*nx*sizeof(double));
 
  double DX=Grid.dx[1];
  double DY=Grid.dx[2];
  double DZ=Grid.dx[0];
  int perm[3][3]={{0,1,2},{1,2,0},{2,0,1}};
// Winkelabhaengige Koeff.
  for(int l=0;l<NMU;++l){ //  Koeff. fuer bilineare Grundflaecheninterp.
    double aM[]={DX/xmu[0][l],DY/xmu[1][l],DZ/xmu[2][l]};
    
    ibase[l]=0;
    ds_upw[l]=DX/xmu[0][l];

    if(((aM[2]<=aM[0])||(NDIM==1))&&((aM[2]<=aM[1])||(NDIM<3))){
      ibase[l]=2;
      ds_upw[l]=DZ/xmu[2][l];  
    }else{
      if((aM[1]<=aM[0])&&(NDIM==3)){
        ibase[l]=1;
        ds_upw[l]=DY/xmu[1][l];
      }
    }
    
    for(int i=0;i<=2;++i){
      a_00[i][l]=(1.0-(aM[perm[i][0]]/aM[perm[i][1]]))*(1.0-(aM[perm[i][0]]/aM[perm[i][2]]));
      a_01[i][l]=(aM[perm[i][0]]/aM[perm[i][2]])*(1.0-(aM[perm[i][0]]/aM[perm[i][1]]));
      a_10[i][l]=(aM[perm[i][0]]/aM[perm[i][1]])*(1.0-(aM[perm[i][0]]/aM[perm[i][2]]));
      a_11[i][l]=(aM[perm[i][0]]/aM[perm[i][1]])*(aM[perm[i][0]]/aM[perm[i][2]]);
    }
  } 



  if (NDIM>1){
#pragma acc parallel loop collapse(6) \
 present(this[:1], x_sbuf[:Nbands*2*2*2*NMU*ny*nz], \
  x_rbuf[:Nbands*2*2*2*NMU*ny*nz], x_oldbuf[:Nbands*2*2*2*NMU*ny*nz])
    for(int band = 0; band < Nbands; band++)
      for(int y = 0; y < 2; y++)
        for(int x = 0; x < 2; x++)
          for(int z = 0; z < 2; z++)
            for(int l = 0; l < NMU; l++)
              for(int i = 0; i < ny*nz; i++) {
                int adr = (((((band*2+y)*2+x)*2+z)*NMU+l)*ny*nz+i);
                x_sbuf[adr] = 0.0;
                x_rbuf[adr] = 0.0;
                x_oldbuf[adr] = 0.0;
              }
  }

  if (NDIM==3){
#pragma acc parallel loop collapse(6) \
 present(this[:1], y_sbuf[:Nbands*2*2*2*NMU*nx*nz], \
  y_rbuf[:Nbands*2*2*2*NMU*nx*nz], y_oldbuf[:Nbands*2*2*2*NMU*nx*nz])
    for(int band = 0; band < Nbands; band++)
      for(int y = 0; y < 2; y++)
        for(int x = 0; x < 2; x++)
          for(int z = 0; z < 2; z++)
            for(int l = 0; l < NMU; l++)
              for(int i = 0; i < nx*nz; i++) {
                int adr = (((((band*2+y)*2+x)*2+z)*NMU+l)*nx*nz+i);
                y_sbuf[adr] = 0.0;
                y_rbuf[adr] = 0.0;
                y_oldbuf[adr] = 0.0;
              }
  }

  {
#pragma acc parallel loop collapse(6) \
 present(this[:1], z_sbuf[:Nbands*2*2*2*NMU*nx*ny], \
  z_rbuf[:Nbands*2*2*2*NMU*nx*ny], z_oldbuf[:Nbands*2*2*2*NMU*nx*ny])
    for(int band = 0; band < Nbands; band++)
      for(int y = 0; y < 2; y++)
        for(int x = 0; x < 2; x++)
          for(int z = 0; z < 2; z++)
            for(int l = 0; l < NMU; l++)
              for(int i = 0; i < nx*ny; i++) {
                int adr = (((((band*2+y)*2+x)*2+z)*NMU+l)*ny*nx+i);
                z_sbuf[adr] = 0.0;
                z_rbuf[adr] = 0.0;
                z_oldbuf[adr] = 0.0;
              }
  }

  memset(numits[0][UP][RIGHT][FWD],0,2*2*2*Nbands*NMU*sizeof(int));
/*
  // init for qrad_tauscale
  int* ranks=new int [cart_sizes[0]];
  MPI_Group MPI_GROUP_WORLD;

  MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);
  MPI_Group ** grp_col = new MPI_Group * [cart_sizes[2]];
  comm_col = new MPI_Comm * [cart_sizes[2]];
  
  for(int j=0; j<cart_sizes[2];j++){
    grp_col[j]=new MPI_Group [cart_sizes[1]];
    comm_col[j] = new MPI_Comm [cart_sizes[1]];
  }
  
  for(int j=0; j< cart_sizes[2]; j++){
    for(int k=0; k< cart_sizes[1]; k++){
      for(int i=0; i< cart_sizes[0]; i++)
        ranks[i]=colranks[j][k][cart_sizes[0]-1-i];
      MPI_Group_incl(MPI_GROUP_WORLD,cart_sizes[0],ranks,&(grp_col[j][k]));
      MPI_Comm_create(MPI_COMM_WORLD,grp_col[j][k],&(comm_col[j][k]));
    }
  }
*/
  delete[] ranks;
// all done...?
  del_i3dim(colranks,0,cart_sizes[2]-1,0,cart_sizes[1]-1,0,cart_sizes[0]-1);

  for(int j = 0;j < cart_sizes[2]; j++) delete [] grp_col[j];
  delete [] grp_col;

  /* if (xl or xh) lies in my x-domain and yl or yh lies in my y-domain
   * then save_col, set limits */

  int gxl = Grid.beg[1];
  int gxh = Grid.end[1];
  int gyl = Grid.beg[2];
  int gyh = Grid.end[2];

  int sl_r[4] = {Grid.gbeg[1],Grid.gend[1],Grid.gbeg[2],Grid.gend[2]};

  int x_range = (((gxl >= sl_r[0])&&(gxl <= sl_r[1]))
      ||((gxh >=sl_r[0])&&(gxh <= sl_r[1])));
  int y_range = (((gyl >= sl_r[2])&&(gyl <= sl_r[3]))
      ||((gyh >=sl_r[2])&&(gyh <=sl_r[2])));

  save_col = 0;
  col_bnd[0] =1;
  col_bnd[1] =0;
  col_bnd[2] =1;
  col_bnd[3] =0;

  num_col = (sl_r[1]-sl_r[0]+1)*(sl_r[3]-sl_r[2]+1);
  avg_col = 1.0/((double) num_col);

  if ((x_range)&&(y_range)){
    //save_col=1;
    col_bnd[0] =min(gxh-sl_r[0],xl+xo);
    col_bnd[1] =max(gxh-sl_r[1],xh);
    col_bnd[2] =min(gyh-sl_r[2],yl+yo);
    col_bnd[3] =max(gyh-sl_r[3],yh);
  }

  I_band = (double*) ACCH::Malloc(nx*ny*sizeof(double));

  tI_n = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  trho = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  tkap = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));
  tB = (double*) ACCH::Malloc(nx*ny*nz*sizeof(double));

  ACCH::UpdateGPU(this, sizeof(RTS));
#pragma acc enter data copyin(this[:1], I_band[:nx*ny], I_o[:nx*ny],Grid[:1], Grid.temp[:Grid.bufsize], \
         Grid.pres[:Grid.bufsize], \
         tr_switch[:nx*ny*nz], lgTe[:nx*ny*nz], lgPe[:nx*ny*nz], \
         rho[:nx*ny*nz], tab_T[:NT], tab_p[:Np], T_ind[:nx*ny*nz], P_ind[:nx*ny*nz], \
         B[:nx*ny*nz], invT_tab[:NT], invP_tab[:Np], \
         B_tab[:Nbands*NT], kap[:nx*ny*nz], kap_tab[:Nbands*NT*Np], \
         Qt[:(ny-yo)*(nx-xo)*(nz-zo)],St[:nx*ny*nz], Jt[:nx*ny*nz],  \
         z_rbuf[:Nbands][:2][:2][:2][:NMU][:nx*ny],sbuf[:ny][:nx], rbuf[:ny][:nx],Qtemp[:ny*nx*nz*2], \
         tI_n[:nx*ny*nz], trho[:nx*ny*nz], tkap[:nx*ny*nz], tB[:nx*ny*nz])

}
void RTS::load_bins(char* kap_name){
  
  std::ifstream fp_rt(kap_name, std::ios::in | std::ios::binary);

  int pt_rhot;
  int junk4;

  if (fp_rt.is_open()) {

    fp_rt.read((char*)&N5000, sizeof(int));
    fp_rt.read((char*)&NT, sizeof(int));
    fp_rt.read((char*)&Np, sizeof(int));
    fp_rt.read((char*)&Nbands, sizeof(int));
    fp_rt.read((char*)&pt_rhot, sizeof(int));
    fp_rt.read((char*)&fullodf, sizeof(int));
    fp_rt.read((char*)&scatter, sizeof(int));
    fp_rt.read((char*)&junk4, sizeof(int)); 
     
    if (pt_rhot!=0){
      fprintf(stdout,"pt_rhot = %i, but rho-T bins are not currently implemented, aborting",pt_rhot);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    if (myrank == 0) {
      fprintf(stdout,"Reading in RTS Header, current settings:\n");
      fprintf(stdout,"Full ODF is %d and scatter is %d \n",fullodf,scatter);
      fprintf(stdout,"RT bins are NT %d Np %d Nbands %d \n",NT,Np,Nbands);
      fprintf(stdout,"Reference bin is %d and coronal back heating bin %d \n",N5000,junk4);
    }
    
    if ((scatter==0)&&(eps_const==0)&&(rttype==1)){
      fprintf(stdout, "rt_type =%i, but scattering bins = %i and photon destruction probability is %e. \n",scatter, rttype,eps_const);
      fprintf(stdout, "For scattering either scattering bins, or constant photon destruction probability are required, aborting. \n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    tab_T = (double*) ACCH::Malloc(NT*sizeof(double));
    tab_p = (double*) ACCH::Malloc(Np*sizeof(double));

    invT_tab = (double*) ACCH::Malloc(NT*sizeof(double));
    invP_tab = (double*) ACCH::Malloc(Np*sizeof(double));
    
    fp_rt.read((char*)&tab_T[0],NT*sizeof(double));
    fp_rt.read((char*)&tab_p[0],Np*sizeof(double));
    
    for (int i=0;i<NT; i++)
      tab_T[i] = tab_T[i]*TENLOG;
     
    for (int j=0;j<Np; j++)
      tab_p[j] = tab_p[j]*TENLOG;

    for (int l=0; l<=NT-2; l++)
      invT_tab[l]= 1./(tab_T[l+1] - tab_T[l]);
    for (int m=0; m<=Np-2; m++)
      invP_tab[m] = 1./ (tab_p[m+1]- tab_p[m]);

    if (N5000){
      kap_5000_tab = (float*) ACCH::Malloc(NT*Np*sizeof(float));
      B_5000_tab = (float*) ACCH::Malloc(NT*sizeof(float));
       
      fp_rt.read((char*)&kap_5000_tab[0],NT*Np*sizeof(float));
      fp_rt.read((char*)&B_5000_tab[0],NT*sizeof(float));
    }
      
    kap_tab = (float*) ACCH::Malloc(Nbands*NT*Np*sizeof(float));
    B_tab = (float*) ACCH::Malloc(Nbands*NT*sizeof(float));
    fp_rt.read((char*)&kap_tab[0],Nbands*NT*Np*sizeof(float));
    fp_rt.read((char*)&B_tab[0],Nbands*NT*sizeof(float));

    if(fullodf){
      nu_tab = f1dim(0,Nlam-1);
      fp_rt.read((char*)&nu_tab[0],Nlam*sizeof(float));

      if (scatter>0){
        acont_pT = f3dim(0,Nlam-1,0,NT-1,0,Np-1);
        kcont_pT = f3dim(0,Nlam-1,0,NT-1,0,Np-1);
        fp_rt.read((char*)&acont_pT[0][0][0],Nlam*NT*Np*sizeof(float));
        fp_rt.read((char*)&kcont_pT[0][0][0],Nlam*NT*Np*sizeof(float));
        
        Npp=1;
        tau_pp_tab = d1dim(0,Npp-1);
        invtau_pp_tab = d1dim(0,Npp-1);

        tau_pp_tab[0] = -99.0;
        invtau_pp_tab[0] = 1.0;

        // If scattering bins and full ODF, but we want to run non-scattering RT,
        // we need to add acont to kappa
        if (rttype==0)
          for (int lam = 0 ;lam < Nlam;lam++)
            for (int bin = 0;bin<Nbin;bin++)
              for (int tt = 0;tt<NT;tt++)
                for (int pp = 0;pp<Np;pp++)
                  kap_tab[((Nbin*lam+bin)*NT+tt)*Np+pp]
                        = log(exp(kap_tab[((Nbin*lam+bin)*NT+tt)*Np+pp])+acont_pT[lam][tt][pp]);
      }
    } else if (scatter>0){
        // output
        Npp = scatter;
        scatter=1;

        sig_tab = f3dim(0,Nbands-1,0,NT-1,0,Np-1);
        abn_tab = f3dim(0,Nbands-1,0,NT-1,0,Np-1);
        fp_rt.read((char*)&abn_tab[0][0][0],Nbands*NT*Np*sizeof(float));
        fp_rt.read((char*)&sig_tab[0][0][0],Nbands*NT*Np*sizeof(float));

        tau_pp_tab = d1dim(0,Npp-1);
        invtau_pp_tab = d1dim(0,Npp-1);

        fp_rt.read((char*)&tau_pp_tab[0],Npp*sizeof(double));

        for (int z=0;z<=Npp-2; z++)
          invtau_pp_tab[z] = 1./(tau_pp_tab[z+1]- tau_pp_tab[z]);
        
        kap_pp_tab = f2dim(0,Nbands-1,0,Npp-1);
        abn_pp_tab = f2dim(0,Nbands-1,0,Npp-1);
        sig_pp_tab = f2dim(0,Nbands-1,0,Npp-1);

        fp_rt.read((char*)&abn_pp_tab[0][0],Nbands*Npp*sizeof(float));
        fp_rt.read((char*)&kap_pp_tab[0][0],Nbands*Npp*sizeof(float));
        fp_rt.read((char*)&sig_pp_tab[0][0],Nbands*Npp*sizeof(float));
    }
        
    if ((scatter==0)&&(rttype==1)){
      Npp=1;
      tau_pp_tab = d1dim(0,Npp-1);
      invtau_pp_tab = d1dim(0,Npp-1);

      tau_pp_tab[0] = -99;
      invtau_pp_tab[0] = 1;
    }

    fp_rt.close();
     
  } else {
    cout << "rt_init: kappa file not found: " << kap_name << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  ACCH::UpdateGPU(tab_T, NT*sizeof(double));
  ACCH::UpdateGPU(tab_p, Np*sizeof(double));
  ACCH::UpdateGPU(invT_tab, NT*sizeof(double));
  ACCH::UpdateGPU(invP_tab, Np*sizeof(double));
  ACCH::UpdateGPU(B_tab, Nbands*NT*sizeof(float));
  ACCH::UpdateGPU(kap_tab, Nbands*NT*Np*sizeof(float));
  if(N5000) {
    ACCH::UpdateGPU(B_5000_tab, NT*sizeof(float));
    ACCH::UpdateGPU(kap_5000_tab, NT*Np*sizeof(float));
  }
}

double RTS::wrapper(int rt_upd,GridData &Grid,RunData &Run,const PhysicsData &Physics){
#pragma acc enter data copyin(Qt[:(ny-yo)*(nx-xo)*(nz-zo)])
//#pragma acc enter data copyin(sbuf[:ny][:nx], rbuf[:ny][:nx],Qtemp[:ny][:nx][:nz][:2]) async
//#pragma acc enter data copyin(this[:1], I_band[:nx*ny], I_o[:nx*ny],Grid[:1], Grid.temp[:Grid.bufsize], \
         Grid.pres[:Grid.bufsize], \
         tr_switch[:nx*ny*nz], lgTe[:nx*ny*nz], lgPe[:nx*ny*nz], \
         rho[:nx*ny*nz], tab_T[:NT], tab_p[:Np], T_ind[:nx*ny*nz], P_ind[:nx*ny*nz], \
         B[:nx*ny*nz], invT_tab[:NT], invP_tab[:Np], \
         B_tab[:Nbands*NT], kap[:nx*ny*nz], kap_tab[:Nbands][:NT][:Np], \
         Qt[:(ny-yo)][:(nx-xo)][:(nz-zo)],St[:nx*ny*nz], Jt[:nx*ny*nz],z_rbuf[:Nbands][:2][:2][:2][:NMU][:nx*ny]) async

#pragma acc parallel loop \
 present(this[:1], I_band[:nx*ny], I_o[:nx*ny])
  for(int i = 0; i < nx*ny; i++) {
    I_band[i] = 0.0;
    I_o[i] = 0.0;
  }
  double DX=Grid.dx[1],DZ=Grid.dx[0],DY=Grid.dx[2];
  int cont_bin = Physics.rt[i_rt_iout];

  const double Temp_TR = Physics.rt[i_rt_tr_tem];
  const double Pres_TR = Physics.rt[i_rt_tr_pre];

  need_I = 0;

  if (Run.NeedsSlice())
    need_I = 1;
  if((Run.iteration%rt_upd)&&(dt_rad>0)){
    if (cont_bin == 2){
    for (int band=Nbands-1;band>=0;--band){
      get_Tau_and_Iout(Grid, Run, Physics,DZ,&B_tab[band*NT],&kap_tab[band*NT],I_band,need_I);
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] +=I_band[y*nx+x];
      }
    }

    if (cont_bin==0){
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_tab,kap_tab,I_band,need_I);
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] +=I_band[y*nx+x];
    }


    if (N5000){
      int I5000_out = 0;
      if ((cont_bin==1)&&(need_I==1))
        I5000_out = 1;
      
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,I_band,I5000_out);
     
      if ((cont_bin==1)&&(need_I==1)) {
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
        for (int y=0;y<ny;y++)
          for (int x=0;x<nx;x++)
            I_o[y*nx+x] += I_band[y*nx+x];
      }
    }

    calc_Qtot_and_Tau(Grid, Run, Physics);
    return dt_rad;
  }


// *****************************************************************
// *        interpolate opacity and Planck function (B)            *
// *****************************************************************

  cState *U=Grid.U;
  const int bufsize = Grid.bufsize;
#pragma acc enter data copyin(U[:bufsize])
  double N = pow(2,NDIM);

  //const double Temp_TR = Physics.rt[i_rt_tr_tem];
  //const double Pres_TR = Physics.rt[i_rt_tr_pre];

  //ACCH::UpdateGPU(Grid.temp, Grid.bufsize*sizeof(double));
  //ACCH::UpdateGPU(Grid.pres, Grid.bufsize*sizeof(double));
  //ACCH::UpdateGPU(U, Grid.bufsize*sizeof(cState));
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Grid[:1], Grid.temp[:Grid.bufsize], \
         Grid.pres[:Grid.bufsize], U[:Grid.bufsize], \
         tr_switch[:nx*ny*nz], lgTe[:nx*ny*nz], lgPe[:nx*ny*nz], \
         rho[:nx*ny*nz], tab_T[:NT], tab_p[:Np], T_ind[:nx*ny*nz], P_ind[:nx*ny*nz], \
         next[1:2]) 
  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
      int y0 = y+yl;
      int x0 = x+xl;
      int off0 = x0*next[1]+y0*next[2];
      int off1 = x0*next[1]+(y0+yo)*next[2];
      int off2 = (x0+xo)*next[1]+y0*next[2];
      int off3 = (x0+xo)*next[1]+(y0+yo)*next[2];
      int xyoff = y*nx*nz+x*nz;
#pragma acc loop vector
      for(int z=0;z<nz;++z){
        int z0 = z+zl;
        int inode[]={off0+z0,off0+z0+zo,off2+z0,off2+z0+zo,off1+z0,off1+z0+zo,off3+z0,off3+z0+zo};
        double Tm=0.0,pm=0.0,rm=0.0;
#pragma acc loop seq
        for(int l=0;l<N;++l){
          Tm+=Grid.temp[inode[l]];
          pm+=Grid.pres[inode[l]];
          rm+=U[inode[l]].d;
        }

        Tm /= N;
        pm /= N;
        rm /= N;

        int ind = xyoff + z;

        //disbale RT if Temp > Temp_TR above the photosphere
        double pswitch = min(max(pm-Pres_TR,0.0),1.0);
        double tswitch = min(max(Temp_TR-Tm,0.0),1.0);
        tr_switch[ind] = (int) max(pswitch,tswitch);
        
        lgTe[ind]=log(Tm);
        lgPe[ind]=log(pm);
        rho[ind] =rm;
    
        lgTe[ind] = min(max(lgTe[ind],(double) tab_T[0]),(double) tab_T[NT-1]);
        lgPe[ind] = min(max(lgPe[ind],(double) tab_p[0]),(double) tab_p[Np-1]);
      }
#pragma acc loop vector 
      for(int z=0;z<nz;++z){

        int ind = xyoff + z;

        // Search table once for all RT bands
        int l=0;
        int m=0;
        if(lgTe[ind]<tab_T[0])
          l=0;
        else if(lgTe[ind]>tab_T[NT-1])
          l=NT-2;
        else {
#pragma acc loop seq 
	for (int li=0; li<=NT-2; li++){
            int lflag = 0;
            if ((lgTe[ind] >= tab_T[li]) && (lgTe[ind] <= tab_T[li+1]))
              lflag = lflag+ 1;
            if(lflag==1)
              l = li;
          }
        }

        if(lgPe[ind]<tab_p[0])
          m=0;
        else if(lgPe[ind]>tab_p[Np-1])
          m=Np-2;
        else {
#pragma acc loop seq 
	for (int mi=0; mi<=Np-2; mi++){
             int mflag = 0;
             if ((lgPe[ind] >= tab_p[mi]) && (lgPe[ind] <= tab_p[mi+1]))
               mflag = mflag+1;
             if(mflag==1)
               m=mi;
          }
        }

        T_ind[ind] = l;
        P_ind[ind] = m;

      }
    }
  }
  
  // Begin RT loop, work band by band,
#pragma acc parallel loop \
 present(this[:1], St[:nx*ny*nz], Jt[:nx*ny*nz])
  for(int i = 0; i < nx*ny*nz; i++) {
    St[i] = 0.0;
    Jt[i] = 0.0;
  }
#pragma acc parallel loop collapse(3) \
 present(this[:1], Qt[:(ny-yo)*(nx-xo)*(nz-zo)])
  for(int y = 0; y < ny-yo; y++)
    for(int x = 0; x < nx-xo; x++)
      for(int z = 0; z < nz-zo; z++) {
        Qt[((y*(nx-xo)+x)*(nz-zo)+z)] = 0.0;
      }
  for(int band=Nbands-1;band>=0;--band){

    if(fullodf){
      nu_ind = band/Nbin;
      bin_ind = band%Nbin;
    }

    if((myrank==0) && (verbose>1)){
      fprintf(stdout,"rt running for band %i of %i \n",band+1,Nbands);
      if (fullodf)
        fprintf(stdout,"nu band = %i, %e and bin %i \n", nu_ind, nu_tab[nu_ind],bin_ind);
    }
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], lgTe[:nx*ny*nz], lgPe[:nx*ny*nz], \
         T_ind[:nx*ny*nz], P_ind[:nx*ny*nz], \
         B[:nx*ny*nz], tab_T[:NT], tab_p[:Np], invT_tab[:NT], invP_tab[:Np], \
         B_tab[:Nbands*NT], kap[:nx*ny*nz], kap_tab[:Nbands*NT*Np], \
         tr_switch[:nx*ny*nz])
    for(int y=0;y<ny;++y){
      for(int x=0;x<nx;++x){
        int xyoff = y*nx*nz + x*nz;
#pragma acc loop vector
        for(int z=0;z<nz;++z){
          int ind = xyoff + z;
      
          int l = T_ind[ind];
          int m = P_ind[ind];
      
          double xt = (lgTe[ind]-tab_T[l])*invT_tab[l];
          double xp = (lgPe[ind]-tab_p[m])*invP_tab[m];
      
          // Interpolate for kappa and B
          B[ind]=exp(xt*B_tab[band*NT+l+1]+(1.-xt)*B_tab[band*NT+l]);
      
          kap[ind] = 
            exp(xt*(xp*kap_tab[(band*NT+l+1)*Np+m+1]+(1.-xp)*kap_tab[(band*NT+l+1)*Np+m])+
            (1.-xt)*(xp*kap_tab[(band*NT+l)*Np+m+1]+(1.-xp)*kap_tab[(band*NT+l)*Np+m]));
    
          // TR switch turn off kappa and B
          kap[ind] *= tr_switch[ind];
          B[ind]   *= tr_switch[ind];
        }
      }
    }


// *****************************************************************
// *    update diffusion-approx. boundary condition at bottom      *
// *****************************************************************
  if(isgbeg[0]==1) {
#pragma acc parallel loop collapse(5) \
 present(this[:1], z_rbuf[:Nbands*2*2*2*NMU*nx*ny], B[:nx*ny*nz])
    for(int YDIR=FWD;YDIR<=BWD;++YDIR) {
      for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR) {
        for(int l=0;l<NMU;++l) {
          for(int y=0;y<ny;++y) {
            for(int x=0;x<nx;++x) {
              //z_rbuf[band][YDIR][XDIR][UP][l][x*ny+y]=B[y*nx*nz+x*nz];
              z_rbuf[((((band*2+YDIR)*2+XDIR)*2+UP)*NMU+l)*nx*ny+(x*ny+y)]=B[y*nx*nz+x*nz];
            }
	  }
	}
      }
    }
  }

// *****************************************************************
// *    no incoming radiation on the top                           *
// *****************************************************************
  if(isgend[0]==1) {
#pragma acc parallel loop collapse(3) \
 present(this[:1], z_rbuf[:Nbands*2*2*2*NMU*nx*ny])
    for(int YDIR=FWD;YDIR<=BWD;++YDIR) {
      for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR) {
        for(int i = 0; i < nx*ny*NMU; i++) {
          //z_rbuf[band][YDIR][XDIR][DOWN][0][i] = 0.0;
          z_rbuf[((((band*2+YDIR)*2+XDIR)*2+DOWN)*NMU*nx*ny+i)] = 0.0;
	}
      }
    }
  }
// *****************************************************************
// *  solve the transfer                                           *
// *****************************************************************

  driver(DZ,DX,DY,band); 

// *****************************************************************
// *  outgoing radiative flux                                      *
// *****************************************************************

  gFr_mean[band] = 0.0;

  if(isgend[0]==1){
    double gFr_mean_reduc = 0.0;
#pragma acc parallel loop collapse(2) \
 present(this[:1], Fz[:nx*ny*nz]) \
 reduction(+:gFr_mean_reduc)
    for(int y=0; y<ny-yo; y++)
      for(int x=0; x<nx-xo; x++) {
        int ind = y*nx*nz + x*nz + nz-1;
        int y_off = yo*nx*nz, x_off = xo*nz;
        gFr_mean_reduc+=Fz[ind] +
                        Fz[ind+y_off] +
                        Fz[ind+x_off] +
                        Fz[ind+y_off+x_off];
      }
    gFr_mean[band]=gFr_mean_reduc*0.25;
  }
 
  // If I am saving 5000A intensity then wipe after last bin. if I am saving the continuum bin then wipe the second last bin.

  if (((cont_bin==1)&&(band==0))||((cont_bin==0)&&(band==1))) {
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
    for(int y=0;y<ny;++y)
      for(int x=0;x<nx;++x)
        I_o[y*nx+x] = 0.0;
  }

  }// end loop over bands

  MPI_Allreduce(&gFr_mean[0],&Fr_mean[0],Nbands,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  F_o = 0.0;
 

  for (int band=0;band<Nbands;++band){
    Fr_mean[band]/=(Grid.gsize[1]*Grid.gsize[2]);
    F_o+=Fr_mean[band];
  }
  
  // If I have the band calculate tau5000 reference grid. If I want 5000 A continuum
  // output then switch that on too.
  
  if (N5000){
    int I5000_out = 0;
    if ((cont_bin==1)&&(need_I==1))
      I5000_out = 1;
    
    get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,I_band,I5000_out);
   
    if ((cont_bin==1)&&(need_I==1)) {
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], I_band[:nx*ny])
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] +=I_band[y*nx+x];
    }
  }
  calc_Qtot_and_Tau(Grid, Run, Physics);
  
  if (Run.NeedsSlice() && Run.RT_HAVG)
    save_1D_avg(Run.path_2D,Run.globiter,Run.time); 

  PGI_COMPARE(&dt_rad, double, 1, "dt_rad", "rt.cc", "RTS::wrapper", 15) 

  return dt_rad;
}

void RTS::calc_Qtot_and_Tau(GridData &Grid, const RunData &Run, const PhysicsData &Physics){

 double tau_min = pow(Physics.rt[i_rt_tau_min],2);
//DP - Tabulated NLTE extension

 cState *U = Grid.U;

// *****************************************************************
// * Limit Qtot if dr_rad<rt_tstep, consider only region with      *
// * rho > rt_rho_min for time step computation. Compute Fout      *
// * from Qrad                                                     *
// *****************************************************************
  dt_rad=0.0;
  double _dt_rad = 0.0;
  double qsum=0.0;
  double _qsum=0.0;
  //ACCH::UpdateGPU(Grid.Tau, Grid.bufsize*sizeof(double));
  //ACCH::UpdateGPU(Grid.Jtot, Grid.bufsize*sizeof(double));
  //ACCH::UpdateGPU(Grid.Stot, Grid.bufsize*sizeof(double));
  //ACCH::UpdateGPU(U, Grid.bufsize*sizeof(cState));
  //ACCH::UpdateGPU(Grid.Qtot, Grid.bufsize*sizeof(double));
#pragma acc enter data copyin(Qt[:(ny-yo)*(nx-xo)*(nz-zo)])
#pragma acc parallel loop gang collapse(2) \
  present(this[:1], Grid[:1], Grid.Tau[:Grid.bufsize], Tau[:nx*ny*nz], \
          Grid.Jtot[:Grid.bufsize], Jt[:nx*ny*nz], Grid.Stot[:Grid.bufsize], \
	  St[:nx*ny*nz],tr_switch[:nx*ny*nz], U[:Grid.bufsize], Grid.Qtot[:Grid.bufsize],Qt[:(ny-yo)*(nx-xo)*(nz-zo)]) \
  copyin(next[1:2]) reduction(+:_qsum) reduction(max:_dt_rad) 
  for(int y = yo; y < ny; y++)
    for(int x = xo; x < nx; x++) {
      int off0 = (x+xl)*next[1]+(y+yl)*next[2];
      int xyoff = y*nx*nz + x*nz;
#pragma acc loop vector reduction(+:_qsum) reduction(max:_dt_rad)
      for(int z = zo; z < nz; z++) {
        int node = off0+z+zl;
	Grid.Tau[node] = tau(z, x, y);
	Grid.Jtot[node] = Jtot(z, x, y);
	Grid.Stot[node] = Stot(z, x, y);
	int ind = xyoff + z;
	double scale = pow(Grid.Tau[node], 2);
	scale = scale/(scale + tau_min)*tr_switch[ind];
        double Qt_step = Qt[(((y-yo)*(nx-xo)+(x-xo))*(nz-zo)+(z-zo))]*scale;
	double inv_dt = fabs(Qt_step)/U[node].e;
	_dt_rad = max(_dt_rad, inv_dt);
	Grid.Qtot[node] = Qt_step;
	_qsum += Qt_step;
      }
    }
  dt_rad = _dt_rad;
  qsum = _qsum;
/*
  ACCH::UpdateCPU(Grid.Tau, Grid.bufsize*sizeof(double));
  //ACCH::UpdateCPU(Grid.Jtot, Grid.bufsize*sizeof(double));
  //ACCH::UpdateCPU(Grid.Stot, Grid.bufsize*sizeof(double));
  //ACCH::UpdateCPU(Grid.Qtot, Grid.bufsize
  ACCH::UpdateCPU(Qt,(ny-yo)*(nx-xo)*(nz-zo)*sizeof(double));
  for(int y = 0; y < ny-yo; y++)
    for(int x = 0; x < nx-xo; x++)
      for(int z = 0; z < nz-zo; z++) {
        double Qtest=Qt[((y*(nx-xo)+x)*(nz-zo)+z)]; 
        //fprintf(stdout,"Qtest: %21.15E\n",Qtest);
        //fprintf(stdout,"tau: %21.15E\n",tau(z, x, y));
      }
   //pcast_compare(Grid.Tau, "double", Grid.bufsize, "Grid.Tau", "rt.cc", "calc_Qtot_and_Tau", 1);
   */
  exchange_single_acc(Grid,Grid.Tau);
  //pcast_compare(Grid.Tau, "double", Grid.bufsize, "Grid.Tau", "rt.cc", "calc_Qtot_and_Tau", 2);

  //ACCH::UpdateGPU(Grid.Tau, Grid.bufsize*sizeof(double));

  double Fqrad;

  MPI_Allreduce(&qsum,&Fqrad,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  Fqrad=-Fqrad*Grid.dx[0]/(Grid.gsize[1]*Grid.gsize[2]);

  if(myrank==0) fprintf(stdout,"RT energy flux: %21.15E %21.15E\n",F_o,Fqrad);
  
  F_o=Fqrad;
  dt_rad=Physics.rt[i_rt_cfl]/dt_rad;
  
}

void RTS::integrate(
  const double c[4]
)
{

  if(NDIM == 3) {
    int off0 = iystep[0]*nx*nz+ixstep[0]*nz+izstep[0];
    int off1 = iystep[1]*nx*nz+ixstep[1]*nz+izstep[1];
    int off2 = iystep[2]*nx*nz+ixstep[2]*nz+izstep[2];
    int off3 = iystep[3]*nx*nz+ixstep[3]*nz+izstep[3];

    const double c0 = c[0];
    const double c1 = c[1];
    const double c2 = c[2];
    const double c3 = c[3];
    const int size = nx*ny*nz;
    int constoff = 1;
    int str11,str21,str31,inustr11,inustr21,inustr31;
    static int szIn1, szIn2;
    static int szCf1, szCf2;
    bool ix, iy, iz;
  ix = (ixstep[0] == 1 && ixstep[1] == 1    &&
        ixstep[2] == 1 && ixstep[3] == 1)   ||
       (ixstep[0] == -1 && ixstep[1] == -1  &&
        ixstep[2] == -1 && ixstep[3] == -1);
  iy = (iystep[0] == 1 && iystep[1] == 1    &&
        iystep[2] == 1 && iystep[3] == 1)   ||
       (iystep[0] == -1 && iystep[1] == -1  &&
        iystep[2] == -1 && iystep[3] == -1);
 iz = (izstep[0] == 1 && izstep[1] == 1    &&
       izstep[2] == 1 && izstep[3] == 1)   ||
      (izstep[0] == -1 && izstep[1] == -1  &&
       izstep[2] == -1 && izstep[3] == -1);

    if(ix){
      constoff = ixstep[0]*ny*nz;
      off0 = constoff+iystep[0]*nz+izstep[0];
      off1 = constoff+iystep[1]*nz+izstep[1];
      off2 = constoff+iystep[2]*nz+izstep[2];
      off3 = constoff+iystep[3]*nz+izstep[3];
      
      str11=ny*nz; str21=nz; str31=1;
      inustr11 = (ny-1)*(nz-1); inustr21=(nz-1); inustr31=1; 
    }else if(iy){
      constoff = iystep[0]*nx*nz;
      off0 = constoff+ixstep[0]*nz+izstep[0];
      off1 = constoff+ixstep[1]*nz+izstep[1];
      off2 = constoff+ixstep[2]*nz+izstep[2];
      off3 = constoff+ixstep[3]*nz+izstep[3];

      str11=nx*nz;  str21=nz; str31=1;
      inustr11 = (nx-1)*(nz-1);inustr21=(nz-1); inustr31=1;
    }else {
      constoff = izstep[0]*nx*ny;
      off0 = constoff+iystep[0]*nx+ixstep[0];
      off1 = constoff+iystep[1]*nx+ixstep[1];
      off2 = constoff+iystep[2]*nx+ixstep[2];
      off3 = constoff+iystep[3]*nx+ixstep[3];

      str11=nx*ny;  str21=nx; str31=1;
      inustr11 = (nx-1)*(ny-1);inustr21=(nx-1); inustr31=1;

    } 
    double * dI_n   = (double*) ACCH::GetDevicePtr(I_n);
    double * dcoeff = (double*) ACCH::GetDevicePtr(coeff);
    double * dIn1   = (double*) ACCH::GetDevicePtr(I_n1);
    double * dCf1   = (double*) ACCH::GetDevicePtr(coeff1);
    double * dCf2   = (double*) ACCH::GetDevicePtr(coeff2);


if (ix || iz) {

#pragma acc parallel loop gang vector tile(1,32, 32) async(1) \
  deviceptr(dI_n, dIn1) 
    for(int b1 = 0; b1 < bound1+1; b1++) { //Z
      for(int b2 = 0; b2 < bound2+1; b2++) { //Y
       for(int b3 = 0; b3 < bound3+1; b3++) { //X
          int b1i = (b1_i-step1) + (b1*step1); //z_i+z*zstep
          int b2i = (b2_i-step2) + (b2*step2); //Y_i+y*ystep
          int b3i = (b3_i-step3) + (b3*step3); //x_i +x*xstep
          
          int ind = b1i*str1 + b2i*str2 + b3i*str3;          // (zi*1)+(yi*nx*nz)       +(xi*nz)

          int ind1 = b1i*str11 + b2i*str21 + b3i*str31;            // (zi*nx*ny)       +(yi*nx)     +(xi*1)
        
          dIn1[ind1] = dI_n[ind];
        }
       }
      }


#pragma acc parallel loop gang vector tile(1,32, 32) async(2) \
  deviceptr(dcoeff,dCf1,dCf2) 
    for(int b1 = 0; b1 < bound1; b1++) { //Z
      for(int b2 = 0; b2 < bound2; b2++) { //Y
       for(int b3 = 0; b3 < bound3; b3++) { //X
          int i_nu_acc = b1*inustr1+ b2*inustr2 + b3*inustr3; //(zi*1)+(yi*(nx-1)(nz-1))+(xi*(nz-1))
          int i_nu_acc1 = b1*inustr11 + b2*inustr21 + b3*inustr31; // (zi*(nx-1)*(ny-1))+(yi*(nx-1))+(xi*1)
          dCf1[i_nu_acc1] = dcoeff[i_nu_acc];
          dCf2[i_nu_acc1] = dcoeff[i_nu_acc+size];
        }
       }
      }

#pragma acc wait
#pragma acc loop seq independent
    for(int b1 = 0; b1 < bound1; b1++) {
      int b1i = b1_i + b1*step1;
      int b1str = b1i * str11;
      int b1inustr = b1*inustr11;
#pragma acc parallel loop collapse(2) async independent deviceptr(dIn1,dCf1,dCf2)
      for(int b2 = 0; b2 < bound2; b2++) {
        for(int b3 = 0; b3 < bound3; b3=b3++) {
#pragma acc cache(dIn1[b1str-constoff:str21],dCf1[b1inustr:inustr21],dCf2[b1inustr:inustr21])
          int b2i = b2_i + b2*step2;
          int b3i = b3_i + b3*step3;
         
          int ind =  b1str + b2i*str21 + b3i*str31;

          int i_nu_acc  = b1inustr + b2*inustr21 +b3*inustr31;

          double I_upw = c0*dIn1[ind-off0] + 
                         c1*dIn1[ind-off1] +
                         c2*dIn1[ind-off2] +
                         c3*dIn1[ind-off3];
 
          dIn1[ind]= I_upw*dCf1[i_nu_acc]+dCf2[i_nu_acc];
        }
      }
    }
#pragma acc wait
#pragma acc parallel loop gang vector tile(1,32, 32) deviceptr(dI_n,dIn1) 
  for(int b1 = 0; b1 < bound1+1; b1++) { //Z
      for(int b2 = 0; b2 < bound2+1; b2++) { //Y
       for(int b3 = 0; b3 < bound3+1; b3++) { //X
          int b1i = (b1_i-step1) + (b1*step1); //z_i+z*zstep
          int b2i = (b2_i-step2) + (b2*step2); //Y_i+y*ystep
          int b3i = (b3_i-step3) + (b3*step3); //x_i +x*xstep
          int ind = b1i*str1 + b2i*str2 + b3i*str3;
          int ind1 = b1i*str11 + b2i*str21 + b3i*str31;
          dI_n[ind] = dIn1[ind1];
        }
       }
      }
//#pragma acc wait
  } else { //end if (ix || iz)
 
#pragma acc loop seq 
    for(int b1 = 0; b1 < bound1; b1++) {
      int b1i = b1_i + b1*step1;
      int b1off = b1i*str1;
      int b1_i_nu = b1*inustr1;
#pragma acc parallel loop gang async \
 deviceptr(dI_n, dcoeff)
      for(int b2 = 0; b2 < bound2; b2++) {
        int b2i = b2_i + b2*step2;
        int b2off = b1off + b2i*str2;
        int b2_i_nu = b1_i_nu + b2*inustr2;
#pragma acc loop vector
        for(int b3 = 0; b3 < bound3; b3++) {
          int b3i = b3_i + b3*step3;
          int ind = b2off + b3i*str3;
          int i_nu_acc = b2_i_nu + b3*inustr3;
          double I_upw = c0*dI_n[ind-off0] +
                         c1*dI_n[ind-off1] +
                         c2*dI_n[ind-off2] +
                         c3*dI_n[ind-off3];
          dI_n[ind]=I_upw*dcoeff[i_nu_acc]+dcoeff[size+i_nu_acc];
        }
      }
    }
#pragma acc wait

 }
//#pragma acc exit data delete(dIn1,dCf1,dCf2)
  } // end DIM 3
/*
  if(NDIM==2){
    for(int i=0;i<4;i++) off[i]=ixstep[i]*nz+izstep[i];
    for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep) {
      int xoff=(xi-xl)*stride[1];
      for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep) {
        int ind=xoff+zi;
        double I_upw=c[0]*ii[ind-off[0]]+c[1]*ii[ind-off[1]]+c[2]*ii[ind-off[2]]+c[3]*ii[ind-off[3]];
        ii[ind]=I_upw*coeff[i_nu]+coeff[nx*ny*nz+i_nu];
        i_nu+=1;
      }
    }
  }

  if(NDIM==1) {
    for(int i=0;i<4;i++) off[i]=izstep[i];
    for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep) {
      double I_upw=c[0]*ii[zi-off[0]]+c[1]*ii[zi-off[1]]+c[2]*ii[zi-off[2]]+c[3]*ii[zi-off[3]];
      ii[zi]=I_upw*coeff[i_nu]+coeff[nx*ny*nz+i_nu];
      i_nu+=1;
    }
  }
  */
}

void RTS::driver(double DZ, double DX, double DY, int band){
  double etime=0.0,atime=0.0,cmp_time1=0.0,cmp_time2=0.0,buf_time=0.0,err_time=0.0,flx_time=0.0,tau_time=0.0; 
  double ttime=MPI_Wtime();
  
  const int stride[2]={nx*nz,nz}; 
  int stepvec[3][4][3] = { {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},
               {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
               {{0,0,1},{0,1,1},{1,0,1},{1,1,1}} };

#pragma acc enter data copyin(this[:1], I_n[:nx*ny*nz], Fz[:nx*ny*nz], Fx[:nx*ny*nz], \
         Fy[:nx*ny*nz], J_band[:nx*ny*nz],rho[:nx*ny*nz],kap[:nx*ny*nz], coeff[:nx*ny*nz*2],B[:nx*ny*nz],sbuf[:ny][:nx], rbuf[:ny][:nx], \
         tI_n[:nx*ny*nz], trho[:nx*ny*nz], tkap[:nx*ny*nz], tB[:nx*ny*nz],Tau[:nx*ny*nz]) //Qt[:(ny-yo)*(nx-xo)*(nz-zo)])

#pragma acc parallel loop \
 present(this[:1], I_n[:nx*ny*nz], Fz[:nx*ny*nz], Fx[:nx*ny*nz], \
         Fy[:nx*ny*nz], J_band[:nx*ny*nz])
  for(int i = 0; i < nx*ny*nz; i++) {
    I_n[i] = 0.0;
    tI_n[i] = 0.0;
    Fz[i] = 0.0;
    Fx[i] = 0.0;
    Fy[i] = 0.0;
    J_band[i] = 0.0;
  }
  double maxerr_up=0.0,maxerr_down=0.0;
// loop over octants & determination of loop direction
  double itavg=0.0;
  double aravg=0.0;

  Transpose_Rho_Kap_B(rho, trho, kap, tkap, B, tB, nx, ny, nz);

// main loop
  for(int ZDIR=UP;ZDIR<=DOWN;++ZDIR){
    int zi_i=(ZDIR==UP)?zl+1:zh-1,zi_f=(ZDIR==UP)?zh:zl,zstep=(ZDIR==UP)?1:-1;
    for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR){
      int xi_i=(XDIR==RIGHT)?xl+1:xh-1,xi_f=(XDIR==RIGHT)?xh:xl,xstep=(XDIR==RIGHT)?1:-1;
      for(int YDIR=FWD;YDIR<=BWD;++YDIR){
        int yi_i=(YDIR==FWD)?yl+1:yh-1,yi_f=(YDIR==FWD)?yh:yl,ystep=(YDIR==FWD)?1:-1;
          for(int l=0;l<NMU;++l){
            double I_min=max(1.0,threshold*Fr_mean[band]/(NMU*pow(2,NDIM)));
            double c[]={a_00[ibase[l]][l],a_01[ibase[l]][l],a_10[ibase[l]][l],a_11[ibase[l]][l]};
        for(int m=0;m<=3;++m){
          ixstep[m]=stepvec[ibase[l]][m][0]*xstep;
          iystep[m]=stepvec[ibase[l]][m][1]*ystep;
          izstep[m]=stepvec[ibase[l]][m][2]*zstep;
        }
//#pragma acc enter data copyin(c[:4])
        //IntegrateSetup(yi_i-yl, xi_i-xl, zi_i-zl, ystep, xstep, zstep);

        double stime=MPI_Wtime();
        //interpol(zi_i,zi_f,zstep,xi_i,xi_f,xstep,yi_i,yi_f,ystep,l,coeff,B);
        Transpose_interpol(
          rho, trho, kap, tkap, B, tB, coeff,
          nx, ny, nz, XDIR, YDIR, ZDIR, l, ixstep, iystep, izstep,
          xi_i-xl, yi_i-yl, zi_i-zl, c, ds_upw
        );
        cmp_time1+=MPI_Wtime()-stime;
        stime=MPI_Wtime(); 
            int rt_iter=0;
            int rt_min_iter=(call_count<2)?0:numits[band][YDIR][XDIR][ZDIR][l]-(!(call_count%3));            
        double gmaxerr=1.0E10;
// Iteration start
        while(gmaxerr>=threshold){
          rt_iter+=1;
          itavg+=1.0;
// read BC from recvbuf
          stime = MPI_Wtime();
          //readbuf(band,l,ZDIR,XDIR,YDIR);
          Transpose_readbuf(
            x_sbuf, x_rbuf, x_oldbuf,
            y_sbuf, y_rbuf, y_oldbuf,
            z_sbuf, z_rbuf, z_oldbuf,
            I_n, tI_n, nx, ny, nz,
            XDIR, YDIR, ZDIR, l, band,
            ixstep, iystep, izstep
          );
          buf_time += MPI_Wtime()-stime;
// loop over the grid points
          stime=MPI_Wtime();
          //integrate(c);
          Transpose_integrate(
            I_n, tI_n, coeff, nx, ny, nz,
            XDIR, YDIR, ZDIR, ixstep, iystep, izstep,
            xi_i-xl, yi_i-yl, zi_i-zl, c
          );
	  cmp_time2 += MPI_Wtime()-stime;
          // write new BC, store old BC in oldbuf
          stime = MPI_Wtime();
          //writebuf(band,l,ZDIR,XDIR,YDIR); 
          Transpose_writebuf(
            x_sbuf, y_sbuf, z_sbuf, I_n, tI_n, nx, ny, nz, XDIR, YDIR, ZDIR, l, band,
            ixstep, iystep, izstep
          );
          buf_time += MPI_Wtime()-stime;
          stime = MPI_Wtime();
          exchange(band, l, ZDIR, XDIR,YDIR);
          etime += MPI_Wtime()-stime;
          if(rt_iter>=rt_min_iter){
            stime = MPI_Wtime();
            double err=error(band,l,ZDIR,XDIR,YDIR,I_min);
            err_time += MPI_Wtime()-stime;
            stime = MPI_Wtime();
            MPI_Allreduce(&err,&gmaxerr,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            aravg += 1.0;
            atime += MPI_Wtime()-stime;
          }
        } // end main loop
        if(ZDIR==UP){
          maxerr_up=max(maxerr_up,gmaxerr);
        }else{
          maxerr_down=max(maxerr_down,gmaxerr);
        }
        numits[band][YDIR][XDIR][ZDIR][l]=rt_iter;
// read BC from recvbuf 
        stime=MPI_Wtime();
        //readbuf(band,l,ZDIR,XDIR,YDIR);
        Transpose_readbuf(
          x_sbuf, x_rbuf, x_oldbuf,
          y_sbuf, y_rbuf, y_oldbuf,
          z_sbuf, z_rbuf, z_oldbuf,
          I_n, tI_n, nx, ny, nz,
          XDIR, YDIR, ZDIR, l, band,
          ixstep, iystep, izstep
        );
        buf_time += MPI_Wtime()-stime;
// increment for F and J 
        stime=MPI_Wtime();
        //flux(l,ZDIR,XDIR,YDIR);
        Transpose_flux(
          I_n, tI_n, Fx, Fy, Fz, J_band,
          wmu, xmu, nx, ny, nz, XDIR, YDIR, ZDIR, l,
          ixstep, iystep, izstep
        );
        
        flx_time += MPI_Wtime()-stime;

      }// end l
      }// end YDIR
    }// end XDIR
  }// end ZDIR 
  double stime=MPI_Wtime();
  tauscale_qrad(band,DX,DY,DZ,B);
  tau_time+=MPI_Wtime()-stime;
  if((myrank==0) && (verbose>1)){
    fprintf(stdout,"rt_driver iter : %f %f \n",aravg/(8.0*NMU),itavg/(8.0*NMU));
    fprintf(stdout,"rt_driver error: %e %e \n",maxerr_up,maxerr_down);
    /*
    int oct=0;
    for(int i1=0;i1<=1;i1++)
      for(int i2=0;i2<=1;i2++)
	for(int i3=0;i3<=1;i3++){
	  oct +=1;
	  for(int l=0;l<NMU;++l)
	    fprintf(stdout,"rt_driver iter : %d %d %d \n",oct,l,numits[0][i1][i2][i3][l]);
	}
    */
  }

  ttime=MPI_Wtime()-ttime;  
  if((myrank==0) && (verbose>2))
    fprintf(stdout,"rt_driver time : %f %f %f %f %f %f %f %f %f %f \n",ttime,cmp_time1,
       cmp_time2,buf_time,err_time,flx_time,tau_time,etime,atime,(etime+atime)/ttime);
   

  call_count+=1;
}

void RTS::interpol_and_integrate(
  const double c[4], double * Ss, int l
)
{
  double ds3=ds_upw[l]*inv3,ds6=ds_upw[l]*inv6;

  const int off0 = iystep[0]*nx*nz+ixstep[0]*nz+izstep[0];
  const int off1 = iystep[1]*nx*nz+ixstep[1]*nz+izstep[1];
  const int off2 = iystep[2]*nx*nz+ixstep[2]*nz+izstep[2];
  const int off3 = iystep[3]*nx*nz+ixstep[3]*nz+izstep[3];

  double * dI_n = (double*) ACCH::GetDevicePtr(I_n);
  double * drho = (double*) ACCH::GetDevicePtr(rho);
  double * dkap = (double*) ACCH::GetDevicePtr(kap);
  double * dSs  = (double*) ACCH::GetDevicePtr(Ss);

  const double c0 = c[0];
  const double c1 = c[1];
  const double c2 = c[2];
  const double c3 = c[3];
 
#pragma acc loop seq
  for(int b1 = 0; b1 < bound1; b1++) {
    int b1i = b1_i + b1*step1;
    int b1off = b1i*str1;
#pragma acc parallel loop gang async \
 deviceptr(dI_n, drho, dkap, dSs)
    for(int b2 = 0; b2 < bound2; b2++) {
      int b2i = b2_i + b2*step2;
      int b2off = b1off + b2i*str2;
#pragma acc loop vector
      for(int b3 = 0; b3 < bound3; b3++) {
        int b3i = b3_i + b3*step3;
        int ind = b2off + b3i*str3;
        double r_upw = c0*drho[ind-off0] + c1*drho[ind-off1] +
                       c2*drho[ind-off2] + c3*drho[ind-off3];
        double k_upw = c0*dkap[ind-off0] + c1*dkap[ind-off1]+
                       c2*dkap[ind-off2] + c3*dkap[ind-off3];
        double S_upw = c0*dSs[ind-off0]  + c1*dSs[ind-off1]+
                       c2*dSs[ind-off2]  + c3*dSs[ind-off3];

        double r0 = drho[ind];
        double k0 = dkap[ind];
        double S0 = dSs[ind];

	double I_upw = c0*dI_n[ind-off0] +
                       c1*dI_n[ind-off1] +
                       c2*dI_n[ind-off2] +
                       c3*dI_n[ind-off3];
 
        double expo, source, dt, w0, w1;
        dt = ds3*(k_upw*r_upw + k0*r0) + ds6*(k0*r_upw+k_upw*r0);
        if(dt > dtau_min2) {
          expo = exp(-dt);
          if (dt > dtau_min){
            w0 = 1.0-expo;
            w1 = w0-dt*expo;
          }else{
            w0 = dt-dt*dt/2.0+dt*dt*dt/6.0;
            w1 = dt*dt/2.0-dt*dt*dt/3.0;
          }
          source=S0*(w0-w1/dt)+S_upw*(w1/dt);
	} else {
          expo = 1.0;
	  source = 0.0;
	}
        dI_n[ind]=I_upw*expo+source;
      }
    }
  }
#pragma acc wait
}

void RTS::interpol(int zi_i,int zi_f,int zstep,int xi_i,int xi_f,int xstep,
           int yi_i,int yi_f,int ystep,int l,double* coeff, double * B)
{

  double ds3=ds_upw[l]*inv3,ds6=ds_upw[l]*inv6;
  double c[]={a_00[ibase[l]][l],a_01[ibase[l]][l],a_10[ibase[l]][l],a_11[ibase[l]][l]};
  double c0 = a_00[ibase[l]][l];
  double c1 = a_01[ibase[l]][l];
  double c2 = a_10[ibase[l]][l];
  double c3 = a_11[ibase[l]][l];
  int zmin=(zi_i<zi_f)?zi_i:zi_f;
  int zmax=(zi_i>zi_f)?zi_i:zi_f;

  double r_upw[zmax+1],k_upw[zmax+1],S_upw[zmax+1],r0[zmax+1],k0[zmax+1],S0[zmax+1];
 // double _r_upw, _k_upw, _S_upw, _r0, _k0, _S0;
  int i_nu=0;

  const int stride[] = {nx*nz, nz, 1};
  const int stride0 = nx*nz;
  const int stride1 = nz;
  const int stride2 = 1;
  const int stride_inter[] = {(nx-1)*(nz-1), nz-1, 1};
  const int stride_inter0 = (nx-1)*(nz-1);
  const int stride_inter1 = nz-1;
  const int stride_inter2 = 1;
  const int off[] = {
    iystep[0]*stride[0]+ixstep[0]*stride[1]+izstep[0],
    iystep[1]*stride[0]+ixstep[1]*stride[1]+izstep[1],
    iystep[2]*stride[0]+ixstep[2]*stride[1]+izstep[2],
    iystep[3]*stride[0]+ixstep[3]*stride[1]+izstep[3]
  };
  const int off0 = iystep[0]*stride[0]+ixstep[0]*stride[1]+izstep[0];
  const int off1 = iystep[1]*stride[0]+ixstep[1]*stride[1]+izstep[1];
  const int off2 = iystep[2]*stride[0]+ixstep[2]*stride[1]+izstep[2];
  const int off3 = iystep[3]*stride[0]+ixstep[3]*stride[1]+izstep[3];
  const int size = nx*ny*nz;

  if(NDIM==3){
#pragma acc parallel loop collapse(3) \
 present(this[:1], rho[:nx*ny*nz], kap[:nx*ny*nz], B[:nx*ny*nz], coeff[:nx*ny*nz*2]) 
    for(int y = 0; y < ny-1; y++) {
      for(int x = 0; x < nx-1; x++) {
        for(int z = 0; z < nz-1; z++) {
          double _r_upw, _k_upw, _S_upw, _r0, _k0, _S0;
          int yi = (yi_i-yl) + y*ystep;
          int xi = (xi_i-xl) + x*xstep;
          int zi = (zi_i-zl) + z*zstep;
          int ind = yi*stride0 + xi*stride1 + zi;
	  int i_nu_ = y*stride_inter0 + x*stride_inter1 + z;
          _r_upw=
            c0*rho[ind-off0]+
            c1*rho[ind-off1]+
            c2*rho[ind-off2]+
            c3*rho[ind-off3];
          
          _k_upw=
            c0*kap[ind-off0]+
            c1*kap[ind-off1]+
            c2*kap[ind-off2]+
            c3*kap[ind-off3];

          _S_upw=
            c0*B[ind-off0]+
            c1*B[ind-off1]+
            c2*B[ind-off2]+
            c3*B[ind-off3];

          _r0=rho[ind];
          _k0=kap[ind];
          _S0=B[ind];

          double dt=ds3*(_k_upw*_r_upw+_k0*_r0)+ds6*(_k0*_r_upw+_k_upw*_r0);
          double expo=exp(-dt);
          double w0,w1;
          if (dt > dtau_min){
            w0=1.0-expo;
            w1=w0-dt*expo;
          }else{
            w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
            w1=dt*dt/2.0-dt*dt*dt/3.0;
          }
          double source=_S0*(w0-w1/dt)+_S_upw*(w1/dt);

          if (dt > dtau_min2){
            coeff[i_nu_] = expo;
            coeff[i_nu_+size] = source;
          }else{
            coeff[i_nu_] = 1.0; 
            coeff[i_nu_+size] = 0.0;
          }
        }
      }
    }
  }

  if(NDIM==1){
    for(int zi=zmin;zi<=zmax;++zi){
      r_upw[zi]=
        c[0]*rho[(zi-izstep[0]-zl)]+
        c[1]*rho[(zi-izstep[1]-zl)]+
        c[2]*rho[(zi-izstep[2]-zl)]+
        c[3]*rho[(zi-izstep[3]-zl)];

      k_upw[zi]=
        c[0]*kap[zi-izstep[0]-zl]+
        c[1]*kap[zi-izstep[1]-zl]+
        c[2]*kap[zi-izstep[2]-zl]+
        c[3]*kap[zi-izstep[3]-zl];

      S_upw[zi]=
        c[0]*B[zi-izstep[0]-zl]+
        c[1]*B[zi-izstep[1]-zl]+
        c[2]*B[zi-izstep[2]-zl]+
        c[3]*B[zi-izstep[3]-zl];

      r0[zi]=rho[zi-zl];
      k0[zi]=kap[zi-zl];
      S0[zi]=B[zi-zl];
    }

    for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
      double dt=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
      double expo=exp(-dt);
      double w0,w1;
      if (dt > dtau_min){
        w0=1.0-expo;
        w1=w0-dt*expo;
      }else{
        w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
        w1=dt*dt/2.0-dt*dt*dt/3.0;
      }
      double source=S0[zi]*(w0-w1/dt)+S_upw[zi]*(w1/dt);

      if (dt > dtau_min2){
        coeff[i_nu] = expo;
        coeff[i_nu+size] = source;
      }else{
        coeff[i_nu] = 1.0; 
        coeff[i_nu+size] = 0.0;
      }
    i_nu+=1;
    }
  }

  if(NDIM==2){
    for(int xi=xi_i;xi!=xi_f+xstep;xi=xi+xstep){
      for(int zi=zmin;zi<=zmax;++zi){
        r_upw[zi]=
          c[0]*rho[(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
          c[1]*rho[(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
          c[2]*rho[(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
          c[3]*rho[(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

        k_upw[zi]=
          c[0]*kap[(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
          c[1]*kap[(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
          c[2]*kap[(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
          c[3]*kap[(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

        S_upw[zi]=
          c[0]*B[(xi-ixstep[0]-xl)*nz+(zi-izstep[0]-zl)]+
          c[1]*B[(xi-ixstep[1]-xl)*nz+(zi-izstep[1]-zl)]+
          c[2]*B[(xi-ixstep[2]-xl)*nz+(zi-izstep[2]-zl)]+
          c[3]*B[(xi-ixstep[3]-xl)*nz+(zi-izstep[3]-zl)];

        r0[zi]=rho[(xi-xl)*nz+(zi-zl)];
        k0[zi]=kap[(xi-xl)*nz+(zi-zl)];
        S0[zi]=B[(xi-xl)*nz+(zi-zl)];
      }
      for(int zi=zi_i;zi!=zi_f+zstep;zi=zi+zstep){
        double dt=ds3*(k_upw[zi]*r_upw[zi]+k0[zi]*r0[zi])+ds6*(k0[zi]*r_upw[zi]+k_upw[zi]*r0[zi]);
        double expo=exp(-dt);
        double w0,w1;
        if (dt > dtau_min){
          w0=1.0-expo;
          w1=w0-dt*expo;
        }else{
          w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
          w1=dt*dt/2.0-dt*dt*dt/3.0;
        }
        double source=S0[zi]*(w0-w1/dt)+S_upw[zi]*(w1/dt);

        if (dt > dtau_min2){
          coeff[i_nu] = expo;
          coeff[i_nu+size] = source;
        }else{
          coeff[i_nu] = 1.0; 
          coeff[i_nu+size] = 0.0;
        }
      i_nu+=1;
      }
    }
  }
}

double RTS::error(int band,int l,int ZDIR,int XDIR,int YDIR,double I_min)
{
  int z0=(ZDIR==UP),z1=(ZDIR==DOWN); 
  int x0=(XDIR==RIGHT),x1=(XDIR==LEFT); 
  double err_max=0.0;
  double _err_max = err_max;
  int adrX = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nz);
  int adrY = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*nx*nz);
  int adrZ = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nx);
   
   real *ysb=&y_sbuf[adrY],*yob=&y_oldbuf[adrY];
   real *xsb=&x_sbuf[adrX],*xob=&x_oldbuf[adrX];
   real *zsb=&z_sbuf[adrZ],*zob=&z_oldbuf[adrZ];
//#pragma acc enter data copyin(err_max)
#pragma acc data present(this[:1], ysb[:nx*nz], yob[:nx*nz],xsb[:ny*nz], xob[:ny*nz],zsb[:ny*nx], zob[:ny*nx]) 
{
  if (NDIM==3){
#pragma acc parallel loop collapse(2) reduction(max:_err_max) \
 present(this[:1], ysb[:nx*nz], yob[:nx*nz]) 
    for(int z=z0;z<nz-z1;++z)
      for(int x=x0;x<nx-x1;++x){
        int ind=z*nx+x;
        _err_max=max(_err_max,fabs(((double) (ysb[ind]-yob[ind]))/max(I_min,(double) yob[ind])));
     }

  }
  if (NDIM>1){
#pragma acc parallel loop collapse(2) reduction(max:_err_max) \
 present(this[:1], xsb[:ny*nz], xob[:ny*nz]) 
    for(int z=z0;z<nz-z1;++z)
      for(int y=0;y<ny;++y){
        int ind=z*ny+y;
        _err_max=max(_err_max,fabs(((double) (xsb[ind]-xob[ind]))/max( I_min,(double)xob[ind])));
      }
//#pragma acc update host(err_max)
  }

#pragma acc parallel loop reduction(max:_err_max) \
 present(this[:1], zsb[:ny*nx], zob[:ny*nx]) 
  for(int ind=0;ind<nx*ny;++ind){
    _err_max=max(_err_max,fabs(((double) (zsb[ind]-zob[ind]))/max(I_min,(double)zob[ind])));
  }
//} //data
  err_max = _err_max;
  return err_max;
} //data
}

void RTS::readbuf(int band,int l,int ZDIR,int XDIR,int YDIR)
{//RHC
    real * ysb = &y_sbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * yrb = &y_rbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * yob = &y_oldbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * xsb = &x_sbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * xrb = &x_rbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * xob = &x_oldbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * zsb = &z_sbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * zrb = &z_rbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
    real * zob = &z_oldbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
#pragma acc enter data copyin(this[:1], I_n[:nx*ny*nz], yrb[:nx*nz],yob[:nx*nz], ysb[:nx*nz],xrb[:ny*nz],xob[:ny*nz], xsb[:ny*nz],zrb[:ny*nx],zob[:ny*nx], zsb[:ny*nx])
  if(NDIM==3){
    int y0=(YDIR==FWD)?0:ny-1;
    int y = y0*nx*nz;

#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], yrb[:nx*nz]) 
    for(int x=0;x<nx;++x)
      for(int z=0;z<nz;++z){
        int ob  = x*nz+z;
        I_n[y+x*nz+z]= (double) yrb[z*nx+x];
        yob[ob] = ysb[ob];
       }
  }

  if(NDIM>1){
    int x0=(XDIR==RIGHT)?0:nx-1;
    int x = x0*nz;
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], xrb[:ny*nz])
    for(int y=0;y<ny;++y)
      for(int z=0;z<nz;++z){
        int ob  = y*nz+z;
        I_n[y*nx*nz+x+z]= (double) xrb[z*ny+y];
        xob[ob] = xsb[ob];
        }
  }
  
  int z0=(ZDIR==UP)?0:nz-1;
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], zrb[:ny*nx]) 
  for(int y=0;y<ny;++y)
    for(int x=0;x<nx;++x){
      int ob  = y*nx+x;
      I_n[y*nx*nz+x*nz+z0]= (double) zrb[x*ny+y];
      zob[ob] = zsb[ob]; 
    }

#pragma acc wait
}

void RTS::writebuf(int band, int l,int ZDIR,int XDIR,int YDIR){
  if (NDIM==3){
    int y0=(YDIR==FWD)?ny-1:0;
    int y = y0*nx*nz;
    real * ysb = &y_sbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], ysb[:nx*nz]) async
    for(int x=0;x<nx;++x)
      for(int z=0;z<nz;++z)
        ysb[z*nx+x]=(real) I_n[y+x*nz+z];
  }

  if (NDIM>1){
    int x0=(XDIR==RIGHT)?nx-1:0;
    int x = x0*nz;
    real * xsb = &x_sbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], xsb[:ny*nz]) async
    for(int y=0;y<ny;++y)
      for(int z=0;z<nz;++z)
        xsb[z*ny+y]=(real) I_n[y*nx*nz+x+z];
  }

  int z0=(ZDIR==UP)?nz-1:0;
  real * zsb = &z_sbuf[((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)];
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], zsb[:ny*nx]) async
  for(int y=0;y<ny;++y)
    for(int x=0;x<nx;++x)
      zsb[x*ny+y]=(real) I_n[y*nx*nz+x*nz+z0];
#pragma acc wait
}

void RTS::exchange(int band,int l,int ZDIR,int XDIR,int YDIR)
{
  MPI_Status s1[2],s2[2],s3[2];
  MPI_Request r1[2],r2[2],r3[2];
  int tag1=1,tag2=2,tag3=3;

  // Update v. z_sbuf an den Endpkten, damit der Eckpkt-Wert vom diagonalen Nachbarn uebertragen wird
  int x0=0,y0=0,z0=0;
  int dest_rk;
  int source_rk;
  real tempxsb[nz];
  real tempzsb[nx];

  int adrX = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nz);
  int adrY = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*nx*nz);
  int adrZ = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nx);
  
  real * ysb = &y_sbuf[adrY];
  real * yrb = &y_rbuf[adrY];
  real * xsb = &x_sbuf[adrX];
  real * xrb = &x_rbuf[adrX];
  real * zsb = &z_sbuf[adrZ];
  real * zrb = &z_rbuf[adrZ];


#pragma acc data present(this[:1],yrb[:nx*nz], ysb[:nx*nz],xrb[:ny*nz], xsb[:ny*nz],zrb[:ny*nx],zsb[:ny*nx]) 
{
  if (NDIM==3){
    // y-direction
    dest_rk=(YDIR==FWD)?rightr[2]:leftr[2];
    source_rk=(YDIR==FWD)?leftr[2]:rightr[2];

    z0=(ZDIR==UP)?nz-1:0;
    x0=(XDIR==RIGHT)?nx-1:0;
    y0=(YDIR==FWD)?0:ny-1;
#pragma acc host_data use_device(yrb,ysb)
{
    MPI_Irecv(yrb,nx*nz,REALTYPE,source_rk,tag2,MPI_COMM_WORLD,r2+0);
    MPI_Isend(ysb,nx*nz,REALTYPE,dest_rk,tag2,MPI_COMM_WORLD,r2+1);

    MPI_Waitall(2,r2,s2);
}  

#pragma acc parallel loop collapse(2) present(this[:1],yrb[:nx*nz], xsb[:ny*nz],\
   zsb[:ny*nx])  
    for(int x=0;x<nx;++x){
      for(int z=0;z<nz;++z){
        xsb[z*ny+y0]=yrb[z*nx+x0];
        zsb[x*ny+y0]=yrb[z0*nx+x];
      }
    }
  
  }

  if (NDIM >1){
    // x-direction
    dest_rk=(XDIR==RIGHT)?rightr[1]:leftr[1];
    source_rk=(XDIR==RIGHT)?leftr[1]:rightr[1];
#pragma acc host_data use_device(xrb,xsb)
{ 
    MPI_Irecv(xrb,nz*ny,REALTYPE,source_rk,tag1,MPI_COMM_WORLD,r1+0);
    MPI_Isend(xsb,nz*ny,REALTYPE,dest_rk,tag1,MPI_COMM_WORLD,r1+1);

    MPI_Waitall(2,r1,s1);
} 
    x0=(XDIR==RIGHT)?0:nx-1;
    z0=(ZDIR==UP)?nz-1:0;
#pragma acc parallel loop present(this[:1],xrb[:ny*nz], zsb[:nx*ny])
    for(int y=0;y<ny;++y)
      zsb[x0*ny+y]=xrb[z0*ny+y];
  }
  
// z-direction
  dest_rk  =(ZDIR==UP)?rightr[0]:leftr[0];
  source_rk=(ZDIR==UP)?leftr[0]:rightr[0];
#pragma acc host_data use_device(zrb,zsb)
{ 
  MPI_Irecv(zrb,nx*ny,REALTYPE,source_rk,tag3,MPI_COMM_WORLD,r3+0);
  MPI_Isend(zsb,nx*ny,REALTYPE,dest_rk,tag3,MPI_COMM_WORLD,r3+1);

  MPI_Waitall(2,r3,s3);
 }

 } //data
}

void RTS::flux(int l,int ZDIR,int XDIR,int YDIR)
{
  double zsign=(ZDIR==UP)?1.0:-1.0,xsign=(XDIR==RIGHT)?1.0:-1.0,ysign=(YDIR==FWD)?1.0:-1.0;
  double c_J = 0.125*wmu[l];
  double c_z = 0.5*PI*zsign*wmu[l]*xmu[2][l];
  double c_x = 0.5*PI*xsign*wmu[l]*xmu[0][l];
  double c_y = 0.5*PI*ysign*wmu[l]*xmu[1][l];

#pragma acc parallel loop collapse(2) gang \
 present(this[:1], I_n[:nx*ny*nz], J_band[:nx*ny*nz], \
         Fz[:nx*ny*nz], Fy[:nx*ny*nz], Fx[:nx*ny*nz])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      int yxind = y*nx*nz + x*nz;
#pragma acc loop vector
      for(int z = 0; z < nz; z++) {
        int ind = yxind + z;
        double tmp = I_n[ind];
        J_band[ind] += c_J*tmp;
        Fz[ind]     += c_z*tmp;
        Fx[ind]     += c_x*tmp;
        Fy[ind]     += c_y*tmp;
      }
    }
  }

}

void RTS::get_Tau_and_Iout(GridData &Grid, const RunData &Run, const PhysicsData &Physics, double DZ, float * B_Iout_tab, float * kap_Iout_tab, double * I_band, int calc_int){
  
  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  
  //double ** sbuf = ACCH::Malloc2D<double>(ny, nx);
  //double ** rbuf = ACCH::Malloc2D<double>(ny, nx);

  double N = pow(2,Grid.NDIM);

  const double Temp_TR = Physics.rt[i_rt_tr_tem];
  const double Pres_TR = Physics.rt[i_rt_tr_pre];


  //ACCH::UpdateGPU(Grid.pres, Grid.bufsize*sizeof(double));
  //ACCH::UpdateGPU(Grid.temp, Grid.bufsize*sizeof(double));
  //ACCH::UpdateGPU(Grid.U, Grid.bufsize*sizeof(cState));
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Grid[:1], Grid.pres[:Grid.bufsize], Grid.temp[:Grid.bufsize], Grid.U[:Grid.bufsize], \
         rho[:nx*ny*nz], tab_T[:NT], tab_p[:Np], B_Iout_tab[:NT], kap_Iout_tab[:NT][:Np], \
	 kap[:nx*ny*nz], B[:nx*ny*nz], invT_tab[:NT], invP_tab[:Np]) \
 copyin(next[1:2])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      int off0 = (x+xl)*next[1]+(y+yl)*next[2];
      int off1 = (x+xl)*next[1]+(y+yo+yl)*next[2];
      int off2 = (x+xo+xl)*next[1]+(y+yl)*next[2];
      int off3 = (x+xo+yl)*next[1]+(y+yo+yl)*next[2];
#pragma acc loop vector
      for(int z = 0; z < nz; z++) {
        int inode[]={off0+z+zl,off0+z+zo+zl,off2+z+zl,off2+z+zo+zl,off1+z+zl,off1+z+zo+zl,off3+z+zl,off3+z+zo+zl};
        double lgP = 0, lgT = 0, rm = 0;
#pragma acc loop seq
        for(int l=0;l<N;++l){
          lgP+=Grid.pres[inode[l]];
          lgT+=Grid.temp[inode[l]];
          rm +=Grid.U[inode[l]].d;
        }
        lgP /= N;
        lgT /= N;
        rm  /= N;

//disbale RT if Temp > Temp_TR above the photosphere

        double pswitch  = min(max(lgP-Pres_TR,0.0),1.0);
        double tswitch  = min(max(Temp_TR-lgT,0.0),1.0);
        double trswitch = max(pswitch,tswitch);

        lgP          = log(lgP);
        lgT          = log(lgT);
        rho[y*nx*nz+x*nz+z] = rm;

        lgT = min(max(lgT, (double) tab_T[0]),(double) tab_T[NT-1]);
        lgP = min(max(lgP, (double) tab_p[0]),(double) tab_p[Np-1]);

        int l=0;
        int m=0;
        if(lgT<tab_T[0])
          l=0;
        else if(lgT>tab_T[NT-1])
          l=NT-2;
        else {
#pragma acc loop seq 
        for (int li=0; li<=NT-2; li++){
            int lflag = 0;
            if ((lgT >= tab_T[li]) && (lgT <= tab_T[li+1])){
              lflag = lflag+1;
            if(lflag==1)
              l = li;
            }
          }
        }

        if(lgP<tab_p[0])
          m=0;
        else if(lgP>tab_p[Np-1])
          m=Np-2;
        else {
#pragma acc loop seq 
        for (int mi=0; mi<=Np-2; mi++){
            int mflag = 0;
            if ((lgP >= tab_p[mi]) && (lgP <= tab_p[mi+1]))
              mflag = mflag+1;
            if(mflag==1)
              m = mi;
          }
        }

        double xt = (lgT-tab_T[l])*invT_tab[l];
        double xp = (lgP-tab_p[m])*invP_tab[m];

// Interpolate for kappa and B
     
        B[y*nx*nz+x*nz+z]=exp(xt*B_Iout_tab[l+1]+(1.-xt)*B_Iout_tab[l]);

        kap[y*nx*nz+x*nz+z]=
                    exp(xt*(xp*kap_Iout_tab[Np*(l+1)+(m+1)]+
                           (1.-xp)*kap_Iout_tab[Np*(l+1)+m])+
                           (1.-xt)*(xp*kap_Iout_tab[Np*l+(m+1)]+
                           (1.-xp)*kap_Iout_tab[(Np*l+m)]));

 // apply tr_switch
 
        B[y*nx*nz+x*nz+z]   *= trswitch;
        kap[y*nx*nz+x*nz+z] *= trswitch;

      }
    }
  }

#pragma acc parallel loop collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], rho[:nx*ny*nz], kap[:nx*ny*nz])
  for(int y = 0; y < ny; y++)
    for(int x = 0; x < nx; x++) {
      Tau[y*nx*nz+x*nz+(nz-1)]=1.0e-12;
#pragma acc loop seq
      for(int z = nz-2; z >= 0; z--) {
        double k0=kap[y*nx*nz+x*nz+z],r0=rho[y*nx*nz+x*nz+z],k_upw=kap[y*nx*nz+x*nz+(z+1)],r_upw=rho[y*nx*nz+x*nz+(z+1)];
        double delta_tau=DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6);
        Tau[y*nx*nz+x*nz+z]=Tau[y*nx*nz+x*nz+(z+1)]+delta_tau;
      }
    }

#pragma acc parallel loop collapse(2) \
 present(this[:1], sbuf[:ny][:nx], rbuf[:ny][:nx], Tau[:nx*ny*nz])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y*nx*nz+x*nz];
    }
  }
  double ctime=MPI_Wtime();
  MPI_Scan(
    ACCH::GetDevicePtr(sbuf[0]), ACCH::GetDevicePtr(rbuf[0]),
    nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
  stime+=MPI_Wtime()-ctime;
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], rbuf[:ny][:nx], Tau[:nx*ny*nz])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      double rbufyx = rbuf[y][x]-Tau[y*nx*nz+x*nz];
#pragma acc loop vector
      for(int z=0;z<nz;++z){
        Tau[y*nx*nz+x*nz+z]+=rbufyx;
      }
    }
  }

  if (calc_int){
  //  Outgoing Intensity at top (Long Characteristics)
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], rbuf[:ny][:nx], sbuf[:ny][:nx], B[:nx*ny*nz], Tau[:nx*ny*nz])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]=0.0;
      sbuf[y][x]=0.0;
      double tmp = 0.0;
#pragma acc loop vector reduction(+:tmp)
      for(int z=1;z<nz;++z){
        double Ss1 = B[y*nx*nz+x*nz+z];
        double Ss2 = B[y*nx*nz+x*nz+z-1];
        double delta_tau=Tau[y*nx*nz+x*nz+z-1]-Tau[y*nx*nz+x*nz+z];
        if(delta_tau>dtau_min){
          double edt=exp(-delta_tau);
          double c1=(1.0-edt)/delta_tau;
          tmp+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y*nx*nz+x*nz+z]);
        }else{
          tmp+=0.5*delta_tau*(Ss1+Ss2);
        }
      }
      sbuf[y][x] = tmp;
    }
  }
  ctime=MPI_Wtime();
  MPI_Allreduce(
    ACCH::GetDevicePtr(sbuf[0]), ACCH::GetDevicePtr(rbuf[0]),
    nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
  atime+=MPI_Wtime()-ctime;
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_band[:ny*nx], rbuf[:ny][:nx])
  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
      I_band[y*nx+x]+=rbuf[y][x];
    }
  }
  }


  //ACCH::Free2D<double>(sbuf, ny, nx);
  //ACCH::Free2D<double>(rbuf, ny, nx);
  //
  ttime=MPI_Wtime()-ttime;
  if((myrank==0)&&(verbose>2)) printf("tau5000 time: %f %f %f %f \n",ttime,stime,atime,(stime+atime)/ttime);
  
}    
void RTS::tauscale_qrad(int band, double DX,double DY,double DZ, double * Ss){

  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  double idx=1.0/DX,idy=1.0/DY,idz=1.0/DZ;

  if(NDIM==1) {
    idx=0.;
    idy=0.;
  }
  if(NDIM==2) {
    idy=0.0;
  }

  //double ** sbuf = ACCH::Malloc2D<double>(ny, nx);
  //double ** rbuf = ACCH::Malloc2D<double>(ny, nx);
  //double **** Qtemp = ACCH::Malloc4D<double>(ny, nx, nz, 2);

#pragma acc enter data copyin(this[:1], Fx[:nx*ny*nz], Fy[:nx*ny*nz], Fz[:nx*ny*nz], \
        I_n[:nx*ny*nz], Qt[:(ny-yo)*(nx-xo)*(nz-zo)], Ss[:nx*ny*nz], \
        Tau[:nx*ny*nz],kap[:nx*ny*nz], rho[:nx*ny*nz],sbuf[:ny][:nx], rbuf[:ny][:nx])
#pragma acc parallel loop collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], kap[:nx*ny*nz], rho[:nx*ny*nz],sbuf[:ny][:nx], rbuf[:ny][:nx])  
  for(int y=0;y<ny;y++){ // loop over RT grid
    for(int x=0;x<nx;x++){
      Tau[y*nx*nz+x*nz+nz-1]=1.0e-12;
      for(int z=nz-2;z>=0;--z){
	double k0=kap[y*nx*nz+x*nz+z],r0=rho[y*nx*nz+x*nz+z],k_upw=kap[y*nx*nz+x*nz+z+1],r_upw=rho[y*nx*nz+x*nz+z+1];
	Tau[y*nx*nz+x*nz+z]=Tau[y*nx*nz+x*nz+z+1] + (DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6));
      }
      
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y*nx*nz+x*nz];
    }
  }

  double ctime=MPI_Wtime(); 
  MPI_Scan(
    ACCH::GetDevicePtr(sbuf[0]), ACCH::GetDevicePtr(rbuf[0]),
    nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
  stime+=MPI_Wtime()-ctime;
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], rbuf[:ny][:nx])
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]-=Tau[y*nx*nz+x*nz];
#pragma acc loop vector
      for(int z=0;z<nz;++z){
        Tau[y*nx*nz+x*nz+z]+=rbuf[y][x];
      }
    }
  }
  if (need_I){
    //  Outgoing Intensity at top (Long Characteristics)
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Tau[:nx*ny*nz], Ss[:nx*ny*nz], \
         rbuf[:ny][:nx], sbuf[:ny][:nx])
    for(int y=0;y<ny;++y){ // loop over RT grid
      for(int x=0;x<nx;++x){
        rbuf[y][x]=0.0;
        sbuf[y][x]=0.0;
#pragma acc loop vector
        for(int z=1;z<nz;++z){
          double Ss1 = Ss[y*nx*nz+x*nz+z];
          double Ss2 = Ss[y*nx*nz+x*nz+z-1];
          double delta_tau=Tau[y*nx*nz+x*nz+z-1]-Tau[y*nx*nz+x*nz+z];
          if(delta_tau>dtau_min){
            double edt=exp(-delta_tau);
            double c1=(1.0-edt)/delta_tau;
            sbuf[y][x]+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y*nx*nz+x*nz+z]);
          }else{
            sbuf[y][x]+=0.5*delta_tau*(Ss1+Ss2);
          }  
        }          
      }
    }
    ctime=MPI_Wtime(); 
    MPI_Allreduce(
      ACCH::GetDevicePtr(sbuf[0]), ACCH::GetDevicePtr(rbuf[0]),
      nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
    ctime+=MPI_Wtime()-ctime;
#pragma acc parallel loop collapse(2) \
 present(this[:1], I_o[:nx*ny], rbuf[:ny][:nx])
    for(int y=0;y<ny;++y){
      for(int x=0;x<nx;++x){
        I_o[y*nx+x]+=rbuf[y][x];
      }
    }

  }
//  radiative energy imbalance
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], I_n[:nx*ny*nz], kap[:nx*ny*nz], rho[:nx*ny*nz], \
         J_band[:nx*ny*nz], St[:nx*ny*nz], Jt[:nx*ny*nz], Ss[:nx*ny*nz])
  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
#pragma acc loop vector
      for(int z=0;z<nz;++z){
        I_n[y*nx*nz+x*nz+z]=kap[y*nx*nz+x*nz+z]*rho[y*nx*nz+x*nz+z]*(J_band[y*nx*nz+x*nz+z]-Ss[y*nx*nz+x*nz+z]);
        St[y*nx*nz+x*nz+z] +=Ss[y*nx*nz+x*nz+z];
        Jt[y*nx*nz+x*nz+z] +=J_band[y*nx*nz+x*nz+z];
      }
    }
  }

// 
  double inv_tau_0=1.0e1;
#pragma acc parallel loop gang collapse(2) \
present(this[:1], Fx[:nx*ny*nz], Fy[:nx*ny*nz], Fz[:nx*ny*nz], \
        I_n[:nx*ny*nz], Qt[:(ny-yo)*(nx-xo)*(nz-zo)], \
        Tau[:nx*ny*nz])
  for(int y=0;y<ny-yo;y++){
    for(int x=0;x<nx-xo;x++){
#pragma acc loop vector
      for(int z=0;z<nz-zo;z++){
        double qf1=
          (
            (
              Fz[y*nx*nz + x*nz + z+zo]+
              Fz[y*nx*nz + (x+xo)*nz + (z+zo)]+
              Fz[(y+yo)*nx*nz + x*nz + (z+zo)]+
              Fz[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fz[y*nx*nz + x*nz + z]+
              Fz[y*nx*nz + (x+xo)*nz + z]+
              Fz[(y+yo)*nx*nz + x*nz + z]+
              Fz[(y+yo)*nx*nz + (x+xo)*nz + z]
            )
          )*idz

          +

          (
            (
              Fx[y*nx*nz + (x+xo)*nz + z]+
              Fx[y*nx*nz + (x+xo)*nz + (z+zo)]+
              Fx[(y+yo)*nx*nz + (x+xo)*nz + z]+
              Fx[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fx[y*nx*nz + x*nz + z]+
              Fx[y*nx*nz + x*nz + (z+zo)]+
              Fx[(y+yo)*nx*nz + x*nz + z]+
              Fx[(y+yo)*nx*nz + x*nz + (z+zo)]
            )
          )*idx

          +

          (
            (
              Fy[(y+yo)*nx*nz + x*nz + z]+
              Fy[(y+yo)*nx*nz + x*nz + (z+zo)]+
              Fy[(y+yo)*nx*nz + (x+xo)*nz + z]+
              Fy[(y+yo)*nx*nz + (x+xo)*nz + (z+zo)]
            )-(
              Fy[y*nx*nz + x*nz + z]+
              Fy[y*nx*nz + x*nz + (z+zo)]+
              Fy[y*nx*nz + (x+xo)*nz + z]+
              Fy[y*nx*nz + (x+xo)*nz + (z+zo)]
            )
          )*idy;
        qf1*=-0.25e0; 
        double qj1 = I_n[y*nx*nz+x*nz+z] +
                     I_n[y*nx*nz+x*nz+(z+zo)] +
                     I_n[y*nx*nz+(x+xo)*nz+z] +
                     I_n[y*nx*nz+(x+xo)*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+x*nz+z] +
                     I_n[(y+yo)*nx*nz+x*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+z] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+(z+zo)];
        qj1*=0.5*PI; 
        double tau_local=tau(z+zo,x+xo,y+yo);
        double weight=exp(-tau_local*inv_tau_0);
        //Qtemp[y][x][z+zo][0]=qf1;
        //Qtemp[y][x][z+zo][1]=qj1;
        int QtmpAdr = (y*nx*nz*2)+(x*nz*2)+((z+zo)*2);
        Qtemp[QtmpAdr]=qf1;
        Qtemp[QtmpAdr+1]=qj1;
        Qt[(((y)*(nx-xo)+(x))*(nz-zo)+(z-zo+1))]+=weight*qj1+(1.0-weight)*qf1;
        
      }
    }
  }
/*
   ACCH::UpdateCPU(Qt,(ny-yo)*(nx-xo)*(nz-zo)*sizeof(double));
   ACCH::UpdateCPU(Qtemp,ny*nx*nz*2*sizeof(double));
   ACCH::UpdateCPU(I_n,ny*nx*nz*sizeof(double));
   for(int y = 0; y < ny-yo; y++)
    for(int x = 0; x < nx-xo; x++)
      for(int z = 0; z < nz-zo; z++) {
        int QtmpAdr = (y*nx*nz*2)+(x*nz*2)+((z+zo)*2);
        //double Qtest=Qt[((y*(nx-xo)+x)*(nz-zo)+(z-zo+1))];
        //double Qtemptest0 = Qtemp[(((y*nx+x)*nz+z+zo)*2+0)];
        //double Qtemptest1 = Qtemp[(((y*nx+x)*nz+z+zo)*2+1)];
        double Qtemptest0 = Qtemp[QtmpAdr];
        double Qtemptest1 = Qtemp[QtmpAdr+1];
        //fprintf(stdout,"Qtest: %21.15E\n",Qtest);
        //fprintf(stdout,"Qtest: %21.5E %21.5E \n",Qtemptest0,Qtemptest1);
        double qj1 = I_n[y*nx*nz+x*nz+z] +
                     I_n[y*nx*nz+x*nz+(z+zo)] +
                     I_n[y*nx*nz+(x+xo)*nz+z] +
                     I_n[y*nx*nz+(x+xo)*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+x*nz+z] +
                     I_n[(y+yo)*nx*nz+x*nz+(z+zo)] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+z] +
                     I_n[(y+yo)*nx*nz+(x+xo)*nz+(z+zo)];
        //fprintf(stdout,"qji: %21.15E\n",qj1);
      } 
*/
   if (save_col){
     ACCH::UpdateGPU3D<double>(Col_out, Nbands, col_nz, col_nvar);
     int col_bnd2 = col_bnd[2];
     int col_bnd3 = col_bnd[3];
     int col_bnd0 = col_bnd[0];
     int col_bnd1 = col_bnd[1];
#pragma acc parallel loop gang collapse(2) \
 present(this[:1], Col_out[:Nbands][:col_nz][:col_nvar], \
         J_band[:nx*ny*nz], Ss[:nx*ny*nz], kap[:nx*ny*nz], abn[:nx*ny*nz], \
         sig[:nx*ny*nz], B[:nx*ny*nz], Tau[:nx*ny*nz], Qtemp[:ny*nx*nz*2])
     for (int y=yl;y<=yh;++y){
       for (int x=xl;x<=xh;++x){
#pragma acc loop vector
         for (int z=zl+zo;z<=zh;++z){
           int ind = (y-yl)*nx*nz + (x-xl)*nz + (z-zl);
           Col_out[band][z-2*zo][0] += J_band[ind]*avg_col;
           Col_out[band][z-2*zo][1] += Ss[ind]*avg_col;
           Col_out[band][z-2*zo][2] += kap[ind]*avg_col;
           Col_out[band][z-2*zo][3] += abn[ind]*avg_col;
           Col_out[band][z-2*zo][4] += sig[ind]*avg_col;
           Col_out[band][z-2*zo][5] += B[ind]*avg_col;
           Col_out[band][z-2*zo][6] += Tau[ind]*avg_col;
           //Col_out[band][z-2*zo][7] += Qtemp[y-yl][x-xl][z-zl][0]*avg_col;
           //Col_out[band][z-2*zo][8] += Qtemp[y-yl][x-xl][z-zl][1]*avg_col;
           Col_out[band][z-2*zo][7] += Qtemp[((((y-yl)*nx+(x-xl))*nz+(z-zl))*2+0)]*avg_col;
           Col_out[band][z-2*zo][8] += Qtemp[((((y-yl)*nx+(x-xl))*nz+(z-zl))*2+1)]*avg_col;
         }
       }
     }
     ACCH::UpdateCPU3D<double>(Col_out, Nbands, col_nz, col_nvar);
   }

 // ACCH::Free4D<double>(Qtemp, ny, nx, nz, 2);
 // ACCH::Free2D<double>(sbuf, ny, nx);
 // ACCH::Free2D<double>(rbuf, ny, nx);
  //
  ttime=MPI_Wtime()-ttime;

  if((myrank==0)&&(verbose>2)) printf("tauscale time: %f %f %f %f \n",ttime,stime,atime,(stime+atime)/ttime);
}

void RTS::save_1D_avg(char * path, int iter,double time){

  // Save 1D avg
  
  static int ini_flag = 1;

  static MPI_Datatype x_subarray;

  if (ini_flag) {
    int array_of_sizes[1];
    int array_of_subsizes[1];
    int array_of_starts[1];

    array_of_sizes[0]=col_nzt;
    array_of_subsizes[0]=col_nz;
    array_of_starts[0]=col_offz;

    MPI_Type_create_subarray(1,array_of_sizes,array_of_subsizes,
                 array_of_starts,MPI_ORDER_FORTRAN,
                 MPI_FLOAT,&x_subarray);
    MPI_Type_commit(&x_subarray);

    ini_flag = 0;
  }

  double *** Col_out_glo = d3dim(0,Nbands-1,0,col_nz-1,0,col_nvar-1);
  memset(Col_out_glo[0][0],0,col_nz*Nbands*col_nvar*sizeof(double));

  MPI_Reduce(Col_out[0][0],Col_out_glo[0][0],col_nvar*Nbands*col_nz,MPI_DOUBLE,MPI_SUM,0,YZ_COMM);

  if(yz_rank==0){ // MPI_Reduce results are only meaningful on rank 0!

    int bufsize = Nbands*col_nz*col_nvar;
    float * iobuf = new float[bufsize];
 
    for(int band=0;band<Nbands;band++)
      for(int ind=0;ind<col_nz;ind++) 
        for(int var=0;var<col_nvar;var++)
          iobuf[band*(col_nz*col_nvar)+var*col_nz+ind] = (float) Col_out_glo[band][ind][var];

    char filename[128];
    sprintf(filename,"%s%s.%06d",path,"RT_mean1D",iter);

    MPI_File fhandle_mpi;

    MPI_File_open(XCOL_COMM,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
          MPI_INFO_NULL,&fhandle_mpi);

    if( xcol_rank == 0 ){
      float header[8];            
      header[0] = (float) Nbands;
      header[1] = (float) col_nvar;
      header[2] = (float) col_nzt;      
      header[3] = (float) time;
      header[4] = (float) col_bnd[0];
      header[5] = (float) col_bnd[1];
      header[6] = (float) col_bnd[2];      
      header[7] = (float) col_bnd[3];
      MPI_File_write(fhandle_mpi,&header[0],8,MPI_FLOAT,MPI_STATUS_IGNORE);
    }
    
    int offset = 8*sizeof(float);
      
    MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,x_subarray,(char*) "native", MPI_INFO_NULL);
    MPI_File_write_all(fhandle_mpi,&iobuf[0],bufsize,MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fhandle_mpi);  
    delete [] iobuf;
  }
  
  del_d3dim(Col_out_glo,0,Nbands-1,0,col_nz-1,0,col_nvar-1);

}

bool CheckDependency(const int step[4])
{
  return step[0] == step[1] &&
         step[0] == step[2] &&
         step[0] == step[3];
}

void Transpose_Rho_Kap_B(
  double * rho, double * trho,
  double * kap, double * tkap,
  double * B,   double * tB,
  const int nx, const int ny, const int nz
)
{
  const int n = nx*ny*nz;

#pragma acc parallel loop collapse(3) \
 present(rho[:n], kap[:n], B[:n], \
         trho[:n], tkap[:n], tB[:n])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      for(int z = 0; z < nz; z++) {
        int ind = y*nx*nz + x*nz + z;
        int tind = z*ny*nx + y*nx + x;
        trho[tind] = rho[ind];
        tkap[tind] = kap[ind];
        tB[tind]   = B[ind];
      }
    }
  }
}

void Transpose_In(
  double * I_n, double * tI_n,
  const int nx, const int ny, const int nz
)
{
  const int n = nx*ny*nz;

#pragma acc parallel loop collapse(3) \
 present(I_n[:nx*ny*nz], tI_n[:nx*ny*nz])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      for(int z = 0; z < nz; z++) {
        int ind = y*nx*nz + x*nz + z;
        int tind = z*ny*nx + y*nx + x;
        tI_n[tind] = I_n[ind];
      }
    }
  }
}

void Untranspose_In(
  double * I_n, double * tI_n,
  const int nx, const int ny, const int nz
)
{
  const int n = nx*ny*nz;

#pragma acc parallel loop collapse(3) \
 present(I_n[:nx*ny*nz], tI_n[:nx*ny*nz])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      for(int z = 0; z < nz; z++) {
        int ind = y*nx*nz + x*nz + z;
        int tind = z*ny*nx + y*nx + x;
        I_n[ind] = tI_n[tind];
      }
    }
  }
}

void Transpose_integrate_kernel(
  double * I_n, double * coeff,
  const double c[4], const int off[4],
  const int bounds[3], const int str[3],
  const int strc[3], const int i[3],
  const int step, const int n
)
{
  GPU_COMPARE(I_n, "double", n, n*sizeof(double), "I_n", "rt.cc", "Transpose_integrate", 1)
  GPU_COMPARE(coeff, "double", n*2, n*2*sizeof(double), "coeff", "rt.cc", "Transpose_integrate", 1)

  const int bound0 = bounds[0];
  const int bound1 = bounds[1];
  const int bound2 = bounds[2];

  const double c0 = c[0];
  const double c1 = c[1];
  const double c2 = c[2];
  const double c3 = c[3];
  const int off0 = off[0];
  const int off1 = off[1];
  const int off2 = off[2];
  const int off3 = off[3];
  const int i0 = i[0];
  const int i1 = i[1];
  const int i2 = i[2];
  const int str0 = str[0];
  const int str1 = str[1];
  const int str2 = str[2];
  const int strc0 = strc[0];
  const int strc1 = strc[1];
  const int strc2 = strc[2];
  #pragma acc data \
  present(I_n[:n], coeff[:n*2])
  {
    #pragma acc loop seq
    for(int b0 = 0; b0 < bound0; b0++) {
      int ind0 = (i0+b0*step)*str0;
      int inu0 = b0*strc0;
      #pragma acc parallel loop gang async\
       default(present)
      for(int b1 = 0; b1 < bound1; b1++) {
        int ind1 = ind0 + (i1+b1)*str1;
        int inu1 = inu0 + b1*strc1;
        #pragma acc loop vector
        for(int b2 = 0; b2 < bound2; b2++) {
          int ind = ind1 + (i2+b2)*str2;
          int inu = inu1 + b2*strc2;
          double I_upw = c0*I_n[ind-off0] +
                         c1*I_n[ind-off1] +
                         c2*I_n[ind-off2] +
                         c3*I_n[ind-off3];
          I_n[ind] = I_upw*coeff[inu]+coeff[inu+n];
        }
      }
    }
    #pragma acc wait
  }
  GPU_COMPARE(I_n, "double", n, n*sizeof(double), "I_n", "rt.cc", "Transpose_integrate", 2)
}

void Transpose_integrate(
  double * I_n, double * tI_n, double * coeff,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR,
  const int ixstep[4], const int iystep[4], const int izstep[4],
  int x_i, int y_i, int z_i,
  const double c[4]
)
{
  int xstep = (x_i < nx/2 ? 1 : -1);
  int ystep = (y_i < ny/2 ? 1 : -1);
  int zstep = (z_i < nz/2 ? 1 : -1);

  int off[4];
  int b_i[3], stride[3], stride2[3], bound[3], step;

  if(CheckDependency(iystep)) { // y pattern
    if(xstep == -1) x_i = 0;
    if(zstep == -1) z_i = 0;
    b_i[0] = y_i; b_i[1] = x_i; b_i[2] = z_i;
    stride[0] = nx*nz; stride[1] = nz; stride[2] = 1;
    stride2[0] = (nx-1)*(nz-1); stride2[1] =  nz-1; stride2[2] = 1;
    bound[0] = ny-1; bound[1] = nx-1; bound[2] = nz-1;
    step = ystep;
    for(int i = 0; i < 4; i++) {
      off[i] = iystep[i]*nx*nz + ixstep[i]*nz + izstep[i];
    }
    Transpose_integrate_kernel(
      I_n, coeff, c, off, bound, stride, stride2, b_i, step, nx*ny*nz
    );
  } else if(CheckDependency(ixstep)) { // x pattern
    if(ystep == -1) y_i = 0;
    if(zstep == -1) z_i = 0;
    b_i[0] = x_i; b_i[1] = y_i; b_i[2] = z_i;
    stride[0] = nz; stride[1] = nx*nz; stride[2] = 1;
    stride2[0] = (ny-1)*(nz-1); stride2[1] = nz-1; stride2[2] = 1;
    bound[0] = nx-1; bound[1] = ny-1; bound[2] = nz-1;
    step = xstep;
    for(int i = 0; i < 4; i++) {
      off[i] = iystep[i]*nx*nz + ixstep[i]*nz + izstep[i];
    }
    Transpose_integrate_kernel(
      I_n, coeff, c, off, bound, stride, stride2, b_i, step, nx*ny*nz
    );
  } else { // z pattern
    if(ystep == -1) y_i = 0;
    if(xstep == -1) x_i = 0;
    b_i[0] = z_i; b_i[1] = y_i; b_i[2] = x_i;
    stride[0] = ny*nx; stride[1] = nx; stride[2] = 1;
    stride2[0] = (ny-1)*(nx-1); stride2[1] = nx-1; stride2[2] = 1;
    bound[0] = nz-1; bound[1] = ny-1; bound[2] = nx-1;
    step = zstep;
    for(int i = 0; i < 4; i++) {
      off[i] = izstep[i]*ny*nx + iystep[i]*nx + ixstep[i];
    }
    Transpose_integrate_kernel(
      tI_n, coeff, c, off, bound, stride, stride2, b_i, step, nx*ny*nz
    );
  }
}
void Transpose_interpol_kernel(
  double * rho, double * kap, double * Ss, double * coeff,
  const double c[4], const int off[4],
  const int bounds[3], const int str[3],
  const int strc[3], const int i[3],
  const int step, const int n,
  const double ds3, const double ds6
)
{
  GPU_COMPARE(rho, "double", n, n*sizeof(double), "rho", "rt.cc", "Transpose_interpol", 1)
  GPU_COMPARE(kap, "double", n, n*sizeof(double), "kap", "rt.cc", "Transpose_interpol", 1)
  GPU_COMPARE(Ss, "double", n, n*sizeof(double), "Ss", "rt.cc", "Transpose_interpol", 1)
  GPU_COMPARE(coeff, "double", n*2, n*2*sizeof(double), "coeff", "rt.cc", "Transpose_interpol", 1)

  const int bound0 = bounds[0];
  const int bound1 = bounds[1];
  const int bound2 = bounds[2];

  const int i0 = i[0];
  const int i1 = i[1];
  const int i2 = i[2];

  const int str0 = str[0];
  const int str1 = str[1];
  const int str2 = str[2];

  const int strc0 = strc[0];
  const int strc1 = strc[1];
  const int strc2 = strc[2];

  const int off0 = off[0];
  const int off1 = off[1];
  const int off2 = off[2];
  const int off3 = off[3];

  const double c0 = c[0];
  const double c1 = c[1];
  const double c2 = c[2];
  const double c3 = c[3];

#pragma acc parallel loop gang collapse(2) \
 present(rho[:n], kap[:n], Ss[:n], coeff[:n*2]) \
 copyin(i[:3], str[:3], strc[:3], off[:4])
  for(int b0 = 0; b0 < bound0; b0++) {
    for(int b1 = 0; b1 < bound1; b1++) {
      int ind1 = (i[0]+b0*step)*str[0] + (b1+i[1])*str[1];
      int inu1 = b0*strc[0] + b1*strc[1];
#pragma acc loop vector
      for(int b2 = 0; b2 < bound2; b2++) {
        int ind = ind1 + (b2+i[2])*str[2];
        int inu = inu1 + b2*strc[2];

        double r_upw = c0*rho[ind-off[0]] +
                       c1*rho[ind-off[1]] +
                       c2*rho[ind-off[2]] +
                       c3*rho[ind-off[3]];
        double k_upw = c0*kap[ind-off[0]] +
                       c1*kap[ind-off[1]] +
                       c2*kap[ind-off[2]] +
                       c3*kap[ind-off[3]];
        double S_upw = c0*Ss[ind-off[0]] +
                       c1*Ss[ind-off[1]] +
                       c2*Ss[ind-off[2]] +
                       c3*Ss[ind-off[3]];
        double r0 = rho[ind];
        double k0 = kap[ind];
        double S0 = Ss[ind];

        double dt = ds3*(k_upw*r_upw+k0*r0) + ds6*(k0*r_upw+k_upw*r0);
        double expo = exp(-dt);
        double w0, w1;
        if (dt > dtau_min) {
          w0=1.0-expo;
          w1=w0-dt*expo;
        } else {
          w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
          w1=dt*dt/2.0-dt*dt*dt/3.0;
        }
        double source=S0*(w0-w1/dt)+S_upw*(w1/dt);

        if (dt > dtau_min2) {
          coeff[inu] = expo;
          coeff[inu+n] = source;
        } else {
          coeff[inu] = 1.0;
          coeff[inu+n] = 0.0;
        }

      }
    }
  }

  GPU_COMPARE(coeff, "double", n*2, n*2*sizeof(double), "coeff", "rt.cc", "Transpose_interpol", 2)
}

void Transpose_interpol(
  double * rho, double * trho,
  double * kap, double * tkap,
  double * Ss, double * tSs, double * coeff,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l,
  const int ixstep[4], const int iystep[4], const int izstep[4],
  int x_i, int y_i, int z_i, const double c[4], const double ds_upw[3]
)
{
  double ds3 = ds_upw[l]*inv3, ds6 = ds_upw[l]*inv6;

  int xstep = (x_i < nx/2 ? 1 : -1);
  int ystep = (y_i < ny/2 ? 1 : -1);
  int zstep = (z_i < nz/2 ? 1 : -1);

  int off[4];
  int b_i[3], stride[3], stride2[3], bound[3], step;

#ifdef PGICOMPARE
  #pragma acc parallel loop present(coeff[:nx*ny*nz*2])
  for(int i = 0; i < nx*ny*nz*2; i++)
    coeff[i] = 0;
#endif


  if(CheckDependency(iystep)) { // y pattern
    if(xstep == -1) x_i = 0;
    if(zstep == -1) z_i = 0;
    b_i[0] = y_i; b_i[1] = x_i; b_i[2] = z_i;
    stride[0] = nx*nz; stride[1] = nz; stride[2] = 1;
    stride2[0] = (nx-1)*(nz-1); stride2[1] =  nz-1; stride2[2] = 1;
    bound[0] = ny-1; bound[1] = nx-1; bound[2] = nz-1;
    step = ystep;
    for(int i = 0; i < 4; i++) {
      off[i] = iystep[i]*nx*nz + ixstep[i]*nz + izstep[i];
    }
    Transpose_interpol_kernel(
      rho, kap, Ss, coeff,
      c, off,
      bound, stride,
      stride2, b_i,
      step, nx*ny*nz,
      ds3, ds6
    );
  } else if(CheckDependency(ixstep)) { // x pattern
    if(ystep == -1) y_i = 0;
    if(zstep == -1) z_i = 0;
    b_i[0] = x_i; b_i[1] = y_i; b_i[2] = z_i;
    stride[0] = nz; stride[1] = nx*nz; stride[2] = 1;
    stride2[0] = (ny-1)*(nz-1); stride2[1] = nz-1; stride2[2] = 1;
    bound[0] = nx-1; bound[1] = ny-1; bound[2] = nz-1;
    step = xstep;
    for(int i = 0; i < 4; i++) {
      off[i] = iystep[i]*nx*nz + ixstep[i]*nz + izstep[i];
    }
    Transpose_interpol_kernel(
      rho, kap, Ss, coeff,
      c, off,
      bound, stride,
      stride2, b_i,
      step, nx*ny*nz,
      ds3, ds6
    );
  } else { // z pattern
    if(ystep == -1) y_i = 0;
    if(xstep == -1) x_i = 0;
    b_i[0] = z_i; b_i[1] = y_i; b_i[2] = x_i;
    stride[0] = ny*nx; stride[1] = nx; stride[2] = 1;
    stride2[0] = (ny-1)*(nx-1); stride2[1] = nx-1; stride2[2] = 1;
    bound[0] = nz-1; bound[1] = ny-1; bound[2] = nx-1;
    step = zstep;
    for(int i = 0; i < 4; i++) {
      off[i] = izstep[i]*ny*nx + iystep[i]*nx + ixstep[i];
    }
    Transpose_interpol_kernel(
      trho, tkap, tSs, coeff,
      c, off,
      bound, stride,
      stride2, b_i,
      step, nx*ny*nz,
      ds3, ds6
    );
  }
}

void Transpose_readbuf_kernel(
  real * xsb, real * xrb, real * xob,
  real * ysb, real * yrb, real * yob,
  real * zsb, real * zrb, real * zob,
  double * I_n, const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR,
  const int strx, const int stry, const int strz
)
{
  int y0=(YDIR==FWD)?0:ny-1;


#pragma acc parallel loop collapse(2) \
 present(yrb[:nx*nz], I_n[:nx*ny*nz]) 
  for(int x=0;x<nx;++x) {
    for(int z=0;z<nz;++z) {
      I_n[y0*stry+x*strx+z*strz]= (double) yrb[z*nx+x];
    }
  }
#pragma acc parallel loop \
 present(yob[:nx*nz], ysb[:nx*nz])
  for(int i = 0; i < nx*nz; i++) {
    yob[i] = ysb[i];
  }
  int x0=(XDIR==RIGHT)?0:nx-1;
#pragma acc parallel loop collapse(2) \
 present(xrb[:ny*nz], I_n[:nx*ny*nz])
  for(int y=0;y<ny;++y) {
    for(int z=0;z<nz;++z) {
      I_n[y*stry+x0*strx+z*strz]= (double) xrb[z*ny+y];
    }
  }
#pragma acc parallel loop \
 present(xob[:ny*nz], xsb[:ny*nz])
  for(int i = 0; i < ny*nz; i++) {
    xob[i] = xsb[i];
  }

  int z0=(ZDIR==UP)?0:nz-1;
#pragma acc parallel loop collapse(2) \
 present(zrb[:ny*nx], I_n[:nx*ny*nz])
  for(int y=0;y<ny;++y) {
    for(int x=0;x<nx;++x) {
      I_n[y*stry+x*strx+z0*strz]= (double) zrb[x*ny+y];
    }
  }
#pragma acc parallel loop \
 present(zob[:ny*nx], zsb[:ny*nx])
  for(int i = 0; i < nx*ny; i++) {
    zob[i] = zsb[i];
  }
}
 
void Transpose_readbuf(
  real * xsb, real * xrb, real * xob,
  real * ysb, real * yrb, real * yob,
  real * zsb, real * zrb, real * zob,
  double * I_n, double * tI_n, const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l, const int band,
  const int ixstep[4], const int iystep[4], const int izstep[4]
)
{
  int adrX = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nz);
  int adrY = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*nx*nz);
  int adrZ = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nx);

   
  if(CheckDependency(izstep)) {
    Transpose_readbuf_kernel(
      &xsb[adrX],
      &xrb[adrX],
      &xob[adrX],
      &ysb[adrY],
      &yrb[adrY],
      &yob[adrY],
      &zsb[adrZ],
      &zrb[adrZ],
      &zob[adrZ],
      tI_n, nx, ny, nz, XDIR, YDIR, ZDIR,
      1, nx, ny*nx
    );
  } else {
    Transpose_readbuf_kernel(
      &xsb[adrX],
      &xrb[adrX],
      &xob[adrX],
      &ysb[adrY],
      &yrb[adrY],
      &yob[adrY],
      &zsb[adrZ],
      &zrb[adrZ],
      &zob[adrZ],
      I_n, nx, ny, nz, XDIR, YDIR, ZDIR,
      nz, nx*nz, 1
    );
  }
}

void Transpose_writebuf_kernel(
  real * xsb, real * ysb, real * zsb, double * I_n,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR,
  const int strx, const int stry, const int strz
)
{
  int y0 = (YDIR==FWD) ? ny - 1 : 0;
#pragma acc parallel loop collapse(2) \
 present(ysb[:nx*nz], I_n[:nx*ny*nz])
  for(int x=0;x<nx;++x) {
    for(int z=0;z<nz;++z) {
      ysb[z*nx+x]=(real) I_n[y0*stry+x*strx+z*strz];
    }
  }
  int x0 = (XDIR==RIGHT) ? nx - 1 : 0;
#pragma acc parallel loop collapse(2) \
 present(xsb[:ny*nz], I_n[:nx*ny*nz])
  for(int y=0;y<ny;++y) {
    for(int z=0;z<nz;++z) {
      xsb[z*ny+y]=(real) I_n[y*stry+x0*strx+z*strz];
    }
  }
  int z0 = (ZDIR==UP) ? nz - 1 : 0;
#pragma acc parallel loop collapse(2) \
 present(zsb[:ny*nx], I_n[:nx*ny*nz])
  for(int y=0;y<ny;++y) {
    for(int x=0;x<nx;++x) {
      zsb[x*ny+y]=(real) I_n[y*stry+x*strx+z0*strz];
    }
  }
}

void Transpose_writebuf(
  real * xsb, real * ysb, real * zsb, double * I_n, double * tI_n,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l, const int band,
  int ixstep[4], int iystep[4], int izstep[4]
)
{
  int adrX = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nz);
  int adrY = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*nx*nz);
  int adrZ = (((((band*2+YDIR)*2+XDIR)*2+ZDIR)*NMU+l)*ny*nx);
  
  if(CheckDependency(izstep)) {
    Transpose_writebuf_kernel(
      &xsb[adrX],
      &ysb[adrY],
      &zsb[adrZ],
      tI_n, nx, ny, nz, XDIR, YDIR, ZDIR,
      1, nx, ny*nx
    );
  } else {
    Transpose_writebuf_kernel(
      &xsb[adrX],
      &ysb[adrY],
      &zsb[adrZ],
      I_n, nx, ny, nz, XDIR, YDIR, ZDIR,
      nz, nx*nz, 1
    );
  }
}

void Transpose_flux_kernel(
  double * I_n, double * Fx, double * Fy, double * Fz, double * J_band,
  const double c_J, const double c_x, const double c_y, const double c_z,
  const int nx, const int ny, const int nz
)
{
#pragma acc parallel loop collapse(2) gang \
 present(I_n[:nx*ny*nz], J_band[:nx*ny*nz], \
         Fz[:nx*ny*nz], Fy[:nx*ny*nz], Fx[:nx*ny*nz])
  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      int yxind = y*nx*nz + x*nz;
#pragma acc loop vector
      for(int z = 0; z < nz; z++) {
        int ind = yxind + z;
        double tmp = I_n[ind];
        J_band[ind] += c_J*tmp;
        Fz[ind]     += c_z*tmp;
        Fx[ind]     += c_x*tmp;
        Fy[ind]     += c_y*tmp;
      }
    }
  }

}

void Transpose_flux(
  double * I_n, double * tI_n, double * Fx, double * Fy, double * Fz, double * J_band,
  double wmu[3], double xmu[3][3],
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l,
  int ixstep[4], int iystep[4], int izstep[4]
)
{
  if(CheckDependency(izstep)) {
    Untranspose_In(I_n, tI_n, nx, ny, nz);
  }
  double zsign=(ZDIR==UP)?1.0:-1.0,xsign=(XDIR==RIGHT)?1.0:-1.0,ysign=(YDIR==FWD)?1.0:-1.0;
  double c_J = 0.125*wmu[l];
  double c_z = 0.5*PI*zsign*wmu[l]*xmu[2][l];
  double c_x = 0.5*PI*xsign*wmu[l]*xmu[0][l];
  double c_y = 0.5*PI*ysign*wmu[l]*xmu[1][l];

  Transpose_flux_kernel(
    I_n, Fx, Fy, Fz, J_band, c_J, c_x, c_y, c_z, nx, ny, nz
  );
}

