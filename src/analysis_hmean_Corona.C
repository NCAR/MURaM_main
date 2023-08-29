#include <mpi.h>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "analysis.H"
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "rt/rt.h"
#include "comm_split.H"
#include "limit_va.H"

using namespace std;

void AnalyzeSolution_VP(const RunData& Run,const GridData& Grid,
		        const PhysicsData& Physics, RTS * rts) {

  static int ini_flag = 1;

  static MPI_Datatype x_subarray;

  if (ini_flag) {
    int array_of_sizes[1];
    int array_of_subsizes[1];
    int array_of_starts[1];

    array_of_sizes[0]=Grid.gsize[0];
    array_of_subsizes[0]=Grid.lsize[0];
    array_of_starts[0]=Grid.beg[0]-Grid.gbeg[0];

    MPI_Type_create_subarray(1,array_of_sizes,array_of_subsizes,
			     array_of_starts,MPI_ORDER_FORTRAN,
			     MPI_FLOAT,&x_subarray);
    MPI_Type_commit(&x_subarray);

    ini_flag = 0;
  }

  register int i,j,k,node,ind,ioff,v;

    //added as part of IO port - SS
  register int kbeg, kend, jbeg, jend, ibeg, iend;
  int stride0,stride1,stride2;

  kbeg = Grid.lbeg[2];
  kend = Grid.lend[2];
  jbeg = Grid.lbeg[1];
  jend = Grid.lend[1];
  ibeg = Grid.lbeg[0];
  iend = Grid.lend[0];
  
  stride0 = Grid.stride[0];
  stride1 = Grid.stride[1];
  stride2 = Grid.stride[2];
  const int bufsize = Grid.bufsize;
  
    double tvar1 = 0.0;
    double tvar2 = 0.0;
    double tvar3 = 0.0;
    double tvar4 = 1.0;
    double tvar5 = 0.0;
    double tvar6 = 0.0;
    double tvar7 = 0.0;
    double tvar8 = 0.0;
    double Qres  = 0.0;
    double Qvis  = 0.0;
    double Qamb  = 0.0;
    
    double Qthin = 0.0;
    double QChr  = 0.0;
    double QCa   = 0.0;
    double QMg   = 0.0;
    double QH    = 0.0;
  //end 

  int next[3];
  double w1[3],w2[3];
  
  for(v=0;v<Grid.NDIM;v++){
    next[v] = Grid.stride[v];
    w1[v]   = 8./(12.*Grid.dx[v]);
    w2[v]   =-1./(12.*Grid.dx[v]);
  }

  char filename[128];

  MPI_File fhandle_mpi;
  int offset;

  //commented as part of IO porting
  //static const cState* U = Grid.U;

  const double b_unit = sqrt(8.0*asin(1.0));

  double ihsz=1.0/((double (Grid.gsize[1]))*(double (Grid.gsize[2])));

  double dn,eps,vx,vy,vz,bx,by,bz,vsqr,idn,rfac;

  static int nvar=50;

  double  *loc, *glo;
  float *iobuf;
  
  int bufsz = nvar*Grid.lsize[0];
 
  loc   = new double [bufsz];
  glo   = new double [bufsz];
  iobuf = new float[bufsz];

#pragma acc enter data copyin(loc[:bufsz]) 

#pragma acc parallel loop   
  for(j=0;j<bufsz;j++){ 
    loc[j] = 0.0;
  }
  
  ioff = Grid.lsize[0];

//added as part of IO porting - SS
  double rt_ext_i_ext_cor_temp = Physics.rt_ext[i_ext_cor];
  double rt_i_rt_tr_tem_temp = Physics.rt[i_rt_tr_tem];
  double rt_i_rt_tr_pre_temp = Physics.rt[i_rt_tr_pre];
  double params_i_param_spitzer_temp = Physics.params[i_param_spitzer];
  double params_i_param_ambipolar_temp = Physics.params[i_param_ambipolar];
  bool diag_temp = Run.diagnostics;

#pragma acc parallel loop collapse(3) gang vector \
present(Grid[:1], Grid.U[:bufsize], Grid.Qthin[:bufsize],   \
       Grid.temp[:bufsize], Grid.pres[:bufsize],  \
       Grid.Tau[:bufsize], Grid.divB[:bufsize],  \
       Grid.sflx[:bufsize], Grid.tvar1[:bufsize],   \
       Grid.tvar2[:bufsize], Grid.tvar3[:bufsize],  \
       Grid.tvar4[:bufsize], Grid.tvar5[:bufsize],  \
       Grid.tvar6[:bufsize], Grid.tvar7[:bufsize],  \
       Grid.tvar8[:bufsize], Grid.Qres[:bufsize],   \
       Grid.Qvis[:bufsize], Grid.Qamb[:bufsize],     \
       Grid.Qtot[:bufsize], loc[:bufsz]) \
private(i,j,k,node, ind, dn, idn, \
     vx, vy, vz, bx, by, bz, vsqr, eps, \
     tvar1, tvar2, tvar3, tvar4, tvar5, \
     tvar6, tvar7, tvar8, Qres, Qvis, Qamb, \
     Qthin, QChr, QCa, QMg, QH, rfac) 
  for(k=kbeg;k<=kend;k++){
  for(j=jbeg;j<=jend;j++){
  for(i=ibeg;i<=iend;i++){

    node = i*stride0+j*stride1+k*stride2;

    ind = i-Grid.ghosts[0];

    dn   = Grid.U[node].d;
    idn  = 1.0/dn;

    vx   = Grid.U[node].M.x;
    vy   = Grid.U[node].M.y;   
    vz   = Grid.U[node].M.z;   
    bx   = Grid.U[node].B.x;
    by   = Grid.U[node].B.y;   
    bz   = Grid.U[node].B.z;

    vsqr = vx*vx+vy*vy+vz*vz;

    eps  = Grid.U[node].e*idn;

    // Output only stuff in corona
    if(rt_ext_i_ext_cor_temp == 1)
      rfac = Grid.Qthin[node] != 0.0;
    else if (rt_ext_i_ext_cor_temp == 2)
      rfac = (Grid.temp[node] > rt_i_rt_tr_tem_temp) && (Grid.pres[node] < rt_i_rt_tr_pre_temp);
    else
      rfac = 0.0;

    // mean state
#pragma acc atomic update
    loc[ind+0*ioff]  += dn;
#pragma acc atomic update
    loc[ind+1*ioff]  += eps;
#pragma acc atomic update
    loc[ind+2*ioff]  += Grid.pres[node];
#pragma acc atomic update
    loc[ind+3*ioff]  += Grid.temp[node];
#pragma acc atomic update
    loc[ind+4*ioff]  += Grid.Tau[node];

    // v, B
#pragma acc atomic update
    loc[ind+5*ioff]  += vx*vx;
#pragma acc atomic update
    loc[ind+6*ioff]  += vy*vy+vz*vz;
#pragma acc atomic update
    loc[ind+7*ioff]  += bx*bx;  
#pragma acc atomic update
    loc[ind+8*ioff]  += (by*by+bz*bz);
#pragma acc atomic update
    loc[ind+9*ioff]  += bx;
#pragma acc atomic update
    loc[ind+10*ioff] += fabs(bx);
#pragma acc atomic update
    loc[ind+11*ioff] += sqrt(by*by+bz*bz);
#pragma acc atomic update
    loc[ind+12*ioff] += Grid.divB[node]*Grid.divB[node];
    
    // vertical MHD + conductive fluxes
#pragma acc atomic update
    loc[ind+13*ioff] += vx*dn;
#pragma acc atomic update
    loc[ind+14*ioff] += fabs(vx*dn);   
#pragma acc atomic update
    loc[ind+15*ioff] += vx*dn*0.5*vsqr;
#pragma acc atomic update
    loc[ind+16*ioff] += vx*(eps*dn+Grid.pres[node]);
#pragma acc atomic update
    loc[ind+17*ioff] += eps+Grid.pres[node]/dn;
#pragma acc atomic update
    loc[ind+18*ioff] += vx*(by*by+bz*bz)-bx*(vy*by+vz*bz);
#pragma acc atomic update
    loc[ind+19*ioff] -= bx*(vy*by+vz*bz);

    if(params_i_param_spitzer_temp > 0.0){
#pragma acc atomic update
      loc[ind+20*ioff] += Grid.sflx[node]*bx;
    }

    tvar1 = 0.0;
    tvar2 = 0.0;
    tvar3 = 0.0;
    tvar4 = 1.0;
    tvar5 = 0.0;
    tvar6 = 0.0;
    tvar7 = 0.0;
    tvar8 = 0.0;
    Qres  = 0.0;
    Qvis  = 0.0;
    Qamb  = 0.0;
    
    if (diag_temp){
      tvar1 = Grid.tvar1[node];
      tvar2 = Grid.tvar2[node];
      tvar3 = Grid.tvar3[node];
      tvar4 = Grid.tvar4[node];
      tvar5 = Grid.tvar5[node];
      tvar6 = Grid.tvar6[node];
      tvar7 = Grid.tvar7[node];
      tvar8 = Grid.tvar8[node];
      Qres  = Grid.Qres[node];
      Qvis  = Grid.Qvis[node];
      if(params_i_param_ambipolar_temp > 0.0)
	Qamb  = Grid.Qamb[node];
    }

    Qthin = 0.0;
    QChr  = 0.0;
    QCa   = 0.0;
    QMg   = 0.0;
    QH    = 0.0;

    if(rt_ext_i_ext_cor_temp >= 1)
      Qthin = Grid.Qthin[node];

    // Radiation quantities
#pragma acc atomic update
    loc[ind+21*ioff] += Grid.Qtot[node]*tvar4; // tvar4 -> efac
#pragma acc atomic update
    loc[ind+22*ioff] += Qthin*tvar4;// tvar4 -> efac
#pragma acc atomic update
    loc[ind+23*ioff] += QChr*tvar4;// tvar4 -> efac
#pragma acc atomic update
    loc[ind+24*ioff] += QCa*tvar4;// tvar4 -> efac
#pragma acc atomic update
    loc[ind+25*ioff] += QMg*tvar4;// tvar4 -> efac
#pragma acc atomic update
    loc[ind+26*ioff] += QH*tvar4;// tvar4 -> efac

    // terms in plasma energy equation
#pragma acc atomic update
    loc[ind+27*ioff] += dn*(eps+0.5*vsqr);
#pragma acc atomic update
    loc[ind+28*ioff] += tvar1; // div(F_adv)
#pragma acc atomic update
    loc[ind+29*ioff] += tvar2; // div(F_cond)
#pragma acc atomic update
    loc[ind+30*ioff] += tvar3; // wlrt
#pragma acc atomic update
    loc[ind+31*ioff] += tvar5; // m*g + damping
#pragma acc atomic update
    loc[ind+32*ioff] += tvar6; // Boris
#pragma acc atomic update
    loc[ind+33*ioff] += tvar7; // Tcheck
#pragma acc atomic update
    loc[ind+34*ioff] += tvar8; // numerical diff
#pragma acc atomic update
    loc[ind+35*ioff] += Qres;
#pragma acc atomic update
    loc[ind+36*ioff] += Qvis;
#pragma acc atomic update
    loc[ind+37*ioff] += Qamb;
    
    // terms in plasma energy equation in Corona only
#pragma acc atomic update
    loc[ind+38*ioff] += rfac; // mask defining corona volume
#pragma acc atomic update
    loc[ind+39*ioff] += rfac*dn*(eps+0.5*vsqr);
#pragma acc atomic update
    loc[ind+40*ioff] += rfac*Grid.Qtot[node]*tvar4; // tvar4 -> efac
#pragma acc atomic update
    loc[ind+41*ioff] += Qthin*tvar4;
#pragma acc atomic update
    loc[ind+42*ioff] += rfac*tvar1; // div(F_adv)
#pragma acc atomic update
    loc[ind+43*ioff] += rfac*tvar3; // wlrt
#pragma acc atomic update
    loc[ind+44*ioff] += rfac*tvar5; // m*g + damping
#pragma acc atomic update
    loc[ind+45*ioff] += rfac*tvar6; // Boris
#pragma acc atomic update
    loc[ind+46*ioff] += rfac*tvar7; // Tcheck
#pragma acc atomic update
    loc[ind+47*ioff] += rfac*tvar8; // numerical diff
#pragma acc atomic update
    loc[ind+48*ioff] += rfac*Qres;
#pragma acc atomic update
    loc[ind+49*ioff] += rfac*Qvis;  
  }
}
}

#pragma acc exit data copyout(loc[:bufsz]) 

  MPI_Reduce(loc,glo,bufsz,MPI_DOUBLE,MPI_SUM,0,YZ_COMM);
  
  if(yz_rank==0){ // MPI_Reduce results are only meaningful on rank 0!
    
    for(ind=0;ind<bufsz;ind++){
      glo[ind] *=  ihsz;
    }

    // B
    for(ind=0;ind<Grid.lsize[0];ind++){
      glo[ind+ 7*ioff] *= b_unit*b_unit;
      glo[ind+ 8*ioff] *= b_unit*b_unit;
      glo[ind+ 9*ioff] *= b_unit;
      glo[ind+10*ioff] *= b_unit;
      glo[ind+11*ioff] *= b_unit;
      glo[ind+12*ioff] *= b_unit*b_unit;
    }
      
    for(ind=0;ind<bufsz;ind++){ 
      iobuf[ind] = (float) glo[ind];
    }

    sprintf(filename,"%s%s.%06d",Run.path_2D,"hmean1D",Run.globiter);
    MPI_File_open(XCOL_COMM,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
		  MPI_INFO_NULL,&fhandle_mpi);

    if( xcol_rank == 0 ){
      float header[4];            
      header[0] = (float) nvar;
      header[1] = (float) Grid.gsize[0];
      header[2] = (float) 1.0;      
      header[3] = (float) Run.time;
      MPI_File_write(fhandle_mpi,header,4,MPI_FLOAT,MPI_STATUS_IGNORE);
    }
    
    offset = 4*sizeof(float);
      
    MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,x_subarray,(char*) "native",
		      MPI_INFO_NULL);
    MPI_File_write_all(fhandle_mpi,iobuf,nvar*Grid.lsize[0],MPI_FLOAT,
		       MPI_STATUS_IGNORE);
      
    MPI_File_close(&fhandle_mpi);    
  }
  
  delete[] loc;
  delete[] glo;
  delete[] iobuf;
  
}
