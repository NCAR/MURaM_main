#include <mpi.h>
#include <fstream>
#include <string.h>
#include <errno.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "eos.H"
#include "comm_split.H"
#include "rt/rt.h"
#include "muramacc.H"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;



const double p1_dmp = 0.8;

#define YZ_GHOST_LOOP(G,i,j) \
  for((j)=(G).lbeg[2]-Grid.ghosts[2];(j)<=(G).lend[2]+Grid.ghosts[2];(j)++) \
  for((i)=(G).lbeg[1]-Grid.ghosts[1];(i)<=(G).lend[1]+Grid.ghosts[1];(i)++)

#define YZ_LOOP(G,i,j) \
  for((j)=(G).lbeg[2];(j)<=(G).lend[2];(j)++) \
  for((i)=(G).lbeg[1];(i)<=(G).lend[1];(i)++)

extern int myrank;

double p_bc, s_bc;
void set_p_s_bc(const GridData&, const PhysicsData&, const RunData&);
#ifdef MURAM_FFTW
extern void potential_ext_fftw3(const GridData&,double*, double*);
extern void potential_ext_fftw3_par(const GridData&,double*, double*);
#else
extern void potential_ext_heffte(const GridData&,double*, double*);
#endif

void SetInitialConditions(GridData& Grid, const PhysicsData& Physics){}

/************************************************************************/
void WriteBackupFile(const char* file,const int iter,const double time) {
  FILE* fh;
  fh=fopen(file,"w");
  fprintf(fh,"%d  %.16lg  %lg  %lg",
	  iter,time,p_bc,s_bc);
  fclose(fh);
}

void ReadBackupFile(const char *file,const int bc,int* iter,double* time) {
  std::ifstream fptr(file,std::ios::in);
  if(fptr) {
    fptr.precision(16);
    if(bc)
      fptr >> *iter >> *time >> p_bc >> s_bc;
    else
      fptr >> *iter >> *time;
    fptr.close();
  }
}

void WriteBClog(const char *file,int iter,double time,RTS *rts){
  static int ini = 1;
  static FILE* fh;
  if(ini){
    fh  = fopen(file,"a");
    ini = 0;
  }
  if(fh){
    fprintf(fh,"%d  %.16lg  %lg  %lg  ||  %lg\n",
	    iter,time,p_bc,s_bc,rts->Fout());
    fflush(fh);
  } else {
    cout << file << " is not open" << endl;
  }
}
/************************************************************************/

void SetBoundaryConditions(const RunData&  Run, GridData& Grid,
                           const PhysicsData& Physics, 
			   const int stage, const int pt_update,RTS *rts){

  static int bnd_ini_flag = 1;

  double ttime;
  ttime=MPI_Wtime();
  
  const int POTENTIAL_BC = (int) Physics.bnd[i_bnd_pot];
  const int OPEN_TOP     = (int) Physics.bnd[i_bnd_top];
  const double B_crit    = Physics.bnd[i_bnd_bcrit]/sqrt(8.0*asin(1.0));
  const double eps_top   = Physics.bnd[i_bnd_eps_top];
  
  const double b_unit = sqrt(8.0*asin(1.0));

  register int i,j,k;

  double c,eps1,eps2;
  double pmin1,pmin2,p1,p2,pmag,pmean1,pmean2;
  double rho1,rho2,s1,smin1,s2,smin2;
  double lbuf[9],gbuf[9];
  double Bh_rms,By_rms,Bx_m,Bz_m,B_eq,vx_rms;

  cState U0,U1,U2; 

  int node,iloc;

  const int localsize       = Grid.lsize[1]*Grid.lsize[2];
  const int localghostsize  = (Grid.lsize[1]+2*Grid.ghosts[1])*
                              (Grid.lsize[2]+2*Grid.ghosts[2]);

  static double *bz0;
  static double *b_ext;
  /**********************************************************************/

  if (bnd_ini_flag == 1){

    if( (p_bc <0.0) or (s_bc < 0.0) )
      set_p_s_bc(Grid,Physics,Run);

/* Allocate regardless so that GPUs do not segfault
    if (POTENTIAL_BC > 0) {
      if( Grid.procs[0] < 8){
	// do all FFTs serial on x_rank = top_rank & y_rank = 0 (MPI parallel in z) 
	if (Grid.is_gend[0]) {
          bz0   = new double[localsize];
	}
      } else {
	// do all FFTs in parallel on top 8 x_ranks & y_rank = 0 (MPI parallel in z)
	// need bz0 on all x_ranks
	bz0   = new double[localsize];
      }
      if (Grid.is_gend[0]) {
        b_ext = new double[8*localghostsize];
      }
    }
*/

    bz0 = new double[localsize];
    b_ext = new double[8*localghostsize];
    
    if(Run.rank == 0){
      cout << "New pressure bnd: mean pressure fixed, fluctuations damped" 
	   << endl;
      cout << "bnd: POTENTIAL_BC = "  << POTENTIAL_BC << endl;
      cout << "bnd: OPEN_TOP     = "  << OPEN_TOP     << endl;
      cout << "bnd: p1_dmp       = "  << p1_dmp       << endl;
      cout << "bnd: p_bc         = "  << p_bc         << endl;
      cout << "bnd: s_bc         = "  << s_bc         << endl;  
    }
    bnd_ini_flag = 0;
  }

  if( (Run.rank == 0) && (stage == 1) )
    WriteBClog("BC.log",Run.globiter,Run.time,rts);

  /**********************************************************************/
 
  if( Grid.is_gbeg[0] ){
    //for(i=0;i<9;i++) lbuf[i] = 0.0;

    double lbuf0, lbuf1, lbuf2, lbuf3, lbuf4, lbuf5, lbuf6, lbuf7, lbuf8;
    lbuf0=0, lbuf1=0, lbuf2=0, lbuf3=0, lbuf4=0, lbuf5=0;
    lbuf6=0, lbuf7=0, lbuf8=0;

    const int kbeg = Grid.lbeg[2], kend = Grid.lend[2];
    const int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
    const int lbeg0 = Grid.beg[0];
    const int bufsize = Grid.bufsize;

#pragma acc parallel loop collapse(2) \
 present(Grid[:1], Grid.U[:bufsize]) \
 private(U1, U2, eps1, eps2, pmin1, pmin2) \
 reduction(+:lbuf0) reduction(+:lbuf1) reduction(+:lbuf2) \
 reduction(+:lbuf3) reduction(+:lbuf4) reduction(+:lbuf5) \
 reduction(+:lbuf6) reduction(+:lbuf7) reduction(+:lbuf8) 
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
      U1 = Grid.U[Grid.node(Grid.lbeg[0],j,k)];
      U2 = Grid.U[Grid.node(Grid.lbeg[0]+1,j,k)];
      
      U1.e -= 0.5*(U1.M.sqr()/U1.d);
      U2.e -= 0.5*(U2.M.sqr()/U2.d);	

      eps1 = U1.e/U1.d;
      eps2 = U2.e/U2.d;
      pmin1 = p_interp(eps1,U1.d); 
      pmin2 = p_interp(eps2,U2.d);

      if(U1.M.x > 0.0){
        lbuf0 += 1.0;
        lbuf1 += U1.B.y*U1.B.y+U1.B.z*U1.B.z;
        lbuf2 += U1.B.x*U1.B.x;
        lbuf3 += U1.B.y;
        lbuf4 += U1.B.z;
      }
      lbuf5 += pmin1;
      lbuf6 += pmin2;
      lbuf7 += U1.M.sqr()/U1.d;
      lbuf8 += pow(U1.M.x/U1.d,2);
    }

    lbuf[0] = lbuf0; 
    lbuf[1] = lbuf1; 
    lbuf[2] = lbuf2; 
    lbuf[3] = lbuf3; 
    lbuf[4] = lbuf4; 
    lbuf[5] = lbuf5; 
    lbuf[6] = lbuf6; 
    lbuf[7] = lbuf7; 
    lbuf[8] = lbuf8; 

    PGI_COMPARE(lbuf, double, 9, "lbuf", "boundary_pdmp_1_fftw3.C", "SetBoundaryConditions", 1)

    MPI_Allreduce(lbuf,gbuf,9,MPI_DOUBLE,MPI_SUM,YZ_COMM);
    Bh_rms = sqrt(gbuf[1]/gbuf[0]);
    By_rms = sqrt(gbuf[2]/gbuf[0]);
    Bx_m   = gbuf[3]/gbuf[0];
    Bz_m   = gbuf[4]/gbuf[0];
    
    c=1.0/(double(Grid.gsize[1]*Grid.gsize[2]));
    pmean1 = gbuf[5]*c;
    pmean2 = gbuf[6]*c;
    B_eq   = sqrt(gbuf[7]*c);
    vx_rms = sqrt(gbuf[8]*c);
    
    if( (yz_rank == 0) && (Run.verbose > 2) && (stage == 1)){
      cout << "Boundary_B : " << B_eq*b_unit << ' ' 
	   << sqrt(Bh_rms*Bh_rms+By_rms*By_rms)*b_unit << ' '
	   << Bh_rms*b_unit  << ' ' << By_rms*b_unit  << ' '
	   << Bx_m*b_unit    << ' ' << Bz_m*b_unit << endl;
    }

    const int kghosts = Grid.ghosts[2];
    const int jghosts = Grid.ghosts[1];
    const bool spitzer = Physics.params[i_param_spitzer] > 0.0;
    const bool ambipolar = Physics.params[i_param_ambipolar] > 0.0;

#pragma acc parallel loop collapse(2) \
 present(Grid[:1], Grid.U[:bufsize], Grid.v_amb[:bufsize], \
         Grid.sflx[:bufsize], xlp3[:N_lp3], xs3[:N_s3], d3_eostab[:N_lp3][:N_s3], eps3_eostab[:N_lp3][:N_s3]) \
 private(U1, U2, eps1, eps2, pmag, pmin1, pmin2, smin1, smin2, s1, s2, \
         c, p1, p2, rho1, rho2)
    for(k=kbeg-kghosts; k<=kend+kghosts; k++)
    for(j=jbeg-jghosts; j<=jend+jghosts; j++) {
      U1 = Grid.U[Grid.node(Grid.lbeg[0],j,k)];
      U2 = Grid.U[Grid.node(Grid.lbeg[0]+1,j,k)];
      
      U1.e -= 0.5*U1.M.sqr()/U1.d;
      U2.e -= 0.5*U2.M.sqr()/U2.d;	
      
      eps1 = U1.e/U1.d;
      eps2 = U2.e/U2.d;
      
      pmag  = 0.5 * U1.B.x*U1.B.x;
      
      pmin1 = p_interp(eps1,U1.d)-pmean1; 
      pmin2 = p_interp(eps2,U2.d)-pmean2;
      
      smin1 = s_interp(eps1,U1.d); 
      smin2 = s_interp(eps2,U2.d);

      if ( U1.M.x >= 0.e0 ){
        s1 = s_bc;
        s2 = s_bc;	  
      } else {     
        s1 = smin1; 
        s2 = smin2;	    
      }
      
      // exponential extrapolation of mean pressure to p_bc
      c   = p_bc/sqrt(pmean1*pmean2);
      p1  = pmean1*c;
      p2  = pmean1*c*c;
      // add damped pressure fluctuation
      p1 += pmin1*p1_dmp;
      p2 += pmin1*p1_dmp*p1_dmp;
      // make sure damping doesn't influence effects due to pmag
      p1 -= pmag*(1.0-p1_dmp);
      p2 -= pmag*(1.0-p1_dmp*p1_dmp);

      // ---- fix footpoint of spot ---------------------
      c = U1.B.x/B_crit;
      c = pow(c,4);
      c = (1.0-c)/(1.0+c);
      
      U1.M *= c;
      U2.M *= c;
      
      // ------------------------------------------------
      
      rho1 = d3_interp(p1,s1);   
      rho2 = d3_interp(p2,s2); 	  	  
      eps1 = eps3_interp(p1,s1); 
      eps2 = eps3_interp(p2,s2);
    
      U1.d = rho1; 
      U2.d = rho2;
     
      U1.e = eps1*rho1 + 0.5*U1.M.sqr()/U1.d;
      U2.e = eps2*rho2 + 0.5*U2.M.sqr()/U2.d;
      
      Grid.U[Grid.node(Grid.lbeg[0]-1,j,k)] = U1;	 
      Grid.U[Grid.node(Grid.lbeg[0]-2,j,k)] = U2;

      // heatflux boundary conditions
      if(spitzer){
	Grid.sflx[Grid.node(Grid.lbeg[0]-1,j,k)] = 0.0;
	Grid.sflx[Grid.node(Grid.lbeg[0]-2,j,k)] = 0.0;
      }

      // ambipolar diffusion boundary conditions
      if(ambipolar){
	Grid.v_amb[Grid.node(Grid.lbeg[0]-1,j,k)] = Vector(0.0,0.0,0.0);
	Grid.v_amb[Grid.node(Grid.lbeg[0]-2,j,k)] = Vector(0.0,0.0,0.0);
      }
    }

    PGI_COMPARE(Grid.U, double, Grid.bufsize*8, "U", "boundary_pdmp_1_fftw3.C",
                "SetBoundaryConditions", 2)
    if(Physics.params[i_param_spitzer] > 0.0) {
      PGI_COMPARE(Grid.sflx, double, Grid.bufsize, "sflx", "boundary_pdmp_1_fftw3.C",
                  "SetBoundaryConditions", 3)
    }
    if(Physics.params[i_param_ambipolar] > 0.0) {
      PGI_COMPARE(Grid.v_amb, double, Grid.bufsize*3, "v_amb", "boundary_pdmp_1_fftw3.C",
                  "SetBoundaryConditions", 4)
    }

  } // endif Grid.is_gbeg[0]
//#pragma acc update self(Grid.U[:Grid.bufsize], Grid.sflx[:Grid.bufsize], Grid.v_amb[:Grid.bufsize])
  //--------POTENTIAL FIELD--------------------------------------------
  //-------------------------------------------------------------------

  if( (POTENTIAL_BC == 1) or (POTENTIAL_BC > 1 and pt_update == 1) ) {
    
    if (Grid.is_gend[0]){
      const int lsize2 = Grid.lsize[2];
      const int lsize1 = Grid.lsize[1];
#pragma acc parallel loop collapse(2) \
 present(Grid[:1], Grid.U[:Grid.bufsize]) \
 private(node) \
 copyout(bz0[:localsize])
      for (k=0;k<lsize2;k++){
	for (j=0;j<lsize1;j++){
	  node = Grid.node(Grid.lend[0],j+Grid.lbeg[1],k+Grid.lbeg[2]);
	  bz0[j+k*lsize1] = Grid.U[node].B.x;
	}
      }
      PGI_COMPARE(bz0, double, localsize, "bz0", "boundary_pdmp_1_fftw3.C",
                  "SetBoundaryConditions", 5)
    }
      
#ifdef MURAM_FFTW
    if( Grid.procs[0] < 8)
    {
      /* do all FFTs serial on x_rank = top_rank & y_rank = 0 (MPI parallel in z) */
      if(Grid.is_gend[0]) 
      {
        potential_ext_fftw3(Grid, bz0, b_ext);
      }
      if( (Run.verbose > 1) && (Run.rank == 0)) cout << "potential update" << endl;
    }
    else
    { 
      potential_ext_fftw3_par(Grid, bz0, b_ext);

      if( (Run.verbose > 1) && (Run.rank == 0)) cout << "potential update parallel" << endl;
    }
#else
    if(Grid.is_gend[0]) 
    {
      potential_ext_heffte(Grid, bz0, b_ext);
    }
    if( (Run.verbose > 1) && (Run.rank == 0)) cout << "potential update" << endl;
#endif


    PGI_COMPARE(b_ext, double, 8*localghostsize, "b_ext", "boundary_pdmp_1_fftw3.C",
                "SetBoundaryConditions", 6)

  }

  if (Grid.is_gend[0]){ 
  //-------------------------------------------------------------------
  // b_ext contains now the local potential field extrapolation---
  //-------------------------------------------------------------------

    const int kbeg = Grid.lbeg[2], kend = Grid.lend[2];
    const int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
    const int bufsize = Grid.bufsize;
    const int kghosts = Grid.ghosts[2];
    const int jghosts = Grid.ghosts[1];
    const int lsize1 = Grid.lsize[1];
    const int lend0 = Grid.lend[0];
    const bool spitzer = Physics.params[i_param_spitzer] > 0.0;
    const bool bnd = Physics.bnd[i_bnd_eps_top] > 0.0;
    const bool ambipolar = Physics.params[i_param_ambipolar] > 0.0;

// BUG cuda segfault when using vector sqr() method, but works in other places in SetBoundary.. very weird
#pragma acc parallel loop collapse(2) \
 present(Grid[:1], Grid.U[:bufsize], Grid.v_amb[:bufsize], \
         Grid.sflx[:bufsize]) \
 copyin(b_ext[:8*localghostsize]) \
 private(iloc, U0, U1, U2)
    for(k=kbeg-kghosts; k<=kend+kghosts; k++)
    for(j=jbeg-jghosts; j<=jend+jghosts; j++) {
      iloc = j + k*(Grid.lsize[1]+2*Grid.ghosts[1]);

      U0 = Grid.U[Grid.node(Grid.lend[0],j,k)];
      U1 = Grid.U[Grid.node(Grid.lend[0],j,k)];
      U2 = Grid.U[Grid.node(Grid.lend[0]-1,j,k)];

      if (OPEN_TOP == 1) {
        if (U1.M.x < 0.0){
          U1.M.x *= -1; 
          U2.M.x *= -1; 
        }
      } else {
        U1.M.x *= -1; 
        U2.M.x *= -1;
      }

      if(eps_top > 0.0){
        //U1.e = eps_top*U1.d + 0.5*U1.M.sqr()/U1.d;
        //U2.e = eps_top*U2.d + 0.5*U2.M.sqr()/U2.d;
        U1.e = eps_top*U1.d + 0.5*(U1.M.x*U1.M.x+U1.M.y*U1.M.y+U1.M.z*U1.M.z)/U1.d;
        U2.e = eps_top*U2.d + 0.5*(U2.M.x*U2.M.x+U2.M.y*U2.M.y+U2.M.z*U2.M.z)/U2.d;
      }

      if( POTENTIAL_BC > 0 ){
        U0.B.y = (double) b_ext[iloc+2*localghostsize];
        U0.B.z = (double) b_ext[iloc+5*localghostsize];

        U1.B.x = (double) b_ext[iloc+0*localghostsize];
        U1.B.y = (double) b_ext[iloc+3*localghostsize];
        U1.B.z = (double) b_ext[iloc+6*localghostsize];

        U2.B.x = (double) b_ext[iloc+1*localghostsize];
        U2.B.y = (double) b_ext[iloc+4*localghostsize];
        U2.B.z = (double) b_ext[iloc+7*localghostsize];
      } else {
        U1.B.y *= -1;
        U2.B.y *= -1;
        U1.B.z *= -1;	
        U2.B.z *= -1;
      }
	
      Grid.U[Grid.node(Grid.lend[0],j,k)]   = U0;	  
      Grid.U[Grid.node(Grid.lend[0]+1,j,k)] = U1;	  
      Grid.U[Grid.node(Grid.lend[0]+2,j,k)] = U2;

      // heatflux boundary conditions
      if(spitzer){
        if(bnd ){
	  // hot plate at top, set sflx symmetric
	  Grid.sflx[Grid.node(Grid.lend[0]+1,j,k)] = Grid.sflx[Grid.node(Grid.lend[0],j,k)];
          Grid.sflx[Grid.node(Grid.lend[0]+2,j,k)] = Grid.sflx[Grid.node(Grid.lend[0]-1,j,k)];
        } else {
          // zero out top 4 layers to ensure zero flux at boundary
          Grid.sflx[Grid.node(Grid.lend[0]-1,j,k)]  = 0.0;
          Grid.sflx[Grid.node(Grid.lend[0],j,k)]    = 0.0;
          Grid.sflx[Grid.node(Grid.lend[0]+1,j,k)]  = 0.0;
          Grid.sflx[Grid.node(Grid.lend[0]+2,j,k)]  = 0.0;
        }   
      }  

      // ambipolar diffusion boundary conditions
      if(ambipolar){
	Grid.v_amb[Grid.node(Grid.lend[0]+1,j,k)] = Vector(0.0,0.0,0.0);
	Grid.v_amb[Grid.node(Grid.lend[0]+2,j,k)] = Vector(0.0,0.0,0.0);
      }
    }
   
    PGI_COMPARE(Grid.U, double, Grid.bufsize*8, "U", "boundary_pdmp_1_fftw3.C",
                "SetBoundaryConditions", 7)
    if(Physics.params[i_param_spitzer] > 0.0) {
      PGI_COMPARE(Grid.sflx, double, Grid.bufsize, "sflx", "boundary_pdmp_1_fftw3.C",
                  "SetBoundaryConditions", 8)
    }
    if(Physics.params[i_param_ambipolar] > 0.0) {
      PGI_COMPARE(Grid.v_amb, double, Grid.bufsize*3, "v_amb", "boundary_pdmp_1_fftw3.C",
                  "SetBoundaryConditions", 9)
    }
 
  }// endif Grid.is_gend[0]   
  
  ttime = MPI_Wtime() - ttime;
  
  if ( (Run.rank == 0) and (Run.verbose  > 2) && (stage == 1) ) cout << "boundary time = " << ttime << endl;
  
}
// ***************************************************************************

void set_p_s_bc(const GridData& Grid, const PhysicsData& Physics, const RunData& Run) {

  register int i,k;
  
  double eps1,eps2;
  double bnd[5],bnd_loc[5];
  double bnd0, bnd1, bnd2, bnd3, bnd4;
  cState U1,U2;
  
  if (Grid.is_gbeg[0] ){
//    for(i=0;i<5;i++)
//      bnd_loc[i] = 0.0;

    bnd0 = 0, bnd1 = 0, bnd2 = 0, bnd3 = 0, bnd4 = 0;

    const int kbeg = Grid.lbeg[2], kend = Grid.lend[2];
    const int ibeg = Grid.lbeg[1], iend = Grid.lend[1];
    const int lbeg0 = Grid.lbeg[0];
    const int bufsize = Grid.bufsize;

#pragma acc parallel loop collapse(2) \
 present(Grid[:1], Grid.U[:bufsize]) \
 private(U1, U2, eps1, eps2) \
 reduction(+:bnd0) reduction(+:bnd1) reduction(+:bnd2) \
 reduction(+:bnd3) reduction(+:bnd4)
    for(k=kbeg; k<=kend; k++)
    for(i=ibeg; i<=iend; i++) {
      U1=Grid.U[Grid.node(Grid.lbeg[0],i,k)];
      U2=Grid.U[Grid.node(Grid.lbeg[0]+1,i,k)];
      U1.e -= 0.5*U1.M.sqr()/U1.d;
      U2.e -= 0.5*U2.M.sqr()/U2.d;      
      eps1 = U1.e/U1.d;
      eps2 = U2.e/U2.d;
      bnd0 += p_interp(eps1,U1.d);
      bnd1 += p_interp(eps2,U2.d);
      if(U1.M.x > 0.0){
	bnd2 += 1.0;
	bnd3 += s_interp(eps1,U1.d);
      }
	bnd4 += s_interp(eps1,U1.d);
    }

    bnd_loc[0] = bnd0;
    bnd_loc[1] = bnd1;
    bnd_loc[2] = bnd2;
    bnd_loc[3] = bnd3;
    bnd_loc[4] = bnd4;

    PGI_COMPARE(bnd_loc, double, 5, "bnd_loc", "boundary_pdmp_1_fftw3.C",
                "set_p_s_bc", 10)

    MPI_Allreduce(bnd_loc,bnd,5,MPI_DOUBLE,MPI_SUM,YZ_COMM);
    bnd[0] /= double(Grid.gsize[1]*Grid.gsize[2]);
    bnd[1] /= double(Grid.gsize[1]*Grid.gsize[2]);

    if(p_bc < 0.0)
      p_bc = bnd[0]*sqrt(bnd[0]/bnd[1]);

    if(s_bc < 0.0){
      if(bnd[2] > 0.0)
	s_bc = bnd[3]/bnd[2];
      else
	s_bc = bnd[4]/double(Grid.gsize[1]*Grid.gsize[2]);
    }
  }

  MPI_Bcast(&p_bc,1,MPI_DOUBLE,0,XCOL_COMM);
  MPI_Bcast(&s_bc,1,MPI_DOUBLE,0,XCOL_COMM);
    
}
