#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "limit_va.H"
#include "src_int_tck.H"
#include "muramacc.H"
#include <iostream>

using std::cout;
using std::endl;
using std::min;
using std::max;

#define OUTER_LOOP(G,i,j,d1,d2) \
  for((j)=(G)[(d2)][0];(j)<=(G)[(d2)][1];(j)++) \
  for((i)=(G)[(d1)][0];(i)<=(G)[(d1)][1];(i)++)

// Hall coefficient is c/(4*pi*e*n_e) , our table is ne so this is the rest.
const double hall_const = 2.99792458e10/(16.0*atan(1)*1.60217733E-19);

//-------------------------------------------------------------
double MHD_Residual(const RunData&  Run, GridData& Grid, 
		    const PhysicsData& Physics) {

  //double time,s_time;
  //static double t_time = 0.0 ,c_time = 0.0 , r_time = 0.0;
  //static int call_count = 0;
  //s_time = MPI_Wtime();

  static int call_count = 1;

  const double CFL_hyp  = 1.8;
  const double avr_max  = 0.25;
  const double eps_hyp  = 0.02; // -> hyperdiffusion used on heatflux, blows up during flare without it, otherwise not needed 

  const double k0_spt   = Physics.params[i_param_spitzer];

  const double Rgas = 8.314e7;
  const double X_H  = 0.7;

  // 1/6 of free streaming limit 0.25*ne*K*T*sqrt(K*T/m_e)
  const double sat_flx = (1.0+X_H)*pow(Rgas,1.5)*sqrt(1836.)/8.0;

  // corresponds to 1.5 rho cs^1.5
  //const double sat_flx  = 4.98277e12;

  cState* U  = (cState*) Grid.U;
  cState* R  = (cState*) Grid.Res;

  register int i,j,k,d,node,d1,d2,d3,ivar,offset,str;
  int strd2,strd3;
  const double eta = Physics.params[i_param_eta];

  const bool ambipolar        = (Physics.params[i_param_ambipolar] > 0.0);
  const bool spitzer          = (Physics.params[i_param_spitzer] > 0.0);
  const bool needs_curlB      = (eta > 0.0);
  const bool need_diagnostics = Run.need_diagnostics;
  
  const int nvar   = Physics.NVAR;
  const int vsize  = Grid.vsize;
  const int v_nvar = Grid.v_nvar;
  const int bufsize = Grid.bufsize;

  double dn,lf,cs2,va2,x2,x4,s,cmax,dx,dy,dz,dxmin,dxmax,dt_res,dt_cfl,p1,p2,p3,vel,mag,ekin,pmag,mflx;
  double dt_amb,amb_diff,c_lim,cmax_hyp,amb_diff_max,amb_vel_max,nu_hyp,nu_hyp_max,rho_i,nu_in;
  double kap_spt,invBB,tt32,f_lim,F_sat,F_spt,hyp_diff;

  double w1,w2,ww0,ww1,ww2;
  double dB2,dB3;

  int stride[3];
  int sz,gh;
  int kbeg, kend, jbeg, jend;

  const double ambfac_max = Physics.params[i_param_ambfac_max];
 
  stride[0] = Grid.stride[0];
  stride[1] = Grid.stride[1];
  stride[2] = Grid.stride[2];

  const int loop_order[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};
  const double del[3][3]     = {{ 1.,0.,0.},{ 0.,1.,0.},{ 0.,0.,1.}};

  const double sign[3] = {-1.,1.,-1.};
  double signd1;

  int bounds[3][2];

  for(d=0;d<3;d++){
    bounds[d][0] = Grid.lbeg[d];
    bounds[d][1] = Grid.lend[d];
  }

  dxmin = Grid.dx[0];
  dxmax = Grid.dx[0]; 
  for(d=1;d<Grid.NDIM;d++){
    dxmin  = min(dxmin,Grid.dx[d]);
    dxmax  = max(dxmax,Grid.dx[d]);
  }

  double var[v_nvar][vsize], flx[v_nvar][vsize], res[v_nvar][vsize], curlBxB[3],
         var0,var1,var2,var3,var4,var5,var6,var7;
  double lf1[vsize], msx[vsize], msy[vsize], msz[vsize], lrt0,lrt1,lrt2, pres,
    sflx[vsize],fcond[vsize],curlB[3], bb[vsize],bbi, vv[vsize],vvi, vv_amb[vsize],
    v_amb[3][vsize], amb_fac[vsize], D_n[vsize], D_ni,temp[vsize], BgradT, hyp_spt;
   
  if(needs_curlB) {
#pragma acc parallel loop present(Grid[:1], Grid.curlB[:bufsize]) 
    for(i=0; i<bufsize; i++) {
      Grid.curlB[i].x = 0;
      Grid.curlB[i].y = 0;
      Grid.curlB[i].z = 0;
    }
  }

  if(spitzer) {
#pragma acc parallel loop present(Grid[:1], Grid.BgradT[:bufsize])
    for(i=0; i<bufsize; i++)
      Grid.BgradT[i] = 0;
  }

  if(ambipolar) {
#pragma acc parallel loop present(Grid[:1], Grid.curlBxB[:bufsize])
    for(i=0; i<bufsize; i++) {
      Grid.curlBxB[i].x = 0;
      Grid.curlBxB[i].y = 0;
      Grid.curlBxB[i].z = 0;
    }
  }

  // Estimate max characteristic velocity from previous time step (not defined at restart, set to zero) 
  if (Run.dt > 0.0) {
    cmax_hyp   = min(Run.CFL,CFL_hyp)*dxmin/Run.dt;
    nu_hyp_max = avr_max/Run.dt;
    hyp_diff   = eps_hyp/Run.dt;
  } else {
    cmax_hyp   = 0.0;
    nu_hyp_max = avr_max*cmax_hyp/(min(Run.CFL,CFL_hyp)*dxmin);
    hyp_diff   = eps_hyp*cmax_hyp/(min(Run.CFL,CFL_hyp)*dxmin);
  }

  double dt_diff_fac;
  if(Grid.NDIM == 1){
    dx=Grid.dx[0];
    dt_diff_fac = 0.5*dx*dx;
  } else if(Grid.NDIM == 2){
    dx=Grid.dx[0];
    dy=Grid.dx[1];
    dt_diff_fac = 0.5/(1./(dx*dx)+1./(dy*dy));
  } else {
    dx=Grid.dx[0];
    dy=Grid.dx[1];
    dz=Grid.dx[2]; 
    dt_diff_fac = 0.5/(1./(dx*dx)+1./(dy*dy)+1./(dz*dz));
  }
  
  cmax = 0.0;
  amb_diff_max = 0.0;
  amb_vel_max  = 0.0;
#pragma acc data  present(Grid[:1], U[:bufsize], R[:bufsize],                 \
   Grid.pres[:bufsize],   Grid.sflx[:bufsize],               \
   Grid.rhoi[:bufsize],   Grid.amb[:bufsize],                \
   Grid.v_amb[:bufsize],  Grid.tvar1[:bufsize],              \
   Grid.tvar2[:bufsize],  Grid.tvar3[:bufsize],              \
   Grid.Qamb[:bufsize],   Grid.curlBxB[:bufsize],            \
   Grid.BgradT[:bufsize], Grid.Rflx[:bufsize],               \
   Grid.curlB[:bufsize],  Grid.R_amb[:bufsize],              \
   Grid.temp[:bufsize]) 
{
#pragma acc loop seq independent
  for (d=0;d<Grid.NDIM;d++){
    d1=loop_order[d][0];
    d2=loop_order[d][1];
    d3=loop_order[d][2];

    dx = Grid.dx[d1];

    w1  = 8./(12.*dx);
    w2  =-1./(12.*dx);

    ww0 = -30./(12.*dx*dx);
    ww1 =  16./(12.*dx*dx);
    ww2 =  -1./(12.*dx*dx);

    p1 = del[d1][0];
    p2 = del[d1][1];
    p3 = del[d1][2];

    sz = bounds[d1][1]-bounds[d1][0]+1+2*Grid.ghosts[d1];
    gh = Grid.ghosts[d1];

    str = stride[d1];
    strd2 = stride[d2];
    strd3 = stride[d3];
    signd1 = sign[d1];
    kbeg = bounds[d3][0];
    kend = bounds[d3][1];
    jbeg = bounds[d2][0];
    jend = bounds[d2][1];
#pragma acc enter data copyin(d1,d2,d3,p1,p2,p3,w1,w2,ww0,ww1,ww2)
if(d == Grid.NDIM-1){

#pragma acc parallel loop collapse(2) gang                   \
 copy(cmax, amb_diff_max, amb_vel_max)                       \
 present(Grid[:1], U[:bufsize], R[:bufsize],                 \
   Grid.pres[:bufsize],   Grid.sflx[:bufsize],               \
   Grid.rhoi[:bufsize],   Grid.amb[:bufsize],                \
   Grid.v_amb[:bufsize],  Grid.tvar1[:bufsize],              \
   Grid.tvar2[:bufsize],  Grid.tvar3[:bufsize],              \
   Grid.Qamb[:bufsize],   Grid.curlBxB[:bufsize],            \
   Grid.BgradT[:bufsize], Grid.Rflx[:bufsize],               \
   Grid.curlB[:bufsize],  Grid.R_amb[:bufsize],              \
   Grid.temp[:bufsize])                                      \
 private(offset,                                             \
   res[:v_nvar][:vsize],               			     \
   flx[:v_nvar][:vsize], var[:v_nvar][:vsize],               \
   D_n[:vsize],v_amb[:3][:vsize],			                     \
   amb_fac[:vsize],               \
   vv_amb[:vsize], vv[:vsize], bb[:vsize],                   \
   sflx[:vsize],                              \
   msx[:vsize], msy[:vsize], msz[:vsize],                    \
   lf1[:vsize])                                              \
 reduction(max:amb_diff_max) reduction(max:amb_vel_max)      \
 reduction(max:cmax)
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
#pragma acc cache(d1,d2,d3,p1,p2,p3,w1,w2,ww0,ww1,ww2)
      offset = j*strd2+k*strd3;
      //time = MPI_Wtime();
      #pragma ivdep
#pragma acc loop vector private(node,rho_i,nu_in)
      for(i=0;i<sz;i++){
	node = offset+i*str;
	var[0][i] = U[node].d;
	var[1][i] = U[node].M.x;
	var[2][i] = U[node].M.y;
	var[3][i] = U[node].M.z;
	var[4][i] = U[node].e;
	var[5][i] = U[node].B.x;
	var[6][i] = U[node].B.y;	
	var[7][i] = U[node].B.z;
	
	sflx[i]   = 0.0;
	vv_amb[i] = 0.0;

      if(spitzer){
        #pragma ivdep
	  sflx[i] = Grid.sflx[node];
      }

      if(ambipolar){
        #pragma ivdep

	  rho_i = Grid.rhoi[node];
	  nu_in = Grid.amb[node];

	  D_n[i]= max(0.0,1.0-rho_i/var[0][i]);
	  
	  amb_fac[i] = 1.0/(rho_i*nu_in);
	  
	  // Switch off for T>1e4 as a saveguard for now
	  amb_fac[i] *= max(0.0,min(1.0,10.0-1e-3*Grid.temp[node]));
	  
	  amb_fac[i] = min(amb_fac[i], ambfac_max);

	  // test profile
	  //amb_fac[i] = 1e12*max(0.0,min(1.0,5.0-1e-3*temp[i]));
	  
	  v_amb[0][i] = Grid.v_amb[node].x*D_n[i];
	  v_amb[1][i] = Grid.v_amb[node].y*D_n[i];
	  v_amb[2][i] = Grid.v_amb[node].z*D_n[i];

	  vv_amb[i] = sqrt(v_amb[0][i]*v_amb[0][i]+v_amb[1][i]*v_amb[1][i]+v_amb[2][i]*v_amb[2][i]);
	}
      }      
      //r_time += MPI_Wtime()-time;

      //time = MPI_Wtime();
#pragma acc loop vector private(node)
      for(i=0;i<sz;i++){
        node = offset+i*str; 
	res[0][i] = 0.0;
 	res[1][i] = 0.0;
	res[2][i] = 0.0;
	res[3][i] = 0.0;
	res[4][i] = 0.0;
	res[5][i] = 0.0;
	res[6][i] = 0.0;
	res[7][i] = 0.0;
       }

#pragma acc loop vector private(var0,var1,var2,var3,var4,var5,var6,var7,node,\
 pres,vel, mag, dn, va2, ekin, pmag, x2, x4, s, cs2) reduction(max:cmax)
      for(i=0;i<sz;i++){
        node = offset+i*str;
        var0 = var[0][i];
 	var1 = var[1][i];
        var2 = var[2][i];
        var3 = var[3][i];
        var4 = var[4][i];
        var5 = var[5][i];
        var6 = var[6][i];
        var7 = var[7][i]; 
      
	vel = var[1+d1][i];
	mag = var[5+d1][i];

	dn   = 1.0/var0;
	
	vv[i] = var1*var1+var2*var2+var3*var3;
	bb[i] = var5*var5+var6*var6+var7*var7;

	va2  = bb[i]*dn;
	ekin = 0.5*var0*vv[i];
	pmag = 0.5*bb[i];

	// approximation of 1/sqrt(1+(va/c)^4)
	x2 = va2*inv_va2max;
	x4 = x2*x2;
	s  = 1.0+x2;
	lf1[i] = s/(s+x4);
        pres   = Grid.pres[node];
	cs2  = 5.0/3.0*pres*dn;

	vv[i] = sqrt(vv[i]);
	cmax = max(cmax,vv[i]+vv_amb[i]+sqrt(max(cs2,lf1[i]*(cs2+va2))));
	
	mflx = vel*var0;

	// Mass, momentum and energy fluxes without Maxwell stress
	flx[0][i] = mflx;
	flx[1][i] = mflx*var1 + p1*pres;
	flx[2][i] = mflx*var2 + p2*pres;
	flx[3][i] = mflx*var3 + p3*pres;
	flx[4][i] = vel*(var4 +ekin+pres);

	//Induction  equation
	flx[5][i] = vel*var5 - var1*mag;
	flx[6][i] = vel*var6 - var2*mag;
	flx[7][i] = vel*var7 - var3*mag;

	// Maxwell stress
	msx[i] = mag*var5 - p1*pmag;
	msy[i] = mag*var6 - p2*pmag;
	msz[i] = mag*var7 - p3*pmag;

      // Ambipolar induction
      if(ambipolar){
	   flx[5][i] += v_amb[d1][i]*var5 - v_amb[0][i]*var[5+d1][i];
	   flx[6][i] += v_amb[d1][i]*var6 - v_amb[1][i]*var[5+d1][i];
	   flx[7][i] += v_amb[d1][i]*var7 - v_amb[2][i]*var[5+d1][i];
      }

      // Heat conduction
      if(spitzer){
	  flx[4][i] += sflx[i]*var[5+d1][i];
        }
      } 	  
 #pragma acc loop vector collapse(2) 
      for (ivar=0;ivar<nvar;ivar++)
	for(i=gh;i<sz-gh;i++)
	  res[ivar][i] -= 
	    w1*(flx[ivar][i+1]-flx[ivar][i-1])+
	    w2*(flx[ivar][i+2]-flx[ivar][i-2]);

      if (need_diagnostics){
	if(spitzer){
#pragma acc loop vector private(node, dn) 
	  for(i=gh;i<sz-gh;i++){
	    node = offset+i*str;
            dn =-(w1*((sflx[i+1]*var[5+d1][i+1])-(sflx[i-1]*var[5+d1][i-1]))+w2*((sflx[i+2]*var[5+d1][i+2])-(sflx[i-2]*var[5+d1][i-2])));
	    Grid.tvar1[node] += res[4][i]-dn;
	    Grid.tvar2[node] += dn;
	  }
	} else {
#pragma acc loop vector private(node) 
	  for(i=gh;i<sz-gh;i++){
	    node = offset+i*str;
	    Grid.tvar1[node] += res[4][i];
	  }
	}
      }
      
      // Special treatment of Lorentz force for Boris Correction
#pragma acc loop vector private(dB2, dB3,lrt0,lrt1,lrt2,curlBxB[:3],curlB[:3]) 
      for(i=gh;i<sz-gh;i++){
#pragma acc cache(curlBxB[:3],curlB[:3])
        node = offset+i*str;
	dB2 = w1*(var[5+d2][i+1]-var[5+d2][i-1])+w2*(var[5+d2][i+2]-var[5+d2][i-2]);
	dB3 = w1*(var[5+d3][i+1]-var[5+d3][i-1])+w2*(var[5+d3][i+2]-var[5+d3][i-2]);

	// Lorentz force
	curlBxB[d1] =-var[5+d2][i]*dB2-var[5+d3][i]*dB3;
	curlBxB[d2] = var[5+d1][i]*dB2;
	curlBxB[d3] = var[5+d1][i]*dB3;

	lrt0 = w1*(msx[i+1]-msx[i-1])+w2*(msx[i+2]-msx[i-2]);
	lrt1 = w1*(msy[i+1]-msy[i-1])+w2*(msy[i+2]-msy[i-2]);
	lrt2 = w1*(msz[i+1]-msz[i-1])+w2*(msz[i+2]-msz[i-2]);

	// Transition from stress to jxB based on lfac
	lf=lf1[i];
	lrt0 = lf*lrt0 + (1.0-lf)*curlBxB[0];
	lrt1 = lf*lrt1 + (1.0-lf)*curlBxB[1];
	lrt2 = lf*lrt2 + (1.0-lf)*curlBxB[2];
	
	res[1][i] += lrt0;
	res[2][i] += lrt1;
	res[3][i] += lrt2;
	res[4][i] += var[1][i]*lrt0+var[2][i]*lrt1+var[3][i]*lrt2;

      // Ambipolar heating + hyperdiffusion
      if(ambipolar){
	  res[4][i] += v_amb[0][i]*lrt0+v_amb[1][i]*lrt1+v_amb[2][i]*lrt2;  
          Grid.curlBxB[node].x += lrt0;
          Grid.curlBxB[node].y += lrt1;
          Grid.curlBxB[node].z += lrt2;
      }

      if(needs_curlB){
          curlB[d1] = 0;
          curlB[d2] = signd1*dB3;
          curlB[d3] =-signd1*dB2;
          Grid.curlB[node].x += curlB[0];
          Grid.curlB[node].y += curlB[1];
          Grid.curlB[node].z += curlB[2];
          R[node].e += eta*(Grid.curlB[node].x*Grid.curlB[node].x+
                              Grid.curlB[node].y*Grid.curlB[node].y+
                              Grid.curlB[node].z*Grid.curlB[node].z);
      }

      if (need_diagnostics){
          Grid.tvar3[node] += var[1][i]*lrt0+var[2][i]*lrt1+var[3][i]*lrt2;
	if(ambipolar)
	    Grid.Qamb[node]  += v_amb[0][i]*lrt0+v_amb[1][i]*lrt1+v_amb[2][i]*lrt2;  
	  }
      }
      if(needs_curlB){
#pragma acc loop vector collapse(2) 
	for(ivar=5;ivar<=7;ivar++){
	  for(i=2;i<sz-2;i++){	    
	    res[ivar][i] += eta*(ww0*var[ivar][i]+
				 ww1*(var[ivar][i+1]+var[ivar][i-1])+
				 ww2*(var[ivar][i+2]+var[ivar][i-2])
				 );				 
	  }
	}	      
      }

      //c_time += MPI_Wtime()-time; 

      //time = MPI_Wtime();
      #pragma ivdep
#pragma acc loop vector private(node) 
      for(i=gh;i<sz-gh;i++){
	node = offset+i*str;
	R[node].d   += res[0][i];
	R[node].M.x += res[1][i];
	R[node].M.y += res[2][i];
	R[node].M.z += res[3][i];
	R[node].e   += res[4][i];
	R[node].B.x += res[5][i];
	R[node].B.y += res[6][i];
	R[node].B.z += res[7][i];

      // use this to temporarily store stuff while building up jxB

      }
      
      if(ambipolar){
#pragma acc loop vector private(node,amb_diff, c_lim, nu_hyp,tt32, kap_spt, \
  f_lim) \
 reduction(max:amb_diff_max) reduction(max:amb_vel_max) 
       for(i=gh;i<sz-gh;i++){
        node = offset+i*str;
      	// now curlBxB is fully computed, compute additional terms in Grid.R_amb
	    amb_diff     = max(1e-20,D_n[i]*amb_fac[i]*bb[i]);	    
	    amb_diff_max = max(amb_diff_max,amb_diff);

	    amb_vel_max = max(amb_vel_max,vv_amb[i]);
	    
	    c_lim  = cmax_hyp - vv[i];
	    f_lim  = 2.0*dt_diff_fac/(dxmin*dxmin);
	    nu_hyp = min(c_lim*c_lim*f_lim/amb_diff,nu_hyp_max);
	    
	    Grid.R_amb[node] += nu_hyp*(amb_fac[i]*Grid.curlBxB[node]-Grid.v_amb[node]);
      }
     }

      if(spitzer){
     // Last sweep, now j, jxB, BgradT are fully computed
#pragma acc loop vector private(node,c_lim,tt32, kap_spt, F_sat, \
 invBB, F_spt, f_lim, nu_hyp) 
       for(i=gh;i<sz-gh;i++){
            node = offset+i*str;
            int np1 = offset+(i+1)*str;
            int np2 = offset+(i+2)*str;
            int nm1 = offset+(i-1)*str;
            int nm2 = offset+(i-2)*str;
            double tempi = Grid.temp[node];
      	// now curlBxB is fully computed, compute additional terms in Grid.R_amb
            Grid.BgradT[node] += (var[5+d1][i]*(w1*(Grid.temp[np1]-Grid.temp[nm1])+w2*(Grid.temp[np2]-Grid.temp[nm2])));
            Grid.Rflx[node]   -= hyp_diff*((sflx[i+2]+sflx[i-2])-4.*(sflx[i+1]+sflx[i-1])+6.*sflx[i]);
 
	    c_lim    = cmax_hyp - vv[i];
	    tt32     = tempi*sqrt(tempi);
	    kap_spt  = k0_spt*tempi*tt32;
	    F_sat    = sat_flx*var[0][i]*tt32;	    
	    invBB    = 1.0/max(1e-100,bb[i]);
	    F_spt    =-kap_spt*Grid.BgradT[node];
	    f_lim    = F_sat/(F_sat+fabs(F_spt*sqrt(invBB)));
	    F_spt   *= invBB*f_lim;
	    kap_spt *= f_lim;
      
	    nu_hyp = min(c_lim*c_lim*var[4][i]/(tempi*kap_spt),nu_hyp_max);

	    Grid.Rflx[node] += nu_hyp*(F_spt-sflx[i]);
	  }
	}
      //r_time += MPI_Wtime()-time;
    }//OUTER_LOOP
}else{
#pragma acc parallel loop collapse(2) gang                   \
 copy(cmax)                       \
 present(Grid[:1], U[:bufsize], R[:bufsize],                 \
   Grid.pres[:bufsize],   Grid.sflx[:bufsize],               \
   Grid.rhoi[:bufsize],   Grid.amb[:bufsize],                \
   Grid.v_amb[:bufsize],  Grid.tvar1[:bufsize],              \
   Grid.tvar2[:bufsize],  Grid.tvar3[:bufsize],              \
   Grid.Qamb[:bufsize],   Grid.curlBxB[:bufsize],            \
   Grid.BgradT[:bufsize], Grid.Rflx[:bufsize],               \
   Grid.curlB[:bufsize],  Grid.R_amb[:bufsize],              \
   Grid.temp[:bufsize])                                      \
 private(offset,                                             \
   res[:v_nvar][:vsize],               			     \
   flx[:v_nvar][:vsize], var[:v_nvar][:vsize],               \
   v_amb[:3][:vsize],			                     \
   vv_amb[:vsize],                   \
   msx[:vsize], msy[:vsize], msz[:vsize],                    \
   lf1[:vsize])                                              \
 reduction(max:cmax)
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {
#pragma acc cache(d1,d2,d3,p1,p2,p3,w1,w2,ww0,ww1,ww2)
      offset = j*strd2+k*strd3;
      //time = MPI_Wtime();
      #pragma ivdep
#pragma acc loop vector private(node,D_ni)
      for(i=0;i<sz;i++){
	node = offset+i*str;
	var[0][i] = U[node].d;
	var[1][i] = U[node].M.x;
	var[2][i] = U[node].M.y;
	var[3][i] = U[node].M.z;
	var[4][i] = U[node].e;
	var[5][i] = U[node].B.x;
	var[6][i] = U[node].B.y;	
	var[7][i] = U[node].B.z;

//	sflx[i]   = 0.0;
	vv_amb[i] = 0.0;

      if(ambipolar){
        #pragma ivdep

	  //rho_i = Grid.rhoi[node];
	  //nu_in = Grid.amb[node];

	  D_ni = max(0.0,1.0-Grid.rhoi[node]/var[0][i]);
	  
	  //amb_fac[i] = 1.0/(rho_i*nu_in);
	  
	  // Switch off for T>1e4 as a saveguard for now
	  //amb_fac[i] *= max(0.0,min(1.0,10.0-1e-3*Grid.temp[node]));
	  
	  //amb_fac[i] = min(amb_fac[i], ambfac_max);

	  // test profile
	  //amb_fac[i] = 1e12*max(0.0,min(1.0,5.0-1e-3*temp[i]));
	  
	  v_amb[0][i] = Grid.v_amb[node].x*D_ni;
	  v_amb[1][i] = Grid.v_amb[node].y*D_ni;
	  v_amb[2][i] = Grid.v_amb[node].z*D_ni;

	  vv_amb[i] = sqrt(v_amb[0][i]*v_amb[0][i]+v_amb[1][i]*v_amb[1][i]+v_amb[2][i]*v_amb[2][i]);
	}
      }      
      //r_time += MPI_Wtime()-time;

      //time = MPI_Wtime();
#pragma acc loop vector private(node)
      for(i=0;i<sz;i++){
        node = offset+i*str; 
	res[0][i] = 0.0;
 	res[1][i] = 0.0;
	res[2][i] = 0.0;
	res[3][i] = 0.0;
	res[4][i] = 0.0;
	res[5][i] = 0.0;
	res[6][i] = 0.0;
	res[7][i] = 0.0;
       }

#pragma acc loop vector private(var0,var1,var2,var3,var4,var5,var6,var7,node,\
 pres,vel, mag, dn, va2, ekin, pmag, x2, x4, s, cs2,bbi,vvi) reduction(max:cmax)
      for(i=0;i<sz;i++){
        node = offset+i*str;
        var0 = var[0][i];
 	var1 = var[1][i];
        var2 = var[2][i];
        var3 = var[3][i];
        var4 = var[4][i];
        var5 = var[5][i];
        var6 = var[6][i];
        var7 = var[7][i]; 
      
	vel = var[1+d1][i];
	mag = var[5+d1][i];

	dn   = 1.0/var0;
	
	vvi= var1*var1+var2*var2+var3*var3;
	bbi = var5*var5+var6*var6+var7*var7;

	va2  = bbi*dn;
	ekin = 0.5*var0*vvi;
	pmag = 0.5*bbi;

	// approximation of 1/sqrt(1+(va/c)^4)
	x2 = va2*inv_va2max;
	x4 = x2*x2;
	s  = 1.0+x2;
	lf1[i] = s/(s+x4);
        pres   = Grid.pres[node];
	cs2  = 5.0/3.0*pres*dn;

	vvi = sqrt(vvi);
	cmax = max(cmax,vvi+vv_amb[i]+sqrt(max(cs2,lf1[i]*(cs2+va2))));
	
	mflx = vel*var0;

	// Mass, momentum and energy fluxes without Maxwell stress
	flx[0][i] = mflx;
	flx[1][i] = mflx*var1 + p1*pres;
	flx[2][i] = mflx*var2 + p2*pres;
	flx[3][i] = mflx*var3 + p3*pres;
	flx[4][i] = vel*(var4 +ekin+pres);

	//Induction  equation
	flx[5][i] = vel*var5 - var1*mag;
	flx[6][i] = vel*var6 - var2*mag;
	flx[7][i] = vel*var7 - var3*mag;

	// Maxwell stress
	msx[i] = mag*var5 - p1*pmag;
	msy[i] = mag*var6 - p2*pmag;
	msz[i] = mag*var7 - p3*pmag;

      // Ambipolar induction
      if(ambipolar){
	   flx[5][i] += v_amb[d1][i]*var5 - v_amb[0][i]*var[5+d1][i];
	   flx[6][i] += v_amb[d1][i]*var6 - v_amb[1][i]*var[5+d1][i];
	   flx[7][i] += v_amb[d1][i]*var7 - v_amb[2][i]*var[5+d1][i];
      }

      // Heat conduction
      if(spitzer){
	  flx[4][i] += Grid.sflx[node]*var[5+d1][i];
        }
      } 	  
 #pragma acc loop vector collapse(2) 
      for (ivar=0;ivar<nvar;ivar++)
	for(i=gh;i<sz-gh;i++)
	  res[ivar][i] -= 
	    w1*(flx[ivar][i+1]-flx[ivar][i-1])+
	    w2*(flx[ivar][i+2]-flx[ivar][i-2]);

      if(spitzer){
#pragma acc loop vector private(node, dn) 
          for(i=gh;i<sz-gh;i++){
            node = offset+i*str;
            int np1 = offset+(i+1)*str;
            int np2 = offset+(i+2)*str;
            int nm1 = offset+(i-1)*str;
            int nm2 = offset+(i-2)*str;
            Grid.Rflx[node]   -= hyp_diff*((Grid.sflx[np2]+Grid.sflx[nm2])-4.*(Grid.sflx[np1]+Grid.sflx[nm1])+6.*Grid.sflx[node]);
            Grid.BgradT[node] += (var[5+d1][i]*(w1*(Grid.temp[np1]-Grid.temp[nm1])+w2*(Grid.temp[np2]-Grid.temp[nm2])));
            if(need_diagnostics){
              dn =-(w1*((Grid.sflx[np1]*var[5+d1][i+1])-(Grid.sflx[nm1]*var[5+d1][i-1]))
                   +w2*((Grid.sflx[np2]*var[5+d1][i+2])-(Grid.sflx[nm2]*var[5+d1][i-2])));
              Grid.tvar1[node] += res[4][i]-dn;
              Grid.tvar2[node] += dn;
            }
          }

      }else if(need_diagnostics){
#pragma acc loop vector private(node) 
	  for(i=gh;i<sz-gh;i++){
	    node = offset+i*str;
	    Grid.tvar1[node] += res[4][i];
	  }
	}
      
      // Special treatment of Lorentz force for Boris Correction
#pragma acc loop vector private(dB2, dB3,lrt0,lrt1,lrt2,curlBxB[:3],curlB[:3]) 
      for(i=gh;i<sz-gh;i++){
#pragma acc cache(curlBxB[:3],curlB[:3])
        node = offset+i*str;
	dB2 = w1*(var[5+d2][i+1]-var[5+d2][i-1])+w2*(var[5+d2][i+2]-var[5+d2][i-2]);
	dB3 = w1*(var[5+d3][i+1]-var[5+d3][i-1])+w2*(var[5+d3][i+2]-var[5+d3][i-2]);

	// Lorentz force
	curlBxB[d1] =-var[5+d2][i]*dB2-var[5+d3][i]*dB3;
	curlBxB[d2] = var[5+d1][i]*dB2;
	curlBxB[d3] = var[5+d1][i]*dB3;

	// curlB
	lrt0 = w1*(msx[i+1]-msx[i-1])+w2*(msx[i+2]-msx[i-2]);
	lrt1 = w1*(msy[i+1]-msy[i-1])+w2*(msy[i+2]-msy[i-2]);
	lrt2 = w1*(msz[i+1]-msz[i-1])+w2*(msz[i+2]-msz[i-2]);

	// Transition from stress to jxB based on lfac
	lf=lf1[i];
	lrt0 = lf*lrt0 + (1.0-lf)*curlBxB[0];
	lrt1 = lf*lrt1 + (1.0-lf)*curlBxB[1];
	lrt2 = lf*lrt2 + (1.0-lf)*curlBxB[2];
	
	res[1][i] += lrt0;
	res[2][i] += lrt1;
	res[3][i] += lrt2;
	res[4][i] += var[1][i]*lrt0+var[2][i]*lrt1+var[3][i]*lrt2;

      // Ambipolar heating + hyperdiffusion
      if(ambipolar){
	  res[4][i] += v_amb[0][i]*lrt0+v_amb[1][i]*lrt1+v_amb[2][i]*lrt2;  
          Grid.curlBxB[node].x += lrt0;
          Grid.curlBxB[node].y += lrt1;
          Grid.curlBxB[node].z += lrt2;
      }

      if(needs_curlB){
          curlB[d1] = 0;
          curlB[d2] = signd1*dB3;
          curlB[d3] =-signd1*dB2;
          Grid.curlB[node].x += curlB[0];
          Grid.curlB[node].y += curlB[1];
          Grid.curlB[node].z += curlB[2];
      }

      if (need_diagnostics){
          Grid.tvar3[node] += var[1][i]*lrt0+var[2][i]*lrt1+var[3][i]*lrt2;
	if(ambipolar)
	    Grid.Qamb[node]  += v_amb[0][i]*lrt0+v_amb[1][i]*lrt1+v_amb[2][i]*lrt2;  
	  }
      }
      if(needs_curlB){
#pragma acc loop vector collapse(2) 
	for(ivar=5;ivar<=7;ivar++){
	  for(i=2;i<sz-2;i++){	    
	    res[ivar][i] += eta*(ww0*var[ivar][i]+
				 ww1*(var[ivar][i+1]+var[ivar][i-1])+
				 ww2*(var[ivar][i+2]+var[ivar][i-2])
				 );				 
	  }
	}	      
      }

      //c_time += MPI_Wtime()-time; 

      //time = MPI_Wtime();
      #pragma ivdep
#pragma acc loop vector private(node) 
      for(i=gh;i<sz-gh;i++){
	node = offset+i*str;
	R[node].d   += res[0][i];
	R[node].M.x += res[1][i];
	R[node].M.y += res[2][i];
	R[node].M.z += res[3][i];
	R[node].e   += res[4][i];
	R[node].B.x += res[5][i];
	R[node].B.y += res[6][i];
	R[node].B.z += res[7][i];

      }
      

      //r_time += MPI_Wtime()-time;
    }//OUTER_LOOP

    } // d!= ndim-1
  }//d
 } //exit data region

  dt_cfl = dxmin/cmax;

  if(ambipolar){
    dt_amb = dt_diff_fac/amb_diff_max;
    if(Run.CFL > 1) dt_amb *= 0.5;
    dt_cfl = min(dt_cfl,dt_amb*Physics.params[i_param_ambipolar]);
    
    if( !(call_count%maxstage) && (Run.verbose > 3) ){
      double sbuf[2],rbuf[2];
      sbuf[0] = amb_diff_max;
      sbuf[1] = amb_vel_max;
      MPI_Reduce(sbuf,rbuf,2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      if( Run.rank == 0 ){
	cout << " Amb_diff_max = " << rbuf[0]      << " cm^2/s"         << endl;
	cout << " Amb_vel_max  = " << rbuf[1]*1e-5 << " km/s"           << endl;
	cout << " Amb_speedup  = " << rbuf[0]*Run.dt/dt_diff_fac/Run.CFL << endl;;
      }
    }
  }
  
  if (needs_curlB) {
    dt_res = dt_diff_fac/eta;
    if(Run.CFL > 1) dt_res *= 0.5;
    dt_cfl = min(dt_cfl,dt_res);
  }

  //t_time += MPI_Wtime()-s_time;
  //call_count += 1;

  //if(Run.rank == 0)
  //cout << "MHD: " << t_time/call_count << ' ' 
  // << c_time/call_count << ' ' 
  // << r_time/call_count << endl;

  PGI_COMPARE(Grid.U, double, Grid.bufsize*8, "U", "mhdres_SR.C", "MHD", 1)
  if(needs_curlB) {
    PGI_COMPARE(Grid.curlB, double, Grid.bufsize*3, "curlB", "mhdres_SR.C", "MHD", 2)
  }
  if(spitzer) {
    PGI_COMPARE(Grid.BgradT, double, Grid.bufsize, "BgradT", "mhdres_SR.C", "MHD", 3)
    PGI_COMPARE(Grid.Rflx, double, Grid.bufsize, "Rflx", "mhdres_SR.C", "MHD", 4)
  }
  if(ambipolar) {
    PGI_COMPARE(Grid.curlBxB, double, Grid.bufsize*3, "curlBxB", "mhdres_SR.C", "MHD", 5)
    PGI_COMPARE(Grid.R_amb, double, Grid.bufsize*3, "R_amb", "mhdres_SR.C", "MHD", 6)
  }
  if(need_diagnostics) {
    PGI_COMPARE(Grid.tvar1, double, Grid.bufsize, "tvar1", "mhdres_SR.C", "MHD", 7)
    PGI_COMPARE(Grid.tvar2, double, Grid.bufsize, "tvar2", "mhdres_SR.C", "MHD", 8)
    PGI_COMPARE(Grid.tvar3, double, Grid.bufsize, "tvar3", "mhdres_SR.C", "MHD", 9)
    if(ambipolar) {
      PGI_COMPARE(Grid.Qamb, double, Grid.bufsize, "Qamb", "mhdres_SR.C", "MHD", 10)
    }
  }
  PGI_COMPARE(&dt_cfl, double, 1, "dt_cfl", "mhdres_SR.C", "MHD", 11)

  call_count += 1;

  return dt_cfl;
}







