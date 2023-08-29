#include <mpi.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include <cmath>
#include <stdlib.h>
#include "eos.H"
#include "limit_va.H"
#include "muramacc.H"
#include <algorithm>
#include <iostream>

using std::min;
using std::max;
using std::cout;
using std::endl;
/*
   - Special Corona version. Allows to set numerical pm in Corona. 
   - Resistive heating in top layer can be omitted.
   - Enhance viscosity for v>vmax_lim.
   - Output of Qres and Qvis for diagnostic reasons.
*/ 

double * hfb;
double * hft;
double * qft;

#define OUTER_LOOP(G,i,j,d1,d2) \
  for((j)=(G)[(d2)][0];(j)<=(G)[(d2)][1];(j)++) \
  for((i)=(G)[(d1)][0];(i)<=(G)[(d1)][1];(i)++)

/********************************************************************/
void TVDlimit(const RunData&  Run, GridData& Grid, 
	      const PhysicsData& Physics, const double dt_tvd){

  //double time,s_time;
  //static double t_time = 0.0 ,c_time = 0.0 , r_time = 0.0;
  //static int call_count = 0;
  //s_time = MPI_Wtime();

  const double rho_min = Physics.tchk[i_tchk_rho_min];
  const double eps_min = Physics.tchk[i_tchk_eps_min];

  const double* tvd_h    = Physics.tvd_h;
  const double* tvd_cs   = Physics.tvd_cs;  
  const double rho_lev   = pow(Physics.tvd[i_tvd_rholev],2);
  const int    rho_log   = (int) Physics.tvd[i_tvd_rholog];
  const double q_rho_max = log(Physics.tvd[i_tvd_qrho]);

  const double visc_coeff_bot = min(1.0,Physics.tvd_visc_bnd[0]);
  const double visc_slope_bot = max(0.0,Physics.tvd_visc_bnd[0]-1.0);
  const double visc_coeff_top = min(1.0,Physics.tvd_visc_bnd[1]);
  const double visc_slope_top = max(0.0,Physics.tvd_visc_bnd[1]-1.0);

  const double eta_coeff_bot = min(1.0,Physics.tvd_eta_bnd[0]);
  const double eta_slope_bot = max(0.0,Physics.tvd_eta_bnd[0]-1.0);
  const double eta_coeff_top = min(1.0,Physics.tvd_eta_bnd[1]);
  const double eta_slope_top = max(0.0,Physics.tvd_eta_bnd[1]-1.0);

  const double B_par_diff = Physics.tvd[i_tvd_Bpar];
  const double vhyp       = Physics.tvd[i_tvd_vhyp]; // Additional Hyperdiffusion in vertical direction
  const double Qdiff_bnd  = Physics.tvd[i_tvd_Qdiff_bnd]; // Value of Qdiff at upper boundary
  const double tvd_pm_v   = Physics.tvd[i_tvd_pm_v]; //
  const double tvd_pm_B   = Physics.tvd[i_tvd_pm_B]; //
  const double vmax_lim   = Physics.tvd[i_tvd_vmax_lim]; // (<1) relative to vlim, (>1) absolut, (0) disable
  const double CME_thresh = Physics.tvd[i_tvd_CME_thresh];

  const int nvar = Physics.NVAR;
  const int vsize = Grid.vsize;
  const int v_nvar = Grid.v_nvar;
  const int bufsize = Grid.bufsize;

  double h_bnd;
  if(Physics.tvd_h_bnd[0] < 1)
    h_bnd = Physics.tvd_h_bnd[0];
  else
    h_bnd = Physics.tvd_h_bnd[0]/double(Grid.gsize[0]);

  const double h_bnd_bot = h_bnd;

  if(Physics.tvd_h_bnd[1] < 1)
    h_bnd = Physics.tvd_h_bnd[1];
  else
    h_bnd = Physics.tvd_h_bnd[1]/double(Grid.gsize[0]);

  const double h_bnd_top = h_bnd;

  static int tvd_ini_flag = 1;

  const bool need_diagnostics = Run.need_diagnostics;
  const bool ambipolar        = (Physics.params[i_param_ambipolar] > 0.0);

  register int i,j,k,node,d,d1,d2,d3,ivar,offset,str,istart;
  int sz,gh;

  int kbeg, kend, jbeg, jend;

  cState* U  = (cState*) Grid.U;
    
  int stride[3];
  int strided1,strided2,strided3;
  //int dim_order[3];
  int dim_order0, dim_order1, dim_order2;
 
  int bounds[3][2];

  //dim_order[0] = 0;
  //dim_order[1] = 1;
  //dim_order[2] = 2;
  dim_order0 = 0;
  dim_order1 = 1;
  dim_order2 = 2;
  
  stride[0] = Grid.stride[0];
  stride[1] = Grid.stride[1];
  stride[2] = Grid.stride[2];

  const int loop_order[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  double dn,vv,bb,cs2,va2,lf,rf,hh,cf,cmax,sl,sr,sl_lim,vsqr_diff,lower,upper,
    slm,x2,x4,s,cfast,CME_mode;

  double idx[3],idxd1;
  double tvd_coeff0,tvd_coeff1,tvd_coeff2,tvd_coeff3;

  double idt_full = 1.0/Run.dt;
  double dt_fac = dt_tvd/Run.dt;

  idx[0] = 1./Grid.dx[0];
  if (Grid.NDIM >= 2) idx[1] = 1./Grid.dx[1];
  if (Grid.NDIM == 3) idx[2] = 1.0/Grid.dx[2];
  tvd_coeff0 = Physics.tvd_coeff[0];
  tvd_coeff1 = Physics.tvd_coeff[1];
  tvd_coeff2 = Physics.tvd_coeff[2];
  tvd_coeff3 = Physics.tvd_coeff[3];

  int need_tvd_coeff;
  if( tvd_coeff0*tvd_coeff1*tvd_coeff2*tvd_coeff3 != 1.0 )
    need_tvd_coeff = 1;
  else
    need_tvd_coeff = 0;

  if (tvd_ini_flag == 1){
    hfb = new double [vsize];
    qft = new double [vsize];
    hft = new double [vsize];

     // parabolic profile, zero for h>2*h_bnd
    for(i=Grid.lbeg[0]-Grid.ghosts[0];i<=Grid.lend[0]+Grid.ghosts[0];i++) {
      hh      = Grid.coord(i,0)/Grid.gxmax[0];
      vv = max(0.0,1.0-0.5*(1.0-hh)/h_bnd_top);
      hft[i]  = min(1.0,vv*vv);  
      vv = max(0.0,1.0-0.5*hh/h_bnd_bot);     
      hfb[i]  = min(1.0,vv*vv);

      if(hh < 1-2*h_bnd_top)
	qft[i] = 1.0;
      else
	qft[i] = Qdiff_bnd;
    }

    if (Run.rank == 0) {
      cout << " *** SPLIT TVD Coronal VERSION *** "           << endl;    
      if(rho_log){
	cout << " *** use log(rho) and log(eps) *** "    << endl;
      }
      cout << "tvd: apply additional 4th order hyperdiff to 0,2,4" << endl;
      cout << "tvd: special CME mode for v > CME_thresh*c_fast   " << endl;
      cout << "tvd: disable mass diffusion correction in corona  " << endl;
      cout << "tvd: change diffusivity based on lfac and rho     " << endl; 
      cout << "tvd: h              = "  << tvd_h[0] << ' ' 
	   << tvd_h[1] << ' ' << tvd_h[2] << ' ' << tvd_h[3] << endl; 
      cout << "tvd: cs             = "  << tvd_cs[0]  << ' ' 
	   << tvd_cs[1] << ' ' << tvd_cs[2] << ' ' << tvd_cs[3] << endl;
      cout << "tvd: tvd_coeff      = "
           << tvd_coeff0 << ' ' << tvd_coeff1 << ' '
           << tvd_coeff2 << ' ' << tvd_coeff3 << endl;
      cout << "tvd: rho_lev        = "  << sqrt(rho_lev)  << endl; 
      cout << "tvd: rho_log        = "  << rho_log        << endl;
      cout << "tvd: q_rho_max      = "  << exp(q_rho_max) << endl;
      cout << "tvd: vhyp           = "  << vhyp         << endl;
      cout << "tvd: B_par_diff     = "  << B_par_diff   << endl;
      cout << "tvd: Qdiff_bnd      = "  << Qdiff_bnd    << endl;
      cout << "tvd: tvd_pm_v       = "  << tvd_pm_v     << endl;
      cout << "tvd: tvd_pm_B       = "  << tvd_pm_B     << endl;
      cout << "tvd: vmax_lim       = "  << vmax_lim     << endl;
      cout << "tvd: CME_thresh     = "  << CME_thresh   << endl;
      cout << "tvd: h_bnd          = " <<  h_bnd_bot << ' ' << h_bnd_top << endl;  
      cout << "tvd: visc_coeff_bnd = " <<  visc_coeff_bot << ' '
	   << visc_coeff_top  << endl;
      cout << "tvd: visc_slope_bnd = " <<  visc_slope_bot << ' '
	   << visc_slope_top << endl;
      cout << "tvd: eta_coeff_bnd  = " <<  eta_coeff_bot << ' '
	   << eta_coeff_top << endl;
      cout << "tvd: eta_slope_bnd  = " <<  eta_slope_bot << ' ' 
	   << eta_slope_top << endl;     
    }
#pragma acc enter data copyin(hfb[:vsize])
#pragma acc enter data copyin(qft[:vsize])
#pragma acc enter data copyin(hft[:vsize])
#pragma acc enter data copyin(tvd_h[:4])
#pragma acc enter data copyin(tvd_cs[:4])
#pragma acc enter data copyin(loop_order[:3][:3])
//#pragma acc enter data copyin(stride[:3])      
//#pragma acc enter data copyin(idx[:3])   
    tvd_ini_flag =0;
  }
  if(vmax_lim < 1.0)
    vsqr_diff = pow(vmax_lim*v_lim,2.0);
  else
    vsqr_diff = pow(vmax_lim,2.0);

  sz = 1;
  for(d=0;d<3;d++){
    bounds[d][0] = Grid.lbeg[d]-Grid.ghosts[d];
    bounds[d][1] = Grid.lend[d]+Grid.ghosts[d];

    sz *= Grid.lsize[d]+2*Grid.ghosts[d];
  }
  
  double var[v_nvar][vsize], slp[v_nvar][vsize], res[v_nvar][vsize], 
    flx[v_nvar][vsize], boris0,boris1,boris2,boris3, vv_amb,
    flx0,flx1,flx2,flx3,flx4,uif0,uif1,uif2,uif3,uif4,min_cz;
  double qrho[vsize], hft1[vsize], hfb1[vsize], cm[vsize],
    tvd_faci, bsqr[vsize], vsqr[vsize], BC2[vsize],
    hyperdiff[vsize], qres, qft1[vsize],
    hyp_e[vsize],hyp_v[vsize],hyp_B[vsize],CZ_fac[vsize];
#pragma acc data present(Grid, U[:bufsize], Grid.pres[:bufsize],                   \
    Grid.Tau[:bufsize], Grid.v_amb[:bufsize],                       \
    hfb[:vsize], hft[:vsize], qft[:vsize],                          \
    Grid.tvar8[:bufsize], Grid.Qres[:bufsize],                      \
    Grid.Qvis[:bufsize], Grid.tvar6[:bufsize],                      \
    Grid.tvar7[:bufsize],tvd_cs[:4],tvd_h[:4],loop_order[:3][:3]) //tvd_coeff[:4])  
{
  /* y direction first to be consistent with vertical boundary */
#pragma acc loop seq independent
  for (d=0;d<Grid.NDIM;d++){
    d1=loop_order[d][0];
    d2=loop_order[d][1];
    d3=loop_order[d][2];
    cmax = 0.975*Grid.dx[d1]/dt_tvd;

    sz = bounds[d1][1]-bounds[d1][0]+1;
    gh = Grid.ghosts[d1];

    str=stride[d1];
    strided1 = stride[d1];
    strided2 = stride[d2];
    strided3 = stride[d3];
    idxd1 = idx[d1];
    kbeg = bounds[d3][0];
    kend = bounds[d3][1];
    jbeg = bounds[d2][0];
    jend = bounds[d2][1];

#pragma acc parallel loop collapse(2) gang                          \
  present(Grid, U[:bufsize], Grid.pres[:bufsize],                   \
    Grid.Tau[:bufsize], Grid.v_amb[:bufsize],                       \
    hfb[:vsize], hft[:vsize], qft[:vsize],                          \
    Grid.tvar8[:bufsize], Grid.Qres[:bufsize],                      \
    Grid.Qvis[:bufsize], Grid.tvar6[:bufsize],                      \
    Grid.tvar7[:bufsize],tvd_cs[:4],tvd_h[:4])                      \
  private(offset, var[:v_nvar][:vsize],                             \
    slp[:v_nvar][:vsize],                                           \
    CZ_fac[:vsize], 		                 		    \
    hfb1[:vsize], hft1[:vsize], qft1[:vsize],                       \
    vsqr[:vsize], bsqr[:vsize], cm[:vsize],                         \
    hyp_e[:vsize], hyp_v[:vsize], hyp_B[:vsize],                    \
    BC2[:vsize], hyperdiff[:vsize], qrho[:vsize],                   \
    flx[:v_nvar][:vsize], 		                            \
    res[:v_nvar][:vsize])				            \
    firstprivate(istart)
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++) {

      offset = j*strided2+k*strided3;
      //time = MPI_Wtime();
      #pragma ivdep
#pragma acc loop vector private(node)
      for(i=0;i<sz;i++){
	node = offset+i*str;
	var[0][i] = max(rho_min,U[node].d);
	var[1][i] = U[node].M.x;
	var[2][i] = U[node].M.y;
	var[3][i] = U[node].M.z;
	var[4][i] = U[node].e;
	var[5][i] = U[node].B.x;
	var[6][i] = U[node].B.y;	
	var[7][i] = U[node].B.z;
	
	//pres[i] = Grid.pres[node];
	CZ_fac[i] = (double) (Grid.Tau[node] > 1.0e-5);

	//vv_amb[i] = 0.0;
	if(d1 == 0){
          hfb1[i] = hfb[i];
          hft1[i] = hft[i];
          qft1[i] = qft[i];
        } else if (d2 == 0){
          hfb1[i] = hfb[j];
          hft1[i] = hft[j];
          qft1[i] = qft[j];
 
         } else if(d3 == 0){
          hfb1[i] = hfb[k];
          hft1[i] = hft[k];
          qft1[i] = qft[k];
        }
      }
      //r_time += MPI_Wtime()-time;

      //time = MPI_Wtime();
#pragma acc loop vector private(node,vv_amb,dn, vv, va2, x2, x4, s, lf, cs2, \
 cfast, CME_mode, rf, hh, cf,dn, vv, bb) 
      for(i=0;i<sz;i++){
        node = offset+i*str;
	dn        = 1.0/var[0][i];
        var[1][i] = var[1][i]*dn;
	var[2][i] = var[2][i]*dn;
	var[3][i] = var[3][i]*dn;
	
	vv        = (var[1][i]*var[1][i]+
		     var[2][i]*var[2][i]+
		     var[3][i]*var[3][i]);
	     
	bb        = (var[5][i]*var[5][i]+
		     var[6][i]*var[6][i]+
		     var[7][i]*var[7][i]);

	var[4][i] = var[4][i]*dn-0.5*vv;

	vsqr[i]   = vv;
	bsqr[i]   = bb;
	
	vv   = sqrt(vv);
	va2  = bb*dn;
	
	x2 = va2*inv_va2max;
	x4 = x2*x2;
	s  = 1.0+x2;
	lf = s/(s+x4);

	cs2  = 5.0/3.0*Grid.pres[node]*dn;

	cfast = sqrt(va2+cs2);

	CME_mode = (1.0-CZ_fac[i])*( (double) (vv >= CME_thresh*cfast) );
	
	dn = var[0][i]*var[0][i];
	rf = min(lf,dn/(dn+rho_lev));
	
	hh   = 1.0-hfb1[i]-hft1[i];
	cf   = hfb1[i]*tvd_cs[0]+(rf*tvd_cs[1]+(1.0-rf)*tvd_cs[2])*hh+hft1[i]*tvd_cs[3];

	cf   = max(cf,CME_mode);
	
	cs2 *= cf*cf;
        if(ambipolar){
          vv_amb = sqrt(Grid.v_amb[node].x*Grid.v_amb[node].x+Grid.v_amb[node].y*Grid.v_amb[node].y+Grid.v_amb[node].z*Grid.v_amb[node].z);
        }else{
          vv_amb = 0.0;
        }	
	cm[i]   = vv_amb+vv+sqrt(max(cs2,lf*(cs2+va2)));
	
	hyp_e[i]  = hfb1[i]*tvd_h[0]+(rf*tvd_h[1]+(1.0-rf)*tvd_h[2])*hh+hft1[i]*tvd_h[3];
	hyp_v[i]  = hfb1[i]*tvd_h[0]+(rf*tvd_h[1]+(1.0-rf)*tvd_h[2]*tvd_pm_v)*hh+hft1[i]*tvd_h[3]*tvd_pm_v;
	hyp_B[i]  = hfb1[i]*tvd_h[0]+(rf*tvd_h[1]+(1.0-rf)*tvd_h[2]*tvd_pm_B)*hh+hft1[i]*tvd_h[3]*tvd_pm_B;

	hyp_e[i] *= (1.0-CME_mode);
	hyp_v[i] *= (1.0-CME_mode);
	
	BC2[i]  = 1.0-lf;
      }
      
      if( (d1 == 0) && (vhyp > 0.0) ){
#pragma acc loop vector 
	for(i=0;i<sz;i++){
	  hyperdiff[i] = vhyp*sqrt(vsqr[i]);	  
	  hyperdiff[i] = 0.5*idxd1*max(0.0,min(hyperdiff[i],cmax-cm[i]));
	}
      }

      if(rho_log){
#pragma acc loop vector 
	for(i=0;i<sz;i++){
	  var[0][i] = log(var[0][i]);
	  var[4][i] = log(var[4][i]);
	}
#pragma acc loop vector 
	for(i=0;i<sz-1;i++)
	  qrho[i] = fabs(var[0][i]-var[0][i+1]);
      } else {
#pragma acc loop vector 
	for(i=0;i<sz-1;i++)
	  qrho[i] = fabs(log(var[0][i])-log(var[0][i+1])); 
      }
	
      // reconstruction slopes
#pragma acc loop vector collapse(2) private(sl, sr, slm, lower, upper)
      for(ivar=0;ivar<nvar;ivar++){
	for(i=1;i<sz-1;i++){
	  sl    = var[ivar][i]-var[ivar][i-1];
	  sr    = var[ivar][i+1]-var[ivar][i];
	  slm   = 0.25*(sl+sr);
	  lower = min(min(sl,sr),slm);
	  upper = max(max(sl,sr),slm);
	  slp[ivar][i]=max(0.0,lower)+min(0.0,upper);
	}
      }

      // more diffusivity at vertical boundaries
      if( (visc_slope_bot > 0.0) or ( visc_slope_top > 0.0) ){
#pragma acc loop vector private(cf)
	for(i=1;i<sz-1;i++){
	  cf =(1.0-visc_slope_bot*hfb1[i])*(1.0-visc_slope_top*hft1[i]);
	  slp[1][i] *= cf; slp[2][i] *= cf; slp[3][i] *= cf;	
	}
      }

      if( (eta_slope_bot > 0.0) or ( eta_slope_top > 0.0) ){
#pragma acc loop vector private(cf) 
	for(i=1;i<sz-1;i++){
	  cf =(1.0-eta_slope_bot*hfb1[i])*(1.0-eta_slope_top*hft1[i]);
	  slp[5][i] *= cf; slp[6][i] *= cf; slp[7][i] *= cf;	
	}
      }

      // more diffusivity around large density jumps
#pragma acc loop vector 
      for(i=1;i<sz-1;i++){
	if(qrho[i] > q_rho_max)
	  slp[0][i] = 0.0;
      }

      // more diffusivity for very fast flows
      if(vmax_lim > 0.0){
#pragma acc loop vector 
	for(i=1;i<sz-1;i++){
	  if(vsqr[i] > vsqr_diff){
	    slp[1][i] = 0.0;
	    slp[2][i] = 0.0;
	    slp[3][i] = 0.0;
	  }
	}
      }

#pragma acc loop vector collapse(2) private(sl, sl_lim, dn, rf, cf) 
      for(ivar=0;ivar<nvar;ivar++){
        for(i=1;i<sz-2;i++){
	    sl     = var[ivar][i+1]-var[ivar][i];
	    sl_lim = var[ivar][i+1]-var[ivar][i]-(slp[ivar][i]+slp[ivar][i+1]);
            if( (ivar >=1) && (ivar <= 3) ){
	      dn     = hyp_v[i];
            }else if( (ivar >=5) && (ivar <= 7) ){
              dn     = hyp_B[i];
            }else{
              dn     = hyp_e[i]; 
            }
	    rf     = fabs(sl_lim)/max(1e-100,fabs(sl));
	    cf     = max(0.0,min(sl_lim,sl))+min(0.0,max(sl_lim,sl));
	    flx[ivar][i] = max(0.0,1.0+dn*(rf-1.0))*cf;
	  }
         }


 #pragma acc loop vector private(tvd_faci)
	for(i=1;i<sz-2;i++){
          tvd_faci =  0.5*idxd1*min(cmax,max(cm[i],cm[i+1]));
	  flx[0][i] *= tvd_faci;
          flx[1][i] *= tvd_faci;
          flx[2][i] *= tvd_faci;
          flx[3][i] *= tvd_faci;
          flx[4][i] *= tvd_faci;
          flx[5][i] *= tvd_faci;
          flx[6][i] *= tvd_faci;
          flx[7][i] *= tvd_faci;
          flx[5+d1][i] *= B_par_diff;
       }  
      

      if(need_tvd_coeff){
#pragma acc loop vector private(cf) 
	for(i=1;i<sz-2;i++){
	  flx[0][i] *= tvd_coeff0;
	  
	  cf = max(tvd_coeff1,visc_coeff_bot*hfb1[i]+visc_coeff_top*hft1[i]);
	  flx[1][i] *= cf;
	  flx[2][i] *= cf;
	  flx[3][i] *= cf;
	  
	  flx[4][i] *= tvd_coeff2;
	  
	  cf = max(tvd_coeff3,eta_coeff_bot*hfb1[i]+eta_coeff_top*hft1[i]);
	  flx[5][i] *= cf;
	  flx[6][i] *= cf;
	  flx[7][i] *= cf; 
	}
      } 

      // additional 4th order hyperdiffusion for HD varibales in the
      // y-direction to damp spurious oscillations
      if( (d1 == 0) && (vhyp > 0.0) ){
	istart = 1;
	if(Grid.is_gbeg[0]) istart += gh;
#pragma acc loop vector private(dn, ivar)
	for(i=istart;i<sz-2;i++){
	  dn = max(hyperdiff[i],hyperdiff[i+1]);
	  ivar = 0;
	  flx[ivar][i] -= dn*(0.125*(var[ivar][i+2]-var[ivar][i-1]) -
			      0.375*(var[ivar][i+1]-var[ivar][i]));
	  ivar = 1;
	  flx[ivar][i] -= dn*(0.125*(var[ivar][i+2]-var[ivar][i-1]) -
			      0.375*(var[ivar][i+1]-var[ivar][i]));
	  ivar = 4;
	  flx[ivar][i] -= dn*(0.125*(var[ivar][i+2]-var[ivar][i-1]) -
			      0.375*(var[ivar][i+1]-var[ivar][i]));
	}
      }
#pragma acc loop vector private(flx0,flx1,flx2,flx3,flx4,uif0,uif1,uif2,uif3,uif4,min_cz)   
     for(i=1;i<sz-2;i++){ 
	 flx0 = flx[0][i];
         flx1 = flx[1][i];
         flx2 = flx[2][i];
         flx3 = flx[3][i];
         flx4 = flx[4][i];

         uif0= 0.5*(var[0][i]+slp[0][i]+
                                var[0][i+1]-slp[0][i+1]);

         uif1= 0.5*(var[1][i]+slp[1][i]+
                                var[1][i+1]-slp[1][i+1]);

         uif2= 0.5*(var[2][i]+slp[2][i]+
                                var[2][i+1]-slp[2][i+1]);

         uif3= 0.5*(var[3][i]+slp[3][i]+
                                var[3][i+1]-slp[3][i+1]);

         uif4= 0.5*(var[4][i]+slp[4][i]+
                                var[4][i+1]-slp[4][i+1]);
       if(rho_log){
           uif0  = exp(uif0);
          flx0 *= uif0;

          uif4  = exp(uif4);
          flx4 *= uif4;
         }

          min_cz = flx0*min(CZ_fac[i],CZ_fac[i+1]);
          flx1 *= uif0;
          flx1 += uif1*min_cz;

          flx2 *= uif0;
          flx2 += uif2*min_cz;

          flx3 *= uif0;
          flx3 += uif3*min_cz;

          flx4 *= uif0;
          flx4 += (uif4*min_cz)+qft1[i]*(uif1*flx1+uif2*flx2+uif3*flx3);

          flx[0][i] = flx0;
          flx[1][i] = flx1;
          flx[2][i] = flx2;
          flx[3][i] = flx3;
          flx[4][i] = flx4;

      }

      // no diffusion for vertical field at boundaries (conserve vertical magnetic flux)
      if( d1 == 0) {
	if( Grid.is_gend[0] ) flx[5][sz-gh-1] = 0.0;
	if( Grid.is_gbeg[0] ) flx[5][gh-1] = 0.0;
      } 

 #pragma acc loop vector collapse(2) 
      for(ivar=0;ivar<nvar;ivar++){
	for(i=gh;i<sz-gh;i++){
	  res[ivar][i] = dt_tvd*(flx[ivar][i]-flx[ivar][i-1]);
	}
      }

      // Remove viscous heating at top boundary
#pragma acc loop vector 
      for(i=gh;i<sz-gh;i++)
	res[4][i] += (1.0-qft1[i])*(res[1][i]*var[1][i] + res[2][i]*var[2][i] + res[3][i]*var[3][i]);

      if(need_diagnostics){
        #pragma ivdep
#pragma acc loop vector private(node) 
	for(i=gh;i<sz-gh;i++){
	  node = offset+i*str;
	  Grid.tvar8[node] += res[4][i]*idt_full;
	}
      }

      // projection of fvisc for consistency with SR treatment
#pragma acc loop vector private(dn,boris0,boris1,boris2,boris3) 
      for(i=gh;i<sz-gh;i++){
	dn = (res[1][i]*var[5][i] + res[2][i]*var[6][i] + res[3][i]*var[7][i])/max(1e-100,bsqr[i]);

        boris0 = -BC2[i]*(res[1][i]-dn*var[5][i]);
        boris1 = -BC2[i]*(res[2][i]-dn*var[6][i]);
        boris2 = -BC2[i]*(res[3][i]-dn*var[7][i]);
        boris3 =  boris0*var[1][i]+boris1*var[2][i]+boris2*var[3][i];

        res[1][i] += boris0;
        res[2][i] += boris1;
        res[3][i] += boris2;
        res[4][i] += boris3;

        if(need_diagnostics){
           node = offset+i*str;
           Grid.tvar6[node] += boris3*idt_full;
        }
      }
      
      // resistive heating
#pragma acc loop vector 
      for(i=1;i<sz-2;i++){
	qrho[i] = 
	  (var[5][i+1]-var[5][i])*flx[5][i] +
	  (var[6][i+1]-var[6][i])*flx[6][i] +
	  (var[7][i+1]-var[7][i])*flx[7][i];
      }

#pragma acc loop vector private(node,qres)      
      for(i=gh;i<sz-gh;i++){
	qres = 0.5*(qrho[i]+qrho[i-1])*qft1[i];
	res[4][i] += dt_tvd*qres;
        if(need_diagnostics){
         node = offset+i*str;
         Grid.Qres[node] += dt_fac*qres;
        }
      }
      
      // viscous heating for diagnostics
      if(need_diagnostics){
#pragma acc loop vector 
	for(i=1;i<sz-2;i++){
	  qrho[i] = 
	    (var[1][i+1]-var[1][i])*flx[1][i] +
	    (var[2][i+1]-var[2][i])*flx[2][i] +
	    (var[3][i+1]-var[3][i])*flx[3][i];
	}

      }
 
      //c_time += MPI_Wtime()-time;

      //time = MPI_Wtime();
      #pragma ivdep
#pragma acc loop vector private(node) 
      for(i=gh;i<sz-gh;i++){
	node = offset+i*str;
	U[node].d   += res[0][i];
	U[node].M.x += res[1][i];
	U[node].M.y += res[2][i];
	U[node].M.z += res[3][i];
	U[node].e   += res[4][i];
	U[node].B.x += res[5][i];
	U[node].B.y += res[6][i];
	U[node].B.z += res[7][i];
	U[node].d    = max(rho_min,U[node].d);
	
        bsqr[i] = U[node].e;
	U[node].e = max(U[node].e,eps_min*U[node].d+0.5*U[node].M.sqr()/U[node].d);
      }

      if (need_diagnostics){ 
        #pragma ivdep
#pragma acc loop vector private(node) 
	for(i=gh;i<sz-gh;i++){
	  node = offset+i*str;     
	  Grid.tvar7[node] += (U[node].e-bsqr[i])*idt_full;
          Grid.Qvis[node] += dt_fac*(0.5*(qrho[i]+qrho[i-1])*qft1[i]);
	}
      }
      //r_time += MPI_Wtime()-time;
    }

    bounds[d1][0] = Grid.lbeg[d1];
    bounds[d1][1] = Grid.lend[d1];
  }
  } //exit data 
//#pragma acc exit data delete(stride[:3],idx[:3])
  //t_time += MPI_Wtime()-s_time;
  //call_count += 1;

  // if(Run.rank == 0)
  // cout << "TVD: " << t_time/call_count << ' ' 
  // << c_time/call_count << ' ' 
  // << r_time/call_count << endl; 
  PGI_COMPARE(Grid.U, double, Grid.bufsize*8, "U", "tvdlimit_SR.C", "TVD", 1)
  if(ambipolar) {
    PGI_COMPARE(Grid.v_amb, double, Grid.bufsize*3, "v_amb", "tvdlimit_SR.C", "TVD", 2)
  }
  if(need_diagnostics) {
    PGI_COMPARE(Grid.tvar6, double, Grid.bufsize, "tvar6", "tvdlimit_SR.C", "TVD", 3)
    PGI_COMPARE(Grid.tvar7, double, Grid.bufsize, "tvar7", "tvdlimit_SR.C", "TVD", 4)
    PGI_COMPARE(Grid.tvar8, double, Grid.bufsize, "tvar8", "tvdlimit_SR.C", "TVD", 5)
    PGI_COMPARE(Grid.Qres, double, Grid.bufsize, "Qres", "tvdlimit_SR.C", "TVD", 6)
    PGI_COMPARE(Grid.Qvis, double, Grid.bufsize, "Qvis", "tvdlimit_SR.C", "TVD", 7)
  }

}

