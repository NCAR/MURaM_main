#ifndef __PHYSICS_INCLUDED__
#define __PHYSICS_INCLUDED__

#include "precision.h"

#define RT_DEFAULT    0
#define RT_NLTE       1

#ifndef WEAK_ENCAPSULATION
#define WEAK_ENCAPSULATION
#endif

#include "cfd.H"

const int i_param_grav       = 0;
const int i_param_va_max     = 1;
const int i_param_va_adjust  = 2;
const int i_param_spitzer    = 3;
const int i_param_eta        = 4;
const int i_param_max_fill   = 5;
const int i_param_ambipolar  = 6;
const int i_param_ambfac_max = 7;
const int i_param_ambvel_max = 8;

const int i_bnd_top    = 0;
const int i_bnd_pot    = 1;
const int i_bnd_bcrit  = 2;
const int i_bnd_eps_top= 3;
const int i_bnd_fem    = 4;

const int i_tvd_rholev     = 0;
const int i_tvd_rholog     = 1;
const int i_tvd_qrho       = 2;
const int i_tvd_Bpar       = 3;
const int i_tvd_vhyp       = 4;
const int i_tvd_Qdiff_bnd  = 5;
const int i_tvd_pm_v       = 6;
const int i_tvd_pm_B       = 7;
const int i_tvd_vmax_lim   = 8;
const int i_tvd_CME_thresh = 9;


const int i_divB_switch = 0;
const int i_divB_itmax  = 1;
const int i_divB_err    = 2;

const int i_tchk_eps_min  = 0;
const int i_tchk_rho_min  = 1;
const int i_tchk_eps_max  = 2;
const int i_tchk_vmax     = 3;

const int i_dmp_switch  = 0;
const int i_dmp_tau_ref = 1;
const int i_dmp_vel_ref = 2;
const int i_dmp_tau_min = 3;

const int i_rt_update  = 0;
const int i_rt_tau_min = 1;
const int i_rt_tr_tem  = 2;
const int i_rt_tstep   = 3;
const int i_rt_cfl     = 4;
const int i_rt_type    = 5;
const int i_rt_epsilon = 6;
const int i_rt_iout    = 7;
const int i_rt_tr_pre  = 8;
const int i_rt_pre_cut = 9;

const int i_ext_cor = 0;

const int i_sl_collect = 0;
const int i_sl_ic      = 1;
const int i_sl_tau     = 2;
const int i_sl_xz      = 3;
const int i_sl_xy      = 4;
const int i_sl_yz      = 5;
 
typedef MHD::cState<real> cState;
typedef Vector3D<real>    Vector;

class PhysicsData {
public:
  int NVAR;
  int layout[8];
  double params[9];
  double bnd[5];
  double tvd_h[4];
  double tvd_cs[4];
  double tvd[10];
  double tvd_h_bnd[2];
  double tvd_visc_bnd[2];
  double tvd_eta_bnd[2];  
  double tvd_coeff[4];
  double divB[3]; 
  double tchk[4];
  double dmp[4];
  double rt[10];
  int rt_ext[1];

  int    slice[6];  
  double tau_lev[20];
  int    xy_lev[20];
  int    yz_lev[20];
  int    xz_lev[20];
  int    tau_var[14];
  int    xz_var[12];
  int    xy_var[12];
  int    yz_var[13];

  Vector g;

  PhysicsData() {
    NVAR = 8;

    layout[0] = 0;
    layout[1] = 1;
    layout[2] = 2;
    layout[3] = 3;
    layout[4] = 4;
    layout[5] = 5;
    layout[6] = 6;
    layout[7] = 7;

    params[i_param_grav]       = 0.0;
    params[i_param_va_max]     = 1e7;
    params[i_param_va_adjust]  = 0;
    params[i_param_spitzer]    = 0.0;
    params[i_param_eta]        = 0.0;
    params[i_param_max_fill]   = 1000.;
    params[i_param_ambipolar]  = 0.0;
    params[i_param_ambfac_max] = 1e12;
    params[i_param_ambvel_max] = 1e7;

    bnd[i_bnd_top]    = 0.0;
    bnd[i_bnd_pot]    = 0.0;
    bnd[i_bnd_bcrit]  = 1.0e10; 
    bnd[i_bnd_eps_top]= 0.0;
    bnd[i_bnd_fem]    = 0.0;
    
    tvd_h[0]        = 2.0;
    tvd_h[1]        = 2.0;
    tvd_h[2]        = 2.0;
    tvd_h[3]        = 2.0;

    tvd_cs[0]         = 0.2;
    tvd_cs[1]         = 0.2;
    tvd_cs[2]         = 0.2;
    tvd_cs[3]         = 0.2;

    tvd_coeff[0]      = 1.0;
    tvd_coeff[1]      = 1.0;
    tvd_coeff[2]      = 1.0;
    tvd_coeff[3]      = 1.0;

    tvd[i_tvd_rholev]     = 1.0e-11;
    tvd[i_tvd_rholog]     = 1.0;
    tvd[i_tvd_qrho]       = 10.0;
    tvd[i_tvd_Bpar]       = 0.2;
    tvd[i_tvd_vhyp]       = 1.0;
    tvd[i_tvd_Qdiff_bnd]  = 1.0;
    tvd[i_tvd_pm_v]       = 1.0;
    tvd[i_tvd_pm_B]       = 1.0;
    tvd[i_tvd_vmax_lim]   = 0.75;
    tvd[i_tvd_CME_thresh] = 100.0;

    divB[i_divB_switch] = 1;
    divB[i_divB_itmax]  = 5; 
    divB[i_divB_err]    = 0.2;

    tchk[i_tchk_eps_min] = 1.0e11;
    tchk[i_tchk_rho_min] = 1.0e-16;
    tchk[i_tchk_eps_max] = 4.0e15;
    tchk[i_tchk_vmax]    = 100.0e5;  

    dmp[i_dmp_switch]  = 0;
    dmp[i_dmp_tau_ref] = 1.0e3;
    dmp[i_dmp_vel_ref] = 1.0e3;
    dmp[i_dmp_tau_min] = 1.0e2;

    rt[i_rt_update]  = 1;
    rt[i_rt_tau_min] = 1.0e-8;
    rt[i_rt_tr_tem]  = 20.0e3;
    rt[i_rt_tr_pre]  = 1e2; 
    rt[i_rt_cfl]     = 0.5;
    rt[i_rt_tstep]   = 0.0;
    rt[i_rt_type]    = 0;
    rt[i_rt_epsilon] = 0.0;  
    rt[i_rt_iout]    = 0;
    rt[i_rt_pre_cut] = 1e2;

    rt_ext[i_ext_cor]       = 0;

    slice[i_sl_collect] = 0;
    slice[i_sl_ic]      = 0;
    slice[i_sl_tau]     = 0;
    slice[i_sl_xz]      = 0;
    slice[i_sl_xy]      = 0;
    slice[i_sl_yz]      = 0;    
  }
  ~PhysicsData() {}

  void Init() {
    g = Vector(-params[i_param_grav],0.,0.);
  }

  void Show() const;
};

#endif
