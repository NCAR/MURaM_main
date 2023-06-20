#include "physics.H"
#include <iostream>

using namespace std;

void PhysicsData::Show() const {

  cout << " ------------ Physics Parameter Settings -------------" << endl;
  cout << "Parameters: " << endl
       << "  param_gravity    = " << params[i_param_grav] << endl
       << "  param_va_max     = " << params[i_param_va_max] << endl
       << "  param_va_adjust  = " << params[i_param_va_adjust]<< endl
       << "  param_max_fill   = " << params[i_param_max_fill]<< endl
       << "  param_spitzer    = " << params[i_param_spitzer]<< endl
       << "  param_eta        = " << params[i_param_eta]<< endl
       << "  param_ambipolar  = " << params[i_param_ambipolar]<< endl
       << "  param_ambfac_max = " << params[i_param_ambfac_max]<< endl
       << "  param_ambvel_max = " << params[i_param_ambvel_max]<< endl;

  cout << "Boundary condition: " << endl
       << "  bnd_top     = " << bnd[i_bnd_top]   << endl
       << "  bnd_pot     = " << bnd[i_bnd_pot]   << endl
       << "  bnd_bcrit   = " << bnd[i_bnd_bcrit] << endl
       << "  bnd_eps_top = " << bnd[i_bnd_eps_top]  << endl
       << "  bnd_fem     = " << bnd[i_bnd_fem]  << endl;

  cout << "Tvd limiter: "  << endl
       << "  tvd_h          = "   << tvd_h[0]          << ' ' <<  tvd_h[1]
       << " "              << tvd_h[2]          << ' ' <<  tvd_h[3]   << endl
       << "  tvd_cs         = " << tvd_cs[0]         << ' ' <<  tvd_cs[1]
       << " "              << tvd_cs[2]         << ' ' <<  tvd_cs[3]  << endl
       << "  tvd_rholev     = " << tvd[i_tvd_rholev] << endl
       << "  tvd_rholog     = " << tvd[i_tvd_rholog] << endl
       << "  tvd_qrho       = " << tvd[i_tvd_qrho]   << endl
       << "  tvd_vhyp       = " << tvd[i_tvd_vhyp]   << endl
       << "  tvd_Qdiff_bnd  = " << tvd[i_tvd_Qdiff_bnd]   << endl
       << "  tvd_pm_v       = " << tvd[i_tvd_pm_v]   << endl
       << "  tvd_pm_B       = " << tvd[i_tvd_pm_B]   << endl
       << "  tvd_vmax_lim   = " << tvd[i_tvd_vmax_lim]   << endl
       << "  tvd_CME_thresh = " << tvd[i_tvd_CME_thresh]   << endl;
  
  cout << "divB cleaning: "<< endl
       << "  divB_switch = " << divB[i_divB_switch] << endl
       << "  divB_itmax  = " << divB[i_divB_itmax]  << endl
       << "  divB_err    = " << divB[i_divB_err]    << endl;

  cout << "Tcheck: " << endl
       << "  tchk_rho_min  = " << tchk[i_tchk_rho_min]  << endl
       << "  tchk_eps_min  = " << tchk[i_tchk_eps_min]  << endl
       << "  tchk_eps_max  = " << tchk[i_tchk_eps_max]  << endl
       << "  tchk_vmax     = " << tchk[i_tchk_vmax]     << endl; 

  cout << "Damp box mode: " << endl
       << "  dmp_switch  = " << dmp[i_dmp_switch]  << endl
       << "  dmp_tau_ref = " << dmp[i_dmp_tau_ref] << endl
       << "  dmp_vel_ref = " << dmp[i_dmp_vel_ref] << endl
       << "  dmp_tau_min = " << dmp[i_dmp_tau_min] << endl;

  cout << "RT: " << endl
       << "  rt_tau_min = " << rt[i_rt_tau_min] << endl
       << "  rt_tr_tem  = " << rt[i_rt_tr_tem] << endl
       << "  rt_tr_pre  = " << rt[i_rt_tr_pre] << endl
       << "  rt_pre_cut = " << rt[i_rt_pre_cut] << endl
       << "  rt_update  = " << rt[i_rt_update] << endl
       << "  rt_tstep   = " << rt[i_rt_tstep] << endl
       << "  rt_cfl     = " << rt[i_rt_cfl] << endl
       << "  rt_type    = " << rt[i_rt_type] << endl
       << "  rt_epsilon = " << rt[i_rt_epsilon] << endl
       << "  rt_iout    = " << rt[i_rt_iout] << endl;

  cout << "RT Extended: " << endl
       << "  rt_ext_cor     = " << rt_ext[i_ext_cor] << endl;

  cout << "Slice: " << endl
       << "  sl_collect = " << slice[i_sl_collect] << endl 
       << "  sl_I_out   = " << slice[i_sl_ic] << endl
       << "  sl_tau     = " << slice[i_sl_tau] << endl
       << "  sl_xz      = " << slice[i_sl_xz] << endl
       << "  sl_xy      = " << slice[i_sl_xy] << endl
       << "  sl_yz      = " << slice[i_sl_yz] << endl;
}
