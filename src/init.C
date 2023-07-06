#include <mpi.h>
#include <iostream>
#include <string.h>
#include "init.H"
#include "precision.h"
#include "dfparser.h"
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "io.H"
#include "comm_split.H"
#include "rt/rt.h"
#include "ACCH.h"

using namespace std;

int Initialize(RunData& Run,GridData& Grid,
	       PhysicsData& Physics,RTS *&rts) {
  char datafile[256] = "parameters.dat";
  char rtype[16] = "double";
  int mode = newrun;
  int i,rank;
  if( sizeof(real)==sizeof(float) ) strcpy(rtype,"float");

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  ACCH::SetGPU(rank%ACCH::GetNumGPU());

  if( rank==0 ) {
    int flag = 0;
    char ratype[16];
    strcpy(ratype,rtype);

    getvar(Run.anlfile,"anlfile","char*",datafile);
    getvar(Run.resfile,"resfile","char*",datafile);
    getvar(Run.backfile,"backfile","char*",datafile);
    getvar(&Run.maxiter,"maxiter","int",datafile);
    getvar(&Run.anlfreq,"anlfreq","int",datafile);
    getvar(&Run.resfreq,"resfreq","int",datafile);
    getvar(&Run.backfreq,"backfreq","int",datafile);
    getvar(&Run.outcad,"outcad",rtype,datafile);
    getvar(&Run.slicefreq,"slicefreq","int",datafile);
    getvar(&Run.dt,"dt",rtype,datafile);
    getvar(&Run.Tmax,"Tmax",rtype,datafile);
    getvar(&Run.CFL,"CFL",rtype,datafile);
    getvar(&Run.CFL_tvd,"CFL_tvd",rtype,datafile);
    getvar(&Run.maxWtime,"maxWtime",rtype,datafile);
    getvar(&Run.comment,"comment","char*",datafile);
    getvar(&Run.verbose,"verbose","int",datafile);
    getvar(&Run.path_3D,"path_3D","char*",datafile);
    getvar(&Run.path_2D,"path_2D","char*",datafile);
    getvar(&Run.eos_name,"eos_name","char*",datafile);
    getvar(&Run.kap_name,"kap_name","char*",datafile);
    getvar(&Run.eos_output,"eos_output","int[14]",datafile);
    getvar(&Run.diagnostics,"diagnostics","int",datafile);
    getvar(&Run.HAVG,"HAVG","int",datafile);
    getvar(&Run.RT_HAVG,"RT_HAVG","int",datafile);
    getvar(&Run.DEM,"DEM","int",datafile);
    getvar(&Run.diag_output,"diag_output","int[11]",datafile);

    getvar_s(&Grid.NDIM,"NDIM","int",datafile);

    if( Grid.NDIM==3 ) {
      strcat(ratype,"[3]");
      getvar(Grid.gxmin,"gxmin",ratype,datafile);
      getvar(Grid.gxmax,"gxmax",ratype,datafile);
      getvar(Grid.gsize,"gsize","int[3]",datafile);
      getvar(Grid.pardim,"pardim","int[3]",datafile);
      getvar(Grid.periods,"periods","int[3]",datafile);
      getvar(Grid.procs,"procs","int[3]",datafile);
      flag = getvar(Grid.ghosts,"ghosts","int[3]",datafile);
    } else if (Grid.NDIM==2) {
      strcat(ratype,"[2]");
      getvar(Grid.gxmin,"gxmin",ratype,datafile);
      getvar(Grid.gxmax,"gxmax",ratype,datafile);
      getvar(Grid.gsize,"gsize","int[2]",datafile);
      getvar(Grid.pardim,"pardim","int[2]",datafile);
      getvar(Grid.periods,"periods","int[2]",datafile);
      getvar(Grid.procs,"procs","int[2]",datafile);
      flag = getvar(Grid.ghosts,"ghosts","int[2]",datafile);
    } else if (Grid.NDIM==1) {
      strcat(ratype,"[1]");
      getvar(Grid.gxmin,"gxmin",ratype,datafile);
      getvar(Grid.gxmax,"gxmax",ratype,datafile);
      getvar(Grid.gsize,"gsize","int[1]",datafile);
      getvar(Grid.pardim,"pardim","int[1]",datafile);
      getvar(Grid.periods,"periods","int[1]",datafile);
      getvar(Grid.procs,"procs","int[1]",datafile);
      flag = getvar(Grid.ghosts,"ghosts","int[1]",datafile);
    }

    if( flag )
      for(i=0;i<Grid.NDIM;i++) Grid.ghosts[i] = 2;
    
    getvar_s(&Physics.params[i_param_grav], "param_gravity",rtype,datafile);
    getvar_s(&Physics.params[i_param_va_max],"param_va_max",rtype,datafile);
    getvar_s(&Physics.params[i_param_va_adjust],"param_va_adjust",rtype,datafile);
    getvar_s(&Physics.params[i_param_spitzer],"param_spitzer",rtype,datafile);
    getvar(&Physics.params[i_param_eta],"param_eta",rtype,datafile);
    getvar(&Physics.params[i_param_max_fill],"param_max_fill",rtype,datafile);
    getvar(&Physics.params[i_param_ambipolar],"param_ambipolar",rtype,datafile);
    getvar(&Physics.params[i_param_ambfac_max],"param_ambfac_max",rtype,datafile);
    getvar(&Physics.params[i_param_ambvel_max],"param_ambvel_max",rtype,datafile);

    if (Physics.params[i_param_ambipolar] > 0)
    {
      cout << "The Current implementation of ambipolar diffusion is out of date and likely incorrect" << endl;
      cout << "The code must be updated to the version presented in Rempel + Przybylski 2021" << endl;
      cout << "Set param_ambipolar = 0" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    getvar_s(&Physics.bnd[i_bnd_top],   "bnd_top",rtype,datafile);
    getvar_s(&Physics.bnd[i_bnd_pot],   "bnd_pot",rtype,datafile);
    getvar_s(&Physics.bnd[i_bnd_bcrit], "bnd_bcrit",rtype,datafile);
    getvar(&Physics.bnd[i_bnd_eps_top], "bnd_eps_top",rtype,datafile);
    getvar(&Physics.bnd[i_bnd_fem],     "bnd_fem",rtype,datafile);

    getvar_s(&Physics.tvd_h,           "tvd_h","double[4]",datafile);
    getvar(&Physics.tvd_cs,            "tvd_cs","double[4]",datafile);
    getvar(&Physics.tvd[i_tvd_rholev], "tvd_rholev",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_rholog], "tvd_rholog",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_qrho],   "tvd_qrho",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_Bpar],   "tvd_Bpar",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_vhyp],   "tvd_vhyp",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_Qdiff_bnd], "tvd_Qdiff_bnd",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_pm_v],      "tvd_pm_v",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_pm_B],      "tvd_pm_B",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_vmax_lim],  "tvd_vmax_lim",rtype,datafile);
    getvar(&Physics.tvd[i_tvd_CME_thresh],"tvd_CME_thresh",rtype,datafile);

    getvar_s(&Physics.tvd_h_bnd,   "tvd_h_bnd",   "double[2]",datafile);
    getvar_s(&Physics.tvd_visc_bnd,"tvd_visc_bnd","double[2]",datafile);
    getvar_s(&Physics.tvd_eta_bnd, "tvd_eta_bnd", "double[2]",datafile);
    getvar(&Physics.tvd_coeff,   "tvd_coeff",   "double[4]",datafile);

    getvar_s(&Physics.divB[i_divB_switch],"divB_switch",rtype,datafile);
    getvar(&Physics.divB[i_divB_itmax], "divB_itmax",rtype,datafile);
    getvar(&Physics.divB[i_divB_err],   "divB_err",rtype,datafile);

    getvar_s(&Physics.tchk[i_tchk_eps_min], "tchk_eps_min",rtype,datafile);
    getvar_s(&Physics.tchk[i_tchk_rho_min], "tchk_rho_min",rtype,datafile);
    getvar_s(&Physics.tchk[i_tchk_eps_max], "tchk_eps_max",rtype,datafile);
    getvar_s(&Physics.tchk[i_tchk_vmax],    "tchk_vmax",rtype,datafile);

    getvar_s(&Physics.dmp[i_dmp_switch], "dmp_switch",rtype,datafile);
    getvar_s(&Physics.dmp[i_dmp_tau_ref],"dmp_tau_ref",rtype,datafile);
    getvar_s(&Physics.dmp[i_dmp_vel_ref],"dmp_vel_ref",rtype,datafile);
    getvar_s(&Physics.dmp[i_dmp_tau_min],"dmp_tau_min",rtype,datafile);

    getvar_s(&Physics.rt[i_rt_update],"rt_update",rtype,datafile);
    getvar(&Physics.rt[i_rt_tau_min],"rt_tau_min",rtype,datafile);
    getvar(&Physics.rt[i_rt_tr_tem],"rt_tr_tem",rtype,datafile);
    getvar(&Physics.rt[i_rt_tr_pre],"rt_tr_pre",rtype,datafile);
    getvar(&Physics.rt[i_rt_pre_cut],"rt_pre_cut",rtype,datafile);
    getvar_s(&Physics.rt[i_rt_tstep], "rt_tstep",rtype,datafile);
    getvar(&Physics.rt[i_rt_cfl],"rt_cfl",rtype,datafile);
    getvar_s(&Physics.rt[i_rt_type],"rt_type",rtype,datafile);
    getvar(&Physics.rt[i_rt_epsilon],"rt_epsilon",rtype,datafile);
    getvar(&Physics.rt[i_rt_iout],"rt_iout",rtype,datafile);
    
    getvar_s(&Physics.rt_ext[i_ext_cor],"ext_cor","int",datafile);

    // If set to off (0) or Chianti (1) then turn line losses off

    getvar(&Physics.slice[i_sl_collect],"sl_collect","int",datafile);
    getvar(&Physics.slice[i_sl_ic],     "sl_I_out","int",datafile); 

    int nsl;
    getvar(&nsl,"sl_tau","int",datafile); 
    if (nsl > 20) nsl = 20;
    Physics.slice[i_sl_tau] = nsl;
    getvar(&Physics.tau_lev,"tau_lev","double[20]",datafile);

    getvar(&nsl,"sl_xz","int",datafile); 
    if (nsl > 20) nsl = 20;

    Physics.slice[i_sl_xz] = nsl;    
    getvar(&Physics.xz_lev,"xz_lev","int[20]",datafile);

    getvar(&nsl,"sl_xy","int",datafile); 
    if (nsl > 20) nsl = 20;
    Physics.slice[i_sl_xy] = nsl;    
    getvar(&Physics.xy_lev,"xy_lev","int[20]",datafile);

    getvar(&nsl,"sl_yz","int",datafile); 
    if (nsl > 20) nsl = 20;
    Physics.slice[i_sl_yz] = nsl;    
    getvar(&Physics.yz_lev,"yz_lev","int[20]",datafile);

    getvar(&Physics.tau_var,"tau_var","int[14]",datafile); 
    getvar(&Physics.xz_var,"xz_var","int[12]",datafile); 
    getvar(&Physics.xy_var,"xy_var","int[12]",datafile); 
    getvar(&Physics.yz_var,"yz_var","int[13]",datafile);    
  }

  MPI_Bcast(&Run,sizeof(Run),MPI_BYTE,0,MPI_COMM_WORLD);

  MPI_Bcast(&Grid.NDIM,  1,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.gxmin,  3,REALTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.gxmax,  3,REALTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.gsize,  3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.ghosts, 3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.pardim, 3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.periods,3,MPI_INT, 0,MPI_COMM_WORLD);
  MPI_Bcast(Grid.procs,3,MPI_INT, 0,MPI_COMM_WORLD);

  MPI_Bcast(&Physics,sizeof(Physics),MPI_BYTE,0,MPI_COMM_WORLD);
  
  int tot_procs=1;
  for (int dim = 0; dim<3;dim++)
    tot_procs*=Grid.procs[dim];

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (tot_procs != nprocs) {

    cout << "!!! ERROR procs set in parameters.dat != number of processors!!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);

  }

  Run.Init(rank);
  Physics.Init();
  Grid.Init(Run,Physics);
  comm_split_init(Run,Grid);

  if( rank==0 ) {
    if( HasBackupFile(Run.backfile) )
      mode = restart;
    else
      mode = newrun;
    
    Run.Show();
    Grid.Show();
    Physics.Show();
  }
  
  rts=rt_new(Grid,Run,Physics);

  if (rts == 0) {

    cout << "!!! ERROR initialising RTS, check rt.cfg !!!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);

  }
    
  MPI_Bcast(&mode,1,MPI_INT,0,MPI_COMM_WORLD);

  IO_Init(Grid);

  return mode;
}
