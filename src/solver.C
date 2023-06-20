#include "physics.H"
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <limits>
#include "solver.H"
#include "divB.H"
#include "exchange.H"
#include "mhd_tvd.H"
#include "limit_va.H"
#include "src_int_tck.H"
#include "rt.H"
#include "analysis.H"
#include "grid.H"
#include "run.H"
#include "io.H"
#include "eos.H"
#include "icbc.H"
#include "rt/rt.h"
#include "ACCH.h"

using std::cout;
using std::endl;
using std::min;
using std::max;

extern double start_time,rst_time;

void Sync_timestep(RunData&, const PhysicsData&, const double, const double, 
		   int*, double*);

extern void Get_Radloss(const RunData&  Run, GridData& Grid,const PhysicsData& Physics);

void ComputeSolution(RunData& Run,GridData& Grid,const PhysicsData& Physics,RTS *rts) {

  int stage,rt_upd;
  double dt_rad,dt_mhd,cfl;
  double rt_time,tvd_time,divB_time,eos_time,bnd_time,
    mhd_time,clock,ini_time,io_time,bc_time,
    current_time,dst_time,int_time,vlim_time,sync_time;
  int wtime_backup;
  double inv_it, inv_time;

  int pt_update;

  double maxWtime,endWtime;

  maxWtime=fabs(Run.maxWtime);

  if(Run.rank == 0) cout << "eos_init ..." << endl;
  eos_init(Grid,Run,Physics);
  if( Run.rank==0 ) cout << "eos_init done"<<endl;

  ACCH::UpdateGPU(Grid.U, Grid.bufsize*sizeof(cState));
  ACCH::UpdateGPU(Grid.U0, Grid.bufsize*sizeof(cState));
  ACCH::UpdateGPU(Grid.Res, Grid.bufsize*sizeof(cState));
  if(Physics.params[i_param_spitzer] > 0.0) {
    ACCH::UpdateGPU(Grid.sflx, Grid.bufsize*sizeof(double));
    ACCH::UpdateGPU(Grid.sflx0, Grid.bufsize*sizeof(double));
    ACCH::UpdateGPU(Grid.Rflx, Grid.bufsize*sizeof(double));
  }
  if(Physics.params[i_param_ambipolar] > 0.0) {
    ACCH::UpdateGPU(Grid.v_amb, Grid.bufsize*sizeof(Vector));
    ACCH::UpdateGPU(Grid.v0_amb, Grid.bufsize*sizeof(Vector));
  }
  ACCH::UpdateGPU(Grid.divB, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.pres, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.Jtot, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.Qtot, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.Stot, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.Qthin, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.rhoi, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.Tau, Grid.bufsize*sizeof(double));
  ACCH::UpdateGPU(Grid.R_amb, Grid.bufsize*sizeof(Vector));

  //TCheck(Run,Grid,Physics);  
  exchange_grid_acc(Grid,Physics,1);
  SetBoundaryConditions(Run,Grid,Physics,0,1,rts);
  MPI_Barrier(MPI_COMM_WORLD);
  bnd_time+=MPI_Wtime()-clock;
  
  ini_time   = MPI_Wtime()-start_time;
  start_time = MPI_Wtime(); 

  eos_time = 0.0;
  bnd_time = 0.0;
  dst_time = 0.0;
  int_time = 0.0;
  mhd_time = 0.0;
  rt_time  = 0.0;
  tvd_time = 0.0;
  divB_time= 0.0;
  io_time  = 0.0;
  vlim_time= 0.0;
  sync_time= 0.0;

  bc_time  = rst_time;
  
  endWtime = 0;

  dt_rad = 0.0;
  cfl    = 0.0;

  rt_upd = 1;

  while( Run.IsValid() and (endWtime == 0) ) {

    if (Run.rank ==0)
      cout << "***** Start new iteration ***********************" << endl;
    
    clock=MPI_Wtime();
    SaveCons(Grid,Physics);
    int_time+=MPI_Wtime()-clock;

    for(stage=1;stage<=maxstage;stage++) {

      clock=MPI_Wtime();
      ConsToPrim(Grid,Physics,Run);
      eos_time+=MPI_Wtime()-clock;

      if(stage == 1){
	clock=MPI_Wtime();
        Adjust_Valf_Max(Run,Grid,Physics);
	vlim_time+=MPI_Wtime()-clock;
      }

      clock=MPI_Wtime();
      if( (stage == 1) or (stage == maxstage) ) {
	if(Physics.rt_ext[i_ext_cor]==1) {
          Get_Radloss(Run,Grid,Physics);
        }
      }
      int_time+=MPI_Wtime()-clock;

      if ( stage == 1 ){ 

	clock=MPI_Wtime();
	dt_rad=rts->wrapper(rt_upd,Grid,Run,Physics);
        rt_time+=MPI_Wtime()-clock;
    
        int needoutput = ((Run.NeedsOutput())&&(Run.outcad==0)||(Run.dt_rem()==0.0)&&(Run.outcad>0));
      
        if(needoutput){
	  clock=MPI_Wtime();
          ACCH::UpdateCPU(Grid.Jtot, Grid.bufsize*sizeof(double));
          ACCH::UpdateCPU(Grid.Qtot, Grid.bufsize*sizeof(double));
          ACCH::UpdateCPU(Grid.Stot, Grid.bufsize*sizeof(double));
          ACCH::UpdateCPU(Grid.ne, Grid.bufsize*sizeof(double));
          ACCH::UpdateCPU(Grid.Qthin, Grid.bufsize*sizeof(double));
          ACCH::UpdateCPU(Grid.rhoi, Grid.bufsize*sizeof(double));
          ACCH::UpdateCPU(Grid.amb, Grid.bufsize*sizeof(double));
          ACCH::UpdateCPU(Grid.Tau, Grid.bufsize*sizeof(double));
	  ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
	  eos_output(Run,Grid,Physics,rts);

	  if (Run.diagnostics){
            ACCH::UpdateCPU(Grid.tvar1, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar2, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar3, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar4, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar5, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar6, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar7, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar8, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.Qres, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.Qvis, Grid.bufsize*sizeof(double));
            if(Physics.params[i_param_ambipolar] > 0.0) {
              ACCH::UpdateCPU(Grid.Qamb, Grid.bufsize*sizeof(double));
            }
            diag_output(Run,Grid,Physics,rts);
          }
          
	  io_time+=MPI_Wtime()-clock;
	  if (Run.rank == 0){
	    if(Run.verbose >0) cout << "Output (EOS) in " << MPI_Wtime()-clock 
		 << " seconds" << endl;
	  }
	}
	
	/* output slices at tau=1 and in yz-plane */ 
	if ( Run.NeedsSlice() ){ 
	  clock=MPI_Wtime(); 
	  if (Physics.slice[i_sl_ic] >0) Iout(Run,Grid,Physics,rts);
	  if (Grid.NDIM > 1){
	    if (Physics.slice[i_sl_tau]>0) {
              ACCH::UpdateCPU(Grid.pres, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.Tau, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
              tau_slice(Run,Grid,Physics,rts);
            } 
	    if (Physics.slice[i_sl_yz]>0) {
              ACCH::UpdateCPU(Grid.pres, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
              yz_slice(Run,Grid,Physics,rts);
            }
	  }
	  if (Grid.NDIM == 3){
	    if (Physics.slice[i_sl_xy]>0) {
              ACCH::UpdateCPU(Grid.pres, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
              xy_slice(Run,Grid,Physics);
            }
	    if (Physics.slice[i_sl_xz]>0) {
              ACCH::UpdateCPU(Grid.pres, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
              ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
              xz_slice(Run,Grid,Physics);
            }
	  }

	  if (Run.DEM) {
            //commented as part of IO computation porting
            //ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
            //ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
	    corona_emission_dem_xyz(Run,Grid,Physics);
          }

	  // Some quantities might not be yeat computed after restart
	  if((Run.iteration > 0)&&(Run.HAVG)) {
           //commented as part of IO computation porting
           /* ACCH::UpdateCPU(Grid.tvar1, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar2, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar3, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar4, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar5, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar6, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar7, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.tvar8, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.divB, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.Qres, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.Qvis, Grid.bufsize*sizeof(double));
            if(Physics.params[i_param_ambipolar] > 0.0) {
              ACCH::UpdateCPU(Grid.Qamb, Grid.bufsize*sizeof(double));
            }
            ACCH::UpdateCPU(Grid.pres, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.Qtot, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.Qthin, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.rhoi, Grid.bufsize*sizeof(double));
            ACCH::UpdateCPU(Grid.Tau, Grid.bufsize*sizeof(double));
	    ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
           */
	    AnalyzeSolution_VP(Run,Grid,Physics,rts);
	  }
         
	  io_time += MPI_Wtime()-clock;    
	  if (Run.rank == 0){
	    if(Run.verbose >0) cout << "slice output in " << MPI_Wtime()-clock 
		 << " seconds" << endl;
	  }
	}
      }
      //END STAGE=1 SPECIAL
      
      if ((stage == maxstage)&&(Run.diagnostics==1))
        Run.need_diagnostics = !((Run.globiter+1)%Run.slicefreq) or !((Run.globiter+1)%Run.resfreq);
      else
        Run.need_diagnostics = 0;
     
      int buf = Grid.bufsize; 

      if(Run.need_diagnostics){
#pragma acc parallel loop \
 present(Grid[:1], Grid.tvar1[:Grid.bufsize], \
  Grid.tvar2[:Grid.bufsize], Grid.tvar3[:Grid.bufsize], \
  Grid.tvar4[:Grid.bufsize], Grid.tvar5[:Grid.bufsize], \
  Grid.tvar6[:Grid.bufsize], Grid.tvar7[:Grid.bufsize], \
  Grid.tvar8[:Grid.bufsize], Grid.Qres[:Grid.bufsize], \
  Grid.Qvis[:Grid.bufsize])
	for(int i = 0; i < Grid.bufsize; i++) {
          Grid.tvar1[i] = 0.0;
	  Grid.tvar2[i] = 0.0;
	  Grid.tvar3[i] = 0.0;
	  Grid.tvar4[i] = 0.0;
	  Grid.tvar5[i] = 0.0;
	  Grid.tvar6[i] = 0.0;
	  Grid.tvar7[i] = 0.0;
	  Grid.tvar8[i] = 0.0;
	  Grid.Qres[i] = 0.0;
	  Grid.Qvis[i] = 0.0;
	}

	if(Physics.params[i_param_ambipolar] > 0.0) {
#pragma acc parallel loop \
 present(Grid[:1], Grid.Qamb[:Grid.bufsize])
	  for(int i = 0; i < Grid.bufsize; i++) {
            Grid.Qamb[i] = 0.0;
	  }
	}
      }

      clock=MPI_Wtime();
      dt_mhd = MHD_Residual(Run,Grid,Physics);       
      mhd_time+=MPI_Wtime()-clock;

      if( stage==1 ){
	clock=MPI_Wtime();
	Sync_timestep(Run,Physics,dt_rad,dt_mhd,&rt_upd,&cfl);
	sync_time+=MPI_Wtime()-clock;
      }

      clock=MPI_Wtime();
      Source_Integrate_Tcheck(Run,Grid,Physics,stage);
      int_time+=MPI_Wtime()-clock;

      clock=MPI_Wtime();
      exchange_grid_acc(Grid,Physics,1);   
      dst_time+=MPI_Wtime()-clock;

      clock=MPI_Wtime();
      SetBoundaryConditions(Run,Grid,Physics,stage,0,rts);
      bnd_time+=MPI_Wtime()-clock;

    }

    if (cfl <= Run.CFL_tvd){
      clock=MPI_Wtime();
      TVDlimit(Run,Grid,Physics,Run.dt);
      tvd_time+=MPI_Wtime()-clock;
    } else {
      clock=MPI_Wtime();
      TVDlimit(Run,Grid,Physics,0.5*Run.dt);
      tvd_time+=MPI_Wtime()-clock;
      clock=MPI_Wtime();
      exchange_grid_acc(Grid,Physics,0);
      dst_time+=MPI_Wtime()-clock;
      clock=MPI_Wtime();
      TVDlimit(Run,Grid,Physics,0.5*Run.dt);
      tvd_time+=MPI_Wtime()-clock;
    }

    if (Physics.divB[i_divB_switch] > 0){
      clock=MPI_Wtime();
      Clean_divB(Run,Grid,Physics);
      divB_time+=MPI_Wtime()-clock;
    }

    clock=MPI_Wtime();
    TCheck(Run,Grid,Physics);
    int_time+=MPI_Wtime()-clock;
    
    clock=MPI_Wtime();
    exchange_grid_acc(Grid,Physics,0); 
    dst_time+=MPI_Wtime()-clock;

    clock=MPI_Wtime();
    pt_update=0;
    if ((int) Physics.bnd[i_bnd_pot] <= 2)
      pt_update=1;
    else 
      if ( !(Run.iteration%( ( (int) Physics.bnd[i_bnd_pot] ) - 1 ) ) )
	pt_update=1;

    SetBoundaryConditions(Run,Grid,Physics,maxstage+1,pt_update,rts);
    bnd_time+=MPI_Wtime()-clock;
    Run.Advance();

    if (Run.rank == 0){
      cout << "****** Finish iteration *************************" << endl;
      cout << "Time [" << Run.globiter << "] = " << Run.time << " seconds" << " maxiter= " << Run.maxiter << endl;
    }

    /* output of some global quantities */
    if ( !(Run.globiter%Run.anlfreq)){ 
      clock=MPI_Wtime();    
	  ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
      ACCH::UpdateCPU(Grid.Tau, Grid.bufsize*sizeof(double));
      AnalyzeSolution(Run,Grid,Physics,rts);
      io_time += MPI_Wtime()-clock;         
    }

    /* 
       Root determines if end of run backup needed (to avoid problems in the
       case that MPI_Wtime is not sufficiently syncronized)    
    */
    wtime_backup = 0;
    if( Run.rank == 0 ){
      current_time  = MPI_Wtime()-start_time+ini_time;
      if (current_time >  maxWtime-1.5*bc_time){
	wtime_backup = 1;
      }
    }

    MPI_Bcast(&wtime_backup,1,MPI_INT,0,MPI_COMM_WORLD); 

    if ((wtime_backup == 1) && (Run.maxWtime > 0) ){
      maxWtime = 1.e7;  
    }

    if ((wtime_backup == 1) && (Run.maxWtime < 0) ){
      endWtime = 1;  
    }   

    int needoutput = (((Run.NeedsOutput())&&(Run.outcad==0)) || ((Run.dt_rem()==0)&&(Run.outcad>0)) || 
		       Run.NeedsBackup() || (wtime_backup == 1));

    if(needoutput){
      clock=MPI_Wtime(); 
      ACCH::UpdateCPU(Grid.U, Grid.bufsize*sizeof(cState));
      ACCH::UpdateCPU(Grid.U0, Grid.bufsize*sizeof(cState));
      if(Physics.params[i_param_spitzer] > 0.0) 
      {
        ACCH::UpdateCPU(Grid.sflx, Grid.bufsize*sizeof(double));
        ACCH::UpdateCPU(Grid.sflx0, Grid.bufsize*sizeof(double));
      }
      ACCH::UpdateCPU(Grid.temp, Grid.bufsize*sizeof(double));
      ACCH::UpdateCPU(Grid.pres, Grid.bufsize*sizeof(double));
      BackupSolution(Run,Grid,Physics);
      ACCH::UpdateGPU(Grid.U0, Grid.bufsize*sizeof(cState));
      bc_time = MPI_Wtime()-clock;
      
      if ((Run.rank == 0)&&(Run.verbose >0))
          cout << "Output(Res)in " << MPI_Wtime()-clock << " seconds" << endl;
      
      io_time += MPI_Wtime()-clock;

    }
 
    if( Run.rank==0 ){
      inv_it     = 1.0/((double)Run.iteration);

      if(Run.verbose>1){
	cout << "Timing Info:" << endl;
	cout << "  Wallclock Time [" << Run.iteration <<  "] = " 
	     << ini_time+MPI_Wtime()-start_time << "    maxWtime: " 
	     << maxWtime << "    Final IO: " 
	     << maxWtime-1.5*bc_time << endl;
	cout << "  Initialization time           = " << ini_time << endl;
	cout << "  Total IO time                 = " << io_time << endl;
	cout << "  Average time step (incl IO)   = " 
	     << (MPI_Wtime()-start_time)*inv_it << endl;
	cout << "  Average time step (excl IO)   = " 
	     << (MPI_Wtime()-start_time-io_time)*inv_it << endl;

	inv_time = 1.0/(MPI_Wtime()-start_time-io_time);

	cout << "Breakdown by routine (without IO):" << endl;
	cout << "  rt_time:   "<< rt_time*inv_it     <<" sec "
	     << 100*rt_time*inv_time << " % "  << endl; 
	cout << "  mhd_time:  "<< mhd_time*inv_it    <<" sec "
	     << 100*mhd_time*inv_time << " % " << endl; 
	cout << "  tvd_time:  "<< tvd_time*inv_it    <<" sec "
	   << 100*tvd_time*inv_time << " % " << endl;
	cout << "  eos_time:  "<< eos_time*inv_it    <<" sec "
	   << 100*eos_time*inv_time << " % " << endl;
	cout << "  bnd_time:  "<< bnd_time*inv_it    <<" sec "
	   << 100*bnd_time*inv_time << " % " << endl;
	cout << "  int_time:  "<< int_time*inv_it    <<" sec "
	   << 100*int_time*inv_time << " % " << endl;  
	cout << "  dst_time:  "<< dst_time*inv_it    <<" sec "
	   << 100*dst_time*inv_time << " % " << endl; 
	cout << "  divB_time: "<< divB_time*inv_it   <<" sec "
	   << 100*divB_time*inv_time << " % " << endl;
	cout << "  vlim_time: "<< vlim_time*inv_it   <<" sec "
           << 100*vlim_time*inv_time << " % " << endl;
        cout << "  sync_time: "<< sync_time*inv_it   <<" sec "
           << 100*sync_time*inv_time << " % " << endl;

	double total_time = rt_time + mhd_time + tvd_time + eos_time + bnd_time + int_time + dst_time + divB_time + vlim_time + sync_time;  
	cout << "  Total:     "<< total_time*inv_it    <<" sec "
           << 100*total_time*inv_time << " % " << endl;
      }
    }   
  }
}
/*******************************************************************************/
inline int invalid_dt(double a){
  int status = 0;

  if( (a != a) or (a <= 0.) or (a == std::numeric_limits<double>::infinity()) )
    status = 1;

  return status;
}
inline int imin(int a, int b) { return a < b ? a : b; }
inline int imax(int a, int b) { return a > b ? a : b; }
/*******************************************************************************/
void Sync_timestep(RunData& Run, const PhysicsData& Physics, const double dt_rad, 
		   const double dt_mhd, int *rt_upd, double *cfl){

  double sbuf[2],rbuf[2];
  double dt,dt_min,dt_min_rt;

  sbuf[0] = dt_rad;
  sbuf[1] = dt_mhd;
  MPI_Allreduce(sbuf,rbuf,2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  dt  = rbuf[0];
  if( invalid_dt(dt) ){
    if(Run.rank ==0)
      cout << "Invalid time step: dt_rad = " << dt << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  dt  = rbuf[1];
  if( invalid_dt(dt) ){
    if(Run.rank ==0)
      cout << "Invalid time step: dt_mhd = " << dt << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  dt_min_rt = Physics.rt[i_rt_tstep]*Run.CFL*rbuf[1];

  if(rbuf[0] < dt_min_rt){
    if( Run.rank== 0 )
      if(Run.verbose >1)
        cout << "Override dt_rad: " << rbuf[0] << ' ' << dt_min_rt << endl;
    rbuf[0] = max(rbuf[0],dt_min_rt);
  }

  dt_min = min(rbuf[0],rbuf[1]);
  dt     = min(rbuf[0],Run.CFL*rbuf[1]);

  if (Run.outcad > 0) {
    double rem =  Run.dt_rem();
    if (rem > 0.0)
      dt = min(dt,rem);
  }

  if( Run.autostep ) Run.dt = dt;

  *cfl = Run.dt/dt_min;
  if (*cfl >= 2.) {
    if(Run.rank ==0)
      cout << "Invalid time step: CFL = " << *cfl << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  
  *rt_upd = (int) (rbuf[0]/Run.dt);
  *rt_upd = imax(1,imin(*rt_upd,Physics.rt[i_rt_update]));

  if( Run.rank==0 ){
    if(Run.verbose >1){ 
      cout << "dt_rad   = " << rbuf[0] << " [" << *rt_upd << ']'<< endl;
      cout << "dt_mhd   = " << rbuf[1] << endl;
    }
    if(Run.verbose >0) 
      cout << "dt [CFL] = " << Run.dt << " [" << *cfl << ']' << endl;
  }

}











