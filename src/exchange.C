#include <mpi.h>
#include "grid.H"
#include "physics.H"
#include "comm_split.H"
#define GLOOP(gs,i,j,d1,d2) \
  for((j)=(gs)[(d2)][0];(j)<=(gs)[(d2)][1];(j)++) \
  for((i)=(gs)[(d1)][0];(i)<=(gs)[(d1)][1];(i)++)

/* 
   Use non-blocking send/recv and overlap communication with writing
   and reading of send/recv buffers!

   contains 3 routines:
   
   exchange_grid(Grid): full exchange of MHD variables + conductive heatflux (when enabled)
   exchange_B(Grid): exchange of B only
   exchange_single(Grid,var): exchange of a single variable

*/ 

void exchange_grid(GridData& Grid, const PhysicsData& Physics, const int needs_sflx_amb_bnd) {

  register int i,j,k,i_nu,buf,dim,d1,d2,d3,v,igh;
  int bfsz;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int perm[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  const int *leftr  = Grid.leftr;
  const int *rightr = Grid.rightr;
  const int myrank  = Grid.rank;
  
  static int size[3],i_next[3];
  static int bfsz_max;
  static int ini_flag=1;

  static double *sndbuf_l, *sndbuf_r, *recbuf_l, *recbuf_r;

  const bool spitzer_enabled   = (Physics.params[i_param_spitzer] > 0.0);
  const bool ambipolar_enabled = (Physics.params[i_param_ambipolar] > 0.0);
  
  if(ini_flag){
   
    for(dim=0;dim<3;dim++){
      size[dim]  = Grid.lsize[dim]+2*Grid.ghosts[dim]; 
    }
    
    i_next[0] = 1;
    i_next[1] = size[0]; 
    i_next[2] = size[0]*size[1];

    bfsz_max=0;
    for(dim=0;dim<Grid.NDIM;dim++){
      d1=perm[dim][0];
      d2=perm[dim][1];
      d3=perm[dim][2];  
      
      bfsz = (Physics.NVAR+spitzer_enabled+ambipolar_enabled*3)*Grid.ghosts[d1]*size[d2]*size[d3];

      bfsz_max = (bfsz_max>bfsz ? bfsz_max : bfsz);
    }

    sndbuf_l = new double[bfsz_max];
    recbuf_l = new double[bfsz_max];
    sndbuf_r = new double[bfsz_max];
    recbuf_r = new double[bfsz_max];
    
    ini_flag=0;
  }

  double* U     = (double*) Grid.U;
  double* v_amb = (double*) Grid.v_amb;

  bool spitzer   = spitzer_enabled and (needs_sflx_amb_bnd == 1);
  bool ambipolar = ambipolar_enabled and (needs_sflx_amb_bnd == 1);

  int gs[3][2];
  for(dim=0;dim<3;dim++){
    gs[dim][0] = Grid.lbeg[dim];
    gs[dim][1] = Grid.lend[dim];
  }

  // optional timing info of data exchange
  //const int timing=0;
  //double ttime,time,btime,data_sent;
  //ttime=MPI_Wtime();
  //btime=0;
  //data_sent=0;

  for(dim=0;dim<Grid.NDIM;dim++){
    d1=perm[dim][0];
    d2=perm[dim][1];
    d3=perm[dim][2];  

    bfsz = (Physics.NVAR+spitzer+ambipolar*3)*Grid.ghosts[d1]*(gs[d2][1]-gs[d2][0]+1)
      *(gs[d3][1]-gs[d3][0]+1);

    //if(timing) data_sent += 2*bfsz*sizeof(double);

    //if(timing) time=MPI_Wtime();
    buf = 0;
    for(igh=0;igh<Grid.ghosts[d1];igh++){
      GLOOP(gs,j,k,d2,d3){
	i=Grid.lbeg[d1]+igh;
	i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	for(v=0;v<Physics.NVAR;v++)
	  sndbuf_l[buf++] = U[i_nu*Physics.NVAR+v];
      }
      if(spitzer){
	GLOOP(gs,j,k,d2,d3){
	  i=Grid.lbeg[d1]+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  sndbuf_l[buf++] = Grid.sflx[i_nu];
	}
      }
      if(ambipolar){
	GLOOP(gs,j,k,d2,d3){
	  i=Grid.lbeg[d1]+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  for(v=0;v<3;v++)
	    sndbuf_l[buf++] = v_amb[i_nu*3+v];
	}
      }
      
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Irecv(recbuf_r,bfsz,MPI_DOUBLE,rightr[d1],0,MPI_COMM_WORLD,rr); 
    MPI_Isend(sndbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 0,MPI_COMM_WORLD,rr+1);

    //if(timing) time=MPI_Wtime();
    buf = 0;
    for(igh=0;igh<Grid.ghosts[d1];igh++){
      GLOOP(gs,j,k,d2,d3){
	i=Grid.lend[d1]-Grid.ghosts[d1]+1+igh;
	i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	for(v=0;v<Physics.NVAR;v++)
	  sndbuf_r[buf++] = U[i_nu*Physics.NVAR+v];
      }
      if(spitzer){
	GLOOP(gs,j,k,d2,d3){
	  i=Grid.lend[d1]-Grid.ghosts[d1]+1+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  sndbuf_r[buf++] = Grid.sflx[i_nu];
	}
      }
      if(ambipolar){
	GLOOP(gs,j,k,d2,d3){
	  i=Grid.lbeg[d1]+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  for(v=0;v<3;v++)
	    sndbuf_r[buf++] = v_amb[i_nu*3+v];
	}
      }
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Irecv(recbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 1,MPI_COMM_WORLD,rl);   
    MPI_Isend(sndbuf_r,bfsz,MPI_DOUBLE,rightr[d1],1,MPI_COMM_WORLD,rl+1);

    MPI_Waitall(2,rr,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gend[d1] or Grid.periods[d1] ) {
      buf = 0;
      for(igh=0;igh<Grid.ghosts[d1];igh++){
	GLOOP(gs,j,k,d2,d3){
	  i=Grid.lend[d1]+1+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  for(v=0;v<Physics.NVAR;v++)  
	    U[i_nu*Physics.NVAR+v] = recbuf_r[buf++];
	}
	if(spitzer){
	  GLOOP(gs,j,k,d2,d3){
	    i=Grid.lend[d1]+1+igh;
	    i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	    Grid.sflx[i_nu] = recbuf_r[buf++];
	  }
	}
	if(ambipolar){
	  GLOOP(gs,j,k,d2,d3){
	    i=Grid.lend[d1]+1+igh;
	    i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	    for(v=0;v<3;v++)
	      v_amb[i_nu*3+v] = recbuf_r[buf++];
	  }
	}
      }
      gs[d1][1]+=Grid.ghosts[d1];
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Waitall(2,rl,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gbeg[d1] or Grid.periods[d1] ) { 
      buf = 0;
      for(igh=0;igh<Grid.ghosts[d1];igh++){
	GLOOP(gs,j,k,d2,d3){
	  i =Grid.lbeg[d1]-Grid.ghosts[d1]+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  for(v=0;v<Physics.NVAR;v++)
	    U[i_nu*Physics.NVAR+v] = recbuf_l[buf++]; 
	}
	if(spitzer){
	  GLOOP(gs,j,k,d2,d3){
	    i =Grid.lbeg[d1]-Grid.ghosts[d1]+igh;
	    i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	    Grid.sflx[i_nu] = recbuf_l[buf++];
	  }
	}
	if(ambipolar){
	  GLOOP(gs,j,k,d2,d3){
	    i=Grid.lend[d1]+1+igh;
	    i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	    for(v=0;v<3;v++)
	      v_amb[i_nu*3+v] = recbuf_l[buf++];
	  }
	}
      }
      gs[d1][0]-=Grid.ghosts[d1];
    }
    //if(timing) btime+=MPI_Wtime()-time;
  }

  //if(timing){ 
  //ttime=MPI_Wtime()-ttime;
  //if(myrank == 0)
  //cout << "exchange: " << ttime << ' ' << btime << ' ' 
  //<< btime/ttime << ' ' << data_sent/ttime/1e6 << endl;
  //}
  
}

//*****************************************************************
void exchange_B(GridData& Grid) {

  const int v0=5;
  const int v1=8;
  const int nvar=8;

  register int i,j,k,i_nu,buf,dim,d1,d2,d3,v,igh;
  int bfsz;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int perm[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  const int *leftr  = Grid.leftr;
  const int *rightr = Grid.rightr;
  const int myrank  = Grid.rank;

  static int size[3],i_next[3];
  static int bfsz_max;
  static int ini_flag=1;

  double* U = (double*) Grid.U;

  static double *sndbuf_l, *sndbuf_r, *recbuf_l, *recbuf_r;

  if(ini_flag){
    for(dim=0;dim<3;dim++){
      size[dim]  = Grid.lsize[dim]+2*Grid.ghosts[dim]; 
    }
    
    i_next[0] = 1;
    i_next[1] = size[0]; 
    i_next[2] = size[0]*size[1];

    bfsz_max=0;
    for(dim=0;dim<Grid.NDIM;dim++){
      d1=perm[dim][0];
      d2=perm[dim][1];
      d3=perm[dim][2];  
      
      bfsz = (v1-v0)*Grid.ghosts[d1]*size[d2]*size[d3];

      bfsz_max = (bfsz_max>bfsz ? bfsz_max : bfsz);

    }
      
    sndbuf_l = new double[bfsz_max];
    recbuf_l = new double[bfsz_max];
    sndbuf_r = new double[bfsz_max];
    recbuf_r = new double[bfsz_max];
    
    ini_flag=0;
  }

  int gs[3][2];
  for(dim=0;dim<3;dim++){
    gs[dim][0] = Grid.lbeg[dim];
    gs[dim][1] = Grid.lend[dim];
  }

  // optional timing info of data exchange
  //const int timing=0;
  //double ttime,time,btime,data_sent;
  //ttime=MPI_Wtime();
  //btime=0;
  //data_sent=0;

  for(dim=0;dim<Grid.NDIM;dim++){
    d1=perm[dim][0];
    d2=perm[dim][1];
    d3=perm[dim][2];  

    bfsz = (v1-v0)*Grid.ghosts[d1]*(gs[d2][1]-gs[d2][0]+1)
      *(gs[d3][1]-gs[d3][0]+1);

    //if(timing) data_sent += 2*bfsz*sizeof(double);

    //if(timing) time=MPI_Wtime();
    buf = 0;
    for(igh=0;igh<Grid.ghosts[d1];igh++){
      GLOOP(gs,j,k,d2,d3){
	i=Grid.lbeg[d1]+igh;
	i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	for(v=v0;v<v1;v++)
	  sndbuf_l[buf++] = U[i_nu*nvar+v];
      }  
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Irecv(recbuf_r,bfsz,MPI_DOUBLE,rightr[d1],0,MPI_COMM_WORLD,rr); 
    MPI_Isend(sndbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 0,MPI_COMM_WORLD,rr+1);

    //if(timing) time=MPI_Wtime();
    buf = 0;
    for(igh=0;igh<Grid.ghosts[d1];igh++){
      GLOOP(gs,j,k,d2,d3){
	i=Grid.lend[d1]-Grid.ghosts[d1]+1+igh;
	i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	for(v=v0;v<v1;v++)
	  sndbuf_r[buf++] = U[i_nu*nvar+v];
      }
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Irecv(recbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 1,MPI_COMM_WORLD,rl);   
    MPI_Isend(sndbuf_r,bfsz,MPI_DOUBLE,rightr[d1],1,MPI_COMM_WORLD,rl+1);

    MPI_Waitall(2,rr,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gend[d1] or Grid.periods[d1] ) {
      buf = 0;
      for(igh=0;igh<Grid.ghosts[d1];igh++){
	GLOOP(gs,j,k,d2,d3){
	  i=Grid.lend[d1]+1+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  for(v=v0;v<v1;v++)  
	    U[i_nu*nvar+v] = recbuf_r[buf++];
	}
      }
      gs[d1][1]+=Grid.ghosts[d1];
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Waitall(2,rl,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gbeg[d1] or Grid.periods[d1] ) { 
      buf = 0;
      for(igh=0;igh<Grid.ghosts[d1];igh++){
	GLOOP(gs,j,k,d2,d3){
	  i =Grid.lbeg[d1]-Grid.ghosts[d1]+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  for(v=v0;v<v1;v++)
	    U[i_nu*nvar+v] = recbuf_l[buf++];
	}
      }
      gs[d1][0]-=Grid.ghosts[d1];
    }
    //if(timing) btime+=MPI_Wtime()-time;
  }

  //if(timing){ 
  //ttime=MPI_Wtime()-ttime;
  //if(myrank == 0)
  //cout << "exchange: " << ttime << ' ' << btime << ' ' 
  //<< btime/ttime << ' ' << data_sent/ttime/1e6 << endl;
  //}

}
//*****************************************************************
void exchange_single(const GridData& Grid, double* var) {

  register int i,j,k,i_nu,buf,dim,d1,d2,d3,igh;
  int bfsz;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int perm[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  const int *leftr  = Grid.leftr;
  const int *rightr = Grid.rightr;
  const int myrank  = Grid.rank;

  static int size[3],i_next[3];
  static int bfsz_max;
  static int ini_flag=1;

  static double *sndbuf_l, *sndbuf_r, *recbuf_l, *recbuf_r;

  if(ini_flag){
    for(dim=0;dim<3;dim++){
      size[dim]  = Grid.lsize[dim]+2*Grid.ghosts[dim]; 
    }
    
    i_next[0] = 1;
    i_next[1] = size[0]; 
    i_next[2] = size[0]*size[1];

    bfsz_max=0;
    for(dim=0;dim<Grid.NDIM;dim++){
      d1=perm[dim][0];
      d2=perm[dim][1];
      d3=perm[dim][2];  
      
      bfsz = Grid.ghosts[d1]*size[d2]*size[d3];

      bfsz_max = (bfsz_max>bfsz ? bfsz_max : bfsz);

    }

    sndbuf_l = new double[bfsz_max];
    recbuf_l = new double[bfsz_max];
    sndbuf_r = new double[bfsz_max];
    recbuf_r = new double[bfsz_max];
    
    ini_flag=0;
  }

  int gs[3][2];
  for(dim=0;dim<3;dim++){
    gs[dim][0] = Grid.lbeg[dim];
    gs[dim][1] = Grid.lend[dim];
  }

  // optional timing info of data exchange
  //const int timing=0;
  //double ttime,time,btime,data_sent;
  //ttime=MPI_Wtime();
  //btime=0;
  //data_sent=0;

  for(dim=0;dim<Grid.NDIM;dim++){
    d1=perm[dim][0];
    d2=perm[dim][1];
    d3=perm[dim][2];  

    bfsz = Grid.ghosts[d1]*(gs[d2][1]-gs[d2][0]+1)
      *(gs[d3][1]-gs[d3][0]+1);

    //if(timing) data_sent += 2*bfsz*sizeof(double);

    //if(timing) time=MPI_Wtime();
    buf = 0;
    for(igh=0;igh<Grid.ghosts[d1];igh++){
      GLOOP(gs,j,k,d2,d3){
	i=Grid.lbeg[d1]+igh;
	i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	sndbuf_l[buf++] = var[i_nu];
      }  
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Irecv(recbuf_r,bfsz,MPI_DOUBLE,rightr[d1],0,MPI_COMM_WORLD,rr); 
    MPI_Isend(sndbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 0,MPI_COMM_WORLD,rr+1);

    //if(timing) time=MPI_Wtime();
    buf = 0;
    for(igh=0;igh<Grid.ghosts[d1];igh++){
      GLOOP(gs,j,k,d2,d3){
	i=Grid.lend[d1]-Grid.ghosts[d1]+1+igh;
	i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	sndbuf_r[buf++] = var[i_nu];
      }
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Irecv(recbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 1,MPI_COMM_WORLD,rl);   
    MPI_Isend(sndbuf_r,bfsz,MPI_DOUBLE,rightr[d1],1,MPI_COMM_WORLD,rl+1);

    MPI_Waitall(2,rr,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gend[d1] or Grid.periods[d1] ) {
      buf = 0;
      for(igh=0;igh<Grid.ghosts[d1];igh++){
	GLOOP(gs,j,k,d2,d3){
	  i=Grid.lend[d1]+1+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  var[i_nu] = recbuf_r[buf++];
	}
      }
      gs[d1][1]+=Grid.ghosts[d1];
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Waitall(2,rl,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gbeg[d1] or Grid.periods[d1] ) { 
      buf = 0;
      for(igh=0;igh<Grid.ghosts[d1];igh++){
	GLOOP(gs,j,k,d2,d3){
	  i =Grid.lbeg[d1]-Grid.ghosts[d1]+igh;
	  i_nu  = i*i_next[d1]+j*i_next[d2]+k*i_next[d3];
	  var[i_nu] = recbuf_l[buf++];
	}
      }
      gs[d1][0]-=Grid.ghosts[d1];
    }
    //if(timing) btime+=MPI_Wtime()-time;
  }

  //if(timing){ 
  //ttime=MPI_Wtime()-ttime;
  //if(myrank == 0)
  //cout << "exchange: " << ttime << ' ' << btime << ' ' 
  //<< btime/ttime << ' ' << data_sent/ttime/1e6 << endl;
  //}

}
//*****************************************************************
void exchange_grid_acc(GridData& Grid, const PhysicsData& Physics, const int needs_sflx_amb_bnd) {

  register int i,j,k,i_nu,buf,buf_sp,buf_amb,dim,d1,d2,d3,v,igh;
  int bfsz;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int perm[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  const int *leftr  = Grid.leftr;
  const int *rightr = Grid.rightr;
  const int myrank  = Grid.rank;
  
  static int size[3],i_next[3];
  static int bfsz_max;
  static int ini_flag=1;

  static double *sndbuf_l, *sndbuf_r, *recbuf_l, *recbuf_r;

  const bool spitzer_enabled   = (Physics.params[i_param_spitzer] > 0.0);
  const bool ambipolar_enabled = (Physics.params[i_param_ambipolar] > 0.0);
  
  if(ini_flag){
   
    for(dim=0;dim<3;dim++){
      size[dim]  = Grid.lsize[dim]+2*Grid.ghosts[dim]; 
    }
    
    i_next[0] = 1;
    i_next[1] = size[0]; 
    i_next[2] = size[0]*size[1];

    bfsz_max=0;
    for(dim=0;dim<Grid.NDIM;dim++){
      d1=perm[dim][0];
      d2=perm[dim][1];
      d3=perm[dim][2];  
      
      bfsz = (Physics.NVAR+spitzer_enabled+ambipolar_enabled*3)*Grid.ghosts[d1]*size[d2]*size[d3];

      bfsz_max = (bfsz_max>bfsz ? bfsz_max : bfsz);
    }

    sndbuf_l = new double[bfsz_max];
    recbuf_l = new double[bfsz_max];
    sndbuf_r = new double[bfsz_max];
    recbuf_r = new double[bfsz_max];
#pragma acc enter data create(sndbuf_l[:bfsz_max])
#pragma acc enter data create(sndbuf_r[:bfsz_max])
#pragma acc enter data create(recbuf_l[:bfsz_max])
#pragma acc enter data create(recbuf_r[:bfsz_max])
    
    ini_flag=0;
  }

  double* U     = (double*) Grid.U;
  double* v_amb = (double*) Grid.v_amb;
  double* sflx  =           Grid.sflx;

  bool spitzer   = spitzer_enabled and (needs_sflx_amb_bnd == 1);
  bool ambipolar = ambipolar_enabled and (needs_sflx_amb_bnd == 1);

  int gs[3][2];
  for(dim=0;dim<3;dim++){
    gs[dim][0] = Grid.lbeg[dim];
    gs[dim][1] = Grid.lend[dim];
  }

  // optional timing info of data exchange
  //const int timing=0;
  //double ttime,time,btime,data_sent;
  //ttime=MPI_Wtime();
  //btime=0;
  //data_sent=0;

  int nvar = Physics.NVAR;
  int bufsize = Grid.bufsize;
#pragma acc loop seq independent
  for(dim=0;dim<Grid.NDIM;dim++){
    d1=perm[dim][0];
    d2=perm[dim][1];
    d3=perm[dim][2];  

    int i_next1 = i_next[d1];
    int i_next2 = i_next[d2];
    int i_next3 = i_next[d3];
    int kbeg = gs[d3][0];
    int kend = gs[d3][1];
    int jbeg = gs[d2][0];
    int jend = gs[d2][1];
    int k_length = kend-kbeg+1;
    int j_length = jend-jbeg+1;
    int lbeg = Grid.lbeg[d1];
    int lend = Grid.lend[d1];
    int ghosts = Grid.ghosts[d1];

    bfsz = (Physics.NVAR+spitzer+ambipolar*3)*ghosts*(gs[d2][1]-gs[d2][0]+1)
      *(gs[d3][1]-gs[d3][0]+1);

    //if(timing) data_sent += 2*bfsz*sizeof(double);

    //if(timing) time=MPI_Wtime();
    buf = k_length*j_length*nvar;
    buf_sp = spitzer*k_length*j_length;
    buf_amb = ambipolar* k_length*j_length*3;
//#pragma acc enter data copyin(i_next1,i_next2,i_next3)

#pragma acc parallel loop collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8],v_amb[:bufsize*3],sflx[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
    for(igh=0;igh<ghosts;igh++){
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
        int buffer = igh*(buf+buf_sp+buf_amb);
	i=lbeg+igh;
	i_nu  = i*i_next1+j*i_next2+k*i_next3;
#pragma acc loop seq
	for(v=0;v<nvar;v++)
	  sndbuf_l[buffer + (k-kbeg)*j_length*nvar + (j-jbeg)*nvar + v] = U[i_nu*nvar+v];

        if(spitzer){
          int buffer_sp = ((igh+1)*buf)+(igh*(buf_sp+buf_amb));
	  sndbuf_l[buffer_sp + (k-kbeg)*j_length + (j-jbeg)] = sflx[i_nu];
	}
        if(ambipolar){
          int buffer_amb = ((igh+1)*(buf+buf_sp))+(igh*(buf_amb));
#pragma acc loop seq
	  for(v=0;v<3;v++)
	    sndbuf_l[buffer_amb + (k-kbeg)*j_length*3 + (j-jbeg)*3 + v] = v_amb[i_nu*3+v];
	}
      }
      
    }
    //if(timing) btime+=MPI_Wtime()-time;

#pragma acc host_data use_device(recbuf_r, sndbuf_l)
{
    MPI_Irecv(recbuf_r,bfsz,MPI_DOUBLE,rightr[d1],0,MPI_COMM_WORLD,rr); 
    MPI_Isend(sndbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 0,MPI_COMM_WORLD,rr+1);
}

    //if(timing) time=MPI_Wtime();
#pragma acc parallel loop collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8],v_amb[:bufsize*3],sflx[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
    for(igh=0;igh<ghosts;igh++){
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
        int buffer = igh*(buf+buf_sp+buf_amb);
	i=lend-ghosts+1+igh;
	i_nu  = i*i_next1+j*i_next2+k*i_next3;
#pragma acc loop seq
	for(v=0;v<nvar;v++)
	  sndbuf_r[buffer + (k-kbeg)*j_length*nvar + (j-jbeg)*nvar + v] = U[i_nu*nvar+v];
      
        if(spitzer){
           int buffer_sp = ((igh+1)*buf)+(igh*(buf_sp+buf_amb));
	   sndbuf_r[buffer_sp + (k-kbeg)*j_length + (j-jbeg)] = sflx[i_nu];
	}
      if(ambipolar){
          int buffer_amb = ((igh+1)*(buf+buf_sp))+(igh*(buf_amb));
#pragma acc loop seq
	  for(v=0;v<3;v++)
	    sndbuf_r[buffer_amb + (k-kbeg)*j_length*3 + (j-jbeg)*3 + v] = v_amb[i_nu*3+v];
	}
      }
    }
    //if(timing) btime+=MPI_Wtime()-time;

#pragma acc host_data use_device(recbuf_l, sndbuf_r)
{
    MPI_Irecv(recbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 1,MPI_COMM_WORLD,rl);   
    MPI_Isend(sndbuf_r,bfsz,MPI_DOUBLE,rightr[d1],1,MPI_COMM_WORLD,rl+1);
}

    MPI_Waitall(2,rr,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gend[d1] or Grid.periods[d1] ) {
#pragma acc parallel loop collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8],v_amb[:bufsize*3],sflx[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
      for(igh=0;igh<ghosts;igh++){
        for(k=kbeg; k<=kend; k++)
        for(j=jbeg; j<=jend; j++) {
          int buffer = igh*(buf+buf_sp+buf_amb);
	  i=lend+1+igh;
	  i_nu  = i*i_next1+j*i_next2+k*i_next3;
#pragma acc loop seq
	  for(v=0;v<nvar;v++)  
	    U[i_nu*nvar+v] = recbuf_r[buffer + (k-kbeg)*j_length*nvar + (j-jbeg)*nvar + v];
	
         if(spitzer){
           int buffer_sp = ((igh+1)*buf)+(igh*(buf_sp+buf_amb));
	    sflx[i_nu] = recbuf_r[buffer_sp + (k-kbeg)*j_length + (j-jbeg)];
	  }
	if(ambipolar){
            int buffer_amb = ((igh+1)*(buf+buf_sp))+(igh*(buf_amb));
#pragma acc loop seq
	    for(v=0;v<3;v++)
	      v_amb[i_nu*3+v] = recbuf_r[buffer_amb + (k-kbeg)*j_length*3 + (j-jbeg)*3 + v];
	  }
	}
      }
      gs[d1][1]+=ghosts;
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Waitall(2,rl,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gbeg[d1] or Grid.periods[d1] ) { 
#pragma acc parallel loop collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8],v_amb[:bufsize*3],sflx[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
      for(igh=0;igh<ghosts;igh++){
        for(k=kbeg; k<=kend; k++)
        for(j=jbeg; j<=jend; j++) {
          int buffer = igh*(buf+buf_sp+buf_amb);
	  i =lbeg-ghosts+igh;
	  i_nu  = i*i_next1+j*i_next2+k*i_next3;
#pragma acc loop seq
	  for(v=0;v<nvar;v++)
	    U[i_nu*nvar+v] = recbuf_l[buffer + (k-kbeg)*j_length*nvar + (j-jbeg)*nvar + v]; 
	
         if(spitzer){
            int buffer_sp = ((igh+1)*buf)+(igh*(buf_sp+buf_amb));
	    sflx[i_nu] = recbuf_l[buffer_sp + (k-kbeg)*j_length + (j-jbeg)];
	  }
	if(ambipolar){
            int buffer_amb = ((igh+1)*(buf+buf_sp))+(igh*(buf_amb));
#pragma acc loop seq
	    for(v=0;v<3;v++)
	      v_amb[i_nu*3+v] = recbuf_l[buffer_amb + (k-kbeg)*j_length*3 + (j-jbeg)*3 + v];
	  }
	}
      }
      gs[d1][0]-=ghosts;
    }
//#pragma acc exit data delete(i_next1,i_next2,i_next3)
    //if(timing) btime+=MPI_Wtime()-time;
  }
  //if(timing){ 
  //ttime=MPI_Wtime()-ttime;
  //if(myrank == 0)
  //cout << "exchange: " << ttime << ' ' << btime << ' ' 
  //<< btime/ttime << ' ' << data_sent/ttime/1e6 << endl;
  //}
  
}

//*****************************************************************
void exchange_B_acc(GridData& Grid) {

  const int v0=5;
  const int v1=8;
  const int nvar=8;

  register int i,j,k,i_nu,buf,dim,d1,d2,d3,v,igh;
  int bfsz;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int perm[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  const int *leftr  = Grid.leftr;
  const int *rightr = Grid.rightr;
  const int myrank  = Grid.rank;

  static int size[3],i_next[3],i_nextd1,i_nextd2,i_nextd3;
  static int bfsz_max;
  static int ini_flag=1;

  double* U = (double*) Grid.U;

  static double *sndbuf_l, *sndbuf_r, *recbuf_l, *recbuf_r;

  if(ini_flag){
    for(dim=0;dim<3;dim++){
      size[dim]  = Grid.lsize[dim]+2*Grid.ghosts[dim]; 
    }
    
    i_next[0] = 1;
    i_next[1] = size[0]; 
    i_next[2] = size[0]*size[1];

    bfsz_max=0;
    for(dim=0;dim<Grid.NDIM;dim++){
      d1=perm[dim][0];
      d2=perm[dim][1];
      d3=perm[dim][2];  
      
      bfsz = (v1-v0)*Grid.ghosts[d1]*size[d2]*size[d3];

      bfsz_max = (bfsz_max>bfsz ? bfsz_max : bfsz);
    }

      sndbuf_l = new double[bfsz_max];
      recbuf_l = new double[bfsz_max];
      sndbuf_r = new double[bfsz_max];
      recbuf_r = new double[bfsz_max];
#pragma acc enter data create(sndbuf_l[:bfsz_max])
#pragma acc enter data create(sndbuf_r[:bfsz_max])
#pragma acc enter data create(recbuf_l[:bfsz_max])
#pragma acc enter data create(recbuf_r[:bfsz_max])
    
    ini_flag=0;
  }

  int gs[3][2];
  for(dim=0;dim<3;dim++){
    gs[dim][0] = Grid.lbeg[dim];
    gs[dim][1] = Grid.lend[dim];
  }

  // optional timing info of data exchange
  //const int timing=0;
  //double ttime,time,btime,data_sent;
  //ttime=MPI_Wtime();
  //btime=0;
  //data_sent=0;

  int bufsize = Grid.bufsize;
#pragma acc data copyin(U[:bufsize*8],sndbuf_l[:bfsz_max], sndbuf_r[:bfsz_max], recbuf_l[:bfsz_max], recbuf_r[:bfsz_max])
{
#pragma acc loop seq independent
  for(dim=0;dim<Grid.NDIM;dim++){
    d1=perm[dim][0];
    d2=perm[dim][1];
    d3=perm[dim][2];  
 
    i_nextd1 = i_next[d1];
    i_nextd2 = i_next[d2];
    i_nextd3 = i_next[d3];

    int kbeg = gs[d3][0];
    int kend = gs[d3][1];
    int jbeg = gs[d2][0];
    int jend = gs[d2][1];
    int k_length = kend-kbeg+1;
    int j_length = jend-jbeg+1;
    int lbeg = Grid.lbeg[d1];
    int lend = Grid.lend[d1];
    int ghosts = Grid.ghosts[d1];

    bfsz = (v1-v0)*Grid.ghosts[d1]*(gs[d2][1]-gs[d2][0]+1)
      *(gs[d3][1]-gs[d3][0]+1);

    //if(timing) data_sent += 2*bfsz*sizeof(double);

    //if(timing) time=MPI_Wtime();
#pragma acc parallel loop gang collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
    for(igh=0;igh<ghosts;igh++){
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
        int buffer = igh*k_length*j_length*(v1-v0);
	i=lbeg+igh;
	i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
#pragma acc loop vector
	for(v=v0;v<v1;v++)
	  sndbuf_l[buffer + (k-kbeg)*j_length*(v1-v0) + (j-jbeg)*(v1-v0) + (v-v0)] = U[i_nu*nvar+v];
      }  
    }
    //if(timing) btime+=MPI_Wtime()-time;

#pragma acc host_data use_device(recbuf_r, sndbuf_l)
{
    MPI_Irecv(recbuf_r,bfsz,MPI_DOUBLE,rightr[d1],0,MPI_COMM_WORLD,rr); 
    MPI_Isend(sndbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 0,MPI_COMM_WORLD,rr+1);
}

    //if(timing) time=MPI_Wtime();
#pragma acc parallel loop gang collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz]) 
    for(igh=0;igh<ghosts;igh++){
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
        int buffer = igh*k_length*j_length*(v1-v0);
	i=lend-ghosts+1+igh;
	i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
#pragma acc loop vector
	for(v=v0;v<v1;v++)
	  sndbuf_r[buffer + (k-kbeg)*j_length*(v1-v0) + (j-jbeg)*(v1-v0) + (v-v0)] = U[i_nu*nvar+v];
      }
    }
    //if(timing) btime+=MPI_Wtime()-time;

#pragma acc host_data use_device(recbuf_l, sndbuf_r)
{
    MPI_Irecv(recbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 1,MPI_COMM_WORLD,rl);   
    MPI_Isend(sndbuf_r,bfsz,MPI_DOUBLE,rightr[d1],1,MPI_COMM_WORLD,rl+1);
}

    MPI_Waitall(2,rr,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gend[d1] or Grid.periods[d1] ) {
#pragma acc parallel loop gang collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
      for(igh=0;igh<ghosts;igh++){
        for(k=kbeg; k<=kend; k++)
        for(j=jbeg; j<=jend; j++) {
          int buffer = igh*k_length*j_length*(v1-v0);
	  i=lend+1+igh;
	  i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
#pragma acc loop vector
	  for(v=v0;v<v1;v++)  
	    U[i_nu*nvar+v] = recbuf_r[buffer + (k-kbeg)*j_length*(v1-v0) + (j-jbeg)*(v1-v0) + (v-v0)];
	}
      }
      gs[d1][1]+=ghosts;
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Waitall(2,rl,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gbeg[d1] or Grid.periods[d1] ) { 
#pragma acc parallel loop gang collapse(3) private(i,i_nu,v) \
 present(U[:bufsize*8], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
      for(igh=0;igh<ghosts;igh++){
        for(k=kbeg; k<=kend; k++)
        for(j=jbeg; j<=jend; j++) {
          int buffer = igh*k_length*j_length*(v1-v0);
	  i =lbeg-ghosts+igh;
	  i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
#pragma acc loop vector
	  for(v=v0;v<v1;v++)
	    U[i_nu*nvar+v] = recbuf_l[buffer + (k-kbeg)*j_length*(v1-v0) + (j-jbeg)*(v1-v0) + (v-v0)];
	}
      }
      gs[d1][0]-=ghosts;
    }
    //if(timing) btime+=MPI_Wtime()-time;
  }
  } //end data
  //if(timing){ 
  //ttime=MPI_Wtime()-ttime;
  //if(myrank == 0)
  //cout << "exchange: " << ttime << ' ' << btime << ' ' 
  //<< btime/ttime << ' ' << data_sent/ttime/1e6 << endl;
  //}

}

//*****************************************************************
void exchange_single_acc(const GridData& Grid, double* var) {

  register int i,j,k,i_nu,buf,dim,d1,d2,d3,igh;
  int bfsz;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int perm[3][3] = {{ 0, 1, 2 },{ 1, 0, 2 },{ 2, 0, 1 }};

  const int *leftr  = Grid.leftr;
  const int *rightr = Grid.rightr;
  const int myrank  = Grid.rank;

  static int size[3],i_next[3],i_nextd1,i_nextd2,i_nextd3;
  static int bfsz_max;
  static int ini_flag=1;

  static double *sndbuf_l, *sndbuf_r, *recbuf_l, *recbuf_r;

  if(ini_flag){
    for(dim=0;dim<3;dim++){
      size[dim]  = Grid.lsize[dim]+2*Grid.ghosts[dim]; 
    }
    
    i_next[0] = 1;
    i_next[1] = size[0]; 
    i_next[2] = size[0]*size[1];

    bfsz_max=0;
    for(dim=0;dim<Grid.NDIM;dim++){
      d1=perm[dim][0];
      d2=perm[dim][1];
      d3=perm[dim][2];  
       
      bfsz = Grid.ghosts[d1]*size[d2]*size[d3];

      bfsz_max = (bfsz_max>bfsz ? bfsz_max : bfsz);

    }
      sndbuf_l = new double[bfsz_max];
      recbuf_l = new double[bfsz_max];
      sndbuf_r = new double[bfsz_max];
      recbuf_r = new double[bfsz_max];
#pragma acc enter data create(sndbuf_l[:bfsz_max])
#pragma acc enter data create(sndbuf_r[:bfsz_max])
#pragma acc enter data create(recbuf_l[:bfsz_max])
#pragma acc enter data create(recbuf_r[:bfsz_max])
    
    ini_flag=0;
  }

  int gs[3][2];
  for(dim=0;dim<3;dim++){
    gs[dim][0] = Grid.lbeg[dim];
    gs[dim][1] = Grid.lend[dim];
  }

  // optional timing info of data exchange
  //const int timing=0;
  //double ttime,time,btime,data_sent;
  //ttime=MPI_Wtime();
  //btime=0;
  //data_sent=0;

  int bufsize = Grid.bufsize;
#pragma acc data copyin(var[:bufsize],sndbuf_l[:bfsz_max], sndbuf_r[:bfsz_max], recbuf_l[:bfsz_max], recbuf_r[:bfsz_max])
{
#pragma acc loop seq independent
  for(dim=0;dim<Grid.NDIM;dim++){
    d1=perm[dim][0];
    d2=perm[dim][1];
    d3=perm[dim][2]; 

    i_nextd1 = i_next[d1];
    i_nextd2 = i_next[d2];
    i_nextd3 = i_next[d3];

    int kbeg = gs[d3][0];
    int kend = gs[d3][1];
    int jbeg = gs[d2][0];
    int jend = gs[d2][1];
    int k_length = kend-kbeg+1;
    int j_length = jend-jbeg+1;
    int lbeg = Grid.lbeg[d1];
    int lend = Grid.lend[d1];
    int ghosts = Grid.ghosts[d1]; 

    bfsz = Grid.ghosts[d1]*(gs[d2][1]-gs[d2][0]+1)
      *(gs[d3][1]-gs[d3][0]+1);

    //if(timing) data_sent += 2*bfsz*sizeof(double);

    //if(timing) time=MPI_Wtime();
    buf = 0;


#pragma acc parallel loop gang vector collapse(3) private(i,i_nu) \
 present(var[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
    for(igh=0;igh< ghosts;igh++){
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
        int buffer = igh*k_length*j_length;
	i=lbeg+igh;
	i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
	sndbuf_l[buffer + (k-kbeg)*j_length + (j-jbeg)] = var[i_nu];
      }  
    }
    //if(timing) btime+=MPI_Wtime()-time;

#pragma acc host_data use_device(recbuf_r, sndbuf_l)
{
    MPI_Irecv(recbuf_r,bfsz,MPI_DOUBLE,rightr[d1],0,MPI_COMM_WORLD,rr); 
    MPI_Isend(sndbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 0,MPI_COMM_WORLD,rr+1);
}

    //if(timing) time=MPI_Wtime();
#pragma acc parallel loop gang vector collapse(3)  private(i,i_nu) \
 present(var[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
    for(igh=0;igh<ghosts;igh++){
      for(k=kbeg; k<=kend; k++)
      for(j=jbeg; j<=jend; j++) {
        int buffer = igh*k_length*j_length;
	i=lend-ghosts+1+igh;
	i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
	sndbuf_r[buffer + (k-kbeg)*j_length + (j-jbeg)] = var[i_nu];
      }
    }
    //if(timing) btime+=MPI_Wtime()-time;

   MPI_Waitall(2,rr,st);
#pragma acc host_data use_device(recbuf_l, sndbuf_r)
{
    MPI_Irecv(recbuf_l,bfsz,MPI_DOUBLE,leftr[d1], 1,MPI_COMM_WORLD,rl);   
    MPI_Isend(sndbuf_r,bfsz,MPI_DOUBLE,rightr[d1],1,MPI_COMM_WORLD,rl+1);
}

   // MPI_Waitall(2,rr,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gend[d1] or Grid.periods[d1] ) {
#pragma acc parallel loop gang vector collapse(3) private(i,i_nu) \
 present(var[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
       for(igh=0;igh<ghosts;igh++){
        for(k=kbeg; k<=kend; k++)
        for(j=jbeg; j<=jend; j++) {
          int buffer = igh*k_length*j_length;
	  i=lend+1+igh;
	  i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
	  var[i_nu] = recbuf_r[buffer + (k-kbeg)*j_length + (j-jbeg)];
	}
      }
      gs[d1][1]+=ghosts;
    }
    //if(timing) btime+=MPI_Wtime()-time;

    MPI_Waitall(2,rl,st);

    //if(timing) time=MPI_Wtime();
    if( !Grid.is_gbeg[d1] or Grid.periods[d1] ) { 
#pragma acc parallel loop gang vector collapse(3) private(i,i_nu) \
 present(var[:bufsize], \
  sndbuf_l[:bfsz], sndbuf_r[:bfsz], recbuf_l[:bfsz], recbuf_r[:bfsz])
       for(igh=0;igh< ghosts;igh++){
        for(k=kbeg; k<=kend; k++)
        for(j=jbeg; j<=jend; j++) {
          int buffer = igh*k_length*j_length;
	  i =lbeg-ghosts+igh;
	  i_nu  = i*i_nextd1+j*i_nextd2+k*i_nextd3;
	  var[i_nu] = recbuf_l[buffer + (k-kbeg)*j_length + (j-jbeg)];
	}
      }
      gs[d1][0]-=ghosts;
    }
    //if(timing) btime+=MPI_Wtime()-time;
  }
  } //acc data present

  //if(timing){ 
  //ttime=MPI_Wtime()-ttime;
  //if(myrank == 0)
  //cout << "exchange: " << ttime << ' ' << btime << ' ' 
  //<< btime/ttime << ' ' << data_sent/ttime/1e6 << endl;
  //}

}

