#include <mpi.h>
#include <iostream>
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "ACCH.h"
using namespace std;

GridData::GridData() {
#ifdef NDIM_fixed
  NDIM = NDIM_fixed;
#else
  NDIM = 3;
#endif

  v_nvar = 8;
  bufsize = 1;
  gnnodes = 1;
  cellvol = 1.;
  volume = 1.;
  vsize = 1;

  for(int i=0;i<3;i++) {
    gxmin[i] = 0.;
    gxmax[i] = 1.;
    gsize[i] = 1;
    ghosts[i] = 0;
    periods[i] = 0;
    pardim[i] = 0;
    procs[i] = 1;
    w1[i] = 0.0;
    w2[i] = 0.0;
  }
  
  U   = NULL;
  U0  = NULL;
  Res = NULL;
  
  //eos
  pres = NULL;
  temp = NULL;

  // RT
  Jtot = NULL;
  Stot = NULL;
  Qtot = NULL;
  Tau = NULL;

  // Thin Losses
  
  Qthin = NULL;

  // Resistive and Viscous Heating
  Qres = NULL;
  Qvis = NULL;

  // Heat Conduction
  sflx0 = NULL;
  sflx = NULL;
  Rflx = NULL;
  BgradT = NULL;

  // DivB Cleaning
  divB = NULL;
  phi = NULL;

  // Current
  curlB = NULL;
  
  // ambipolar diffusion
  v_amb  = NULL;
  v0_amb = NULL;
  R_amb  = NULL;
  curlBxB = NULL;
  Qamb   = NULL;
  
  // Electron Number
  ne = NULL;

  // Collisional Frequencies
  amb = NULL;
  rhoi = NULL;

  // Temporary arrays for stuff
  tvar1 = NULL;
  tvar2 = NULL;
  tvar3 = NULL;
  tvar4 = NULL;
  tvar5 = NULL;
  tvar6 = NULL;
  tvar7 = NULL;
  tvar8 = NULL;
  
}

GridData::~GridData() {
    ACCH::Free(U, bufsize*sizeof(cState));
    ACCH::Free(U0, bufsize*sizeof(cState));
    ACCH::Free(Res, bufsize*sizeof(cState));
    
    ACCH::Free(pres, bufsize*sizeof(double));
    ACCH::Free(temp, bufsize*sizeof(double));

    ACCH::Free(ne, bufsize*sizeof(double));
    ACCH::Free(amb, bufsize*sizeof(double));
    ACCH::Free(rhoi, bufsize*sizeof(double));

    ACCH::Free(Qtot, bufsize*sizeof(double));
    ACCH::Free(Stot, bufsize*sizeof(double));
    ACCH::Free(Jtot, bufsize*sizeof(double));
    ACCH::Free(Tau, bufsize*sizeof(double));

    ACCH::Free(Qthin, bufsize*sizeof(double));

    ACCH::Free(Qres, bufsize*sizeof(double));
    ACCH::Free(Qvis, bufsize*sizeof(double));

    ACCH::Free(divB, bufsize*sizeof(double));
    ACCH::Free(phi, bufsize*sizeof(double));

    ACCH::Free(sflx0, bufsize*sizeof(double));
    ACCH::Free(sflx, bufsize*sizeof(double));
    ACCH::Free(Rflx, bufsize*sizeof(double));
    ACCH::Free(BgradT, bufsize*sizeof(double));

    ACCH::Free(curlB, bufsize*sizeof(Vector));

    ACCH::Free(v_amb, bufsize*sizeof(Vector));
    ACCH::Free(v0_amb, bufsize*sizeof(Vector));
    ACCH::Free(R_amb, bufsize*sizeof(Vector));
    ACCH::Free(curlBxB, bufsize*sizeof(Vector));

    ACCH::Free(Qamb, bufsize*sizeof(double));

    ACCH::Free(tvar1, bufsize*sizeof(double));
    ACCH::Free(tvar2, bufsize*sizeof(double));
    ACCH::Free(tvar3, bufsize*sizeof(double));
    ACCH::Free(tvar4, bufsize*sizeof(double));
    ACCH::Free(tvar5, bufsize*sizeof(double));
    ACCH::Free(tvar6, bufsize*sizeof(double));
    ACCH::Free(tvar7, bufsize*sizeof(double));
    ACCH::Free(tvar8, bufsize*sizeof(double));

    ACCH::Delete(this, sizeof(GridData));
}

void GridData::Init(const RunData &Run,const PhysicsData &Physics) {
  MPI_Datatype MPI_STATE;

  MPI_Type_contiguous(sizeof(cState),MPI_BYTE,&MPI_STATE);
  MPI_Type_commit(&MPI_STATE);
  
  int reorder = 0,ndim=3;
  MPI_Cart_create(MPI_COMM_WORLD,3,procs,periods,reorder,&cart_comm);
  MPI_Cart_coords(cart_comm,Run.rank,ndim,lrank);

  // who am I?
  rank = Run.rank;
 
  // find my neighbors
  int ierr[3];
  for (int dim=0 ; dim<3 ; dim++)
    ierr[dim]=MPI_Cart_shift(cart_comm,dim,1,&leftr[dim],&rightr[dim]);

  // Calculate size of local grid nx/ncores_x etc
  lsize[0] = (int) gsize[0]/procs[0];
  lsize[1] = (int) gsize[1]/procs[1];
  lsize[2] = (int) gsize[2]/procs[2];

  // My beginning point in the grid
  beg[0] = ghosts[0]+lrank[0]*lsize[0];
  beg[1] = ghosts[1]+lrank[1]*lsize[1];
  beg[2] = ghosts[2]+lrank[2]*lsize[2];

  // any remaining points
  int remx = gsize[0]-lsize[0]*procs[0];
  int remy = gsize[1]-lsize[1]*procs[1];
  int remz = gsize[2]-lsize[2]*procs[2];

  // distribute remaining points
  for (int ii=0;ii<min(remx,lrank[0]);ii++)
    beg[0]+=1;
  if (lrank[0]<remx)
    lsize[0]+=1;

  for (int ii=0;ii<min(remy,lrank[1]);ii++)
    beg[1]+=1;
  if (lrank[1]<remy)
    lsize[1]+=1;

  for (int ii=0;ii<min(remz,lrank[2]);ii++)
    beg[2]+=1;
  if (lrank[2]<remz)
    lsize[2]+=1;

  // Strides
  stride[0] = 1;
  stride[1] = (lsize[0]+2*ghosts[0]);
  stride[2] = stride[1]*(lsize[1]+2*ghosts[1]);
  
  // Beginning of local grid
  lbeg[0] = ghosts[0];
  lbeg[1] = ghosts[1];
  lbeg[2] = ghosts[2];

  // Beginning of full grid
  gbeg[0] = ghosts[0];
  gbeg[1] = ghosts[1];
  gbeg[2] = ghosts[2];

  // Am I the lower boundary
  if (lrank[0]==0)
    is_gbeg[0] = 1;
  else
    is_gbeg[0] = 0;

  if (lrank[1]==0)
    is_gbeg[1] = 1;
  else
    is_gbeg[1] = 0;

  if (lrank[2]==0)
    is_gbeg[2] = 1;
  else
    is_gbeg[2] = 0;
  
  // My end point in the grid
  end[0] = beg[0] + lsize[0] - 1;
  end[1] = beg[1] + lsize[1] - 1;
  end[2] = beg[2] + lsize[2] - 1;

  // My local end point
  lend[0] = lsize[0] + ghosts[0]-1;
  lend[1] = lsize[1] + ghosts[1]-1;
  lend[2] = lsize[2] + ghosts[2]-1;

  // Full grid end point
  gend[0] = gsize[0]+ghosts[0]-1;
  gend[1] = gsize[1]+ghosts[1]-1;
  gend[2] = gsize[2]+ghosts[2]-1;

  // Am I the upper boundary
  if (lrank[0]==procs[0]-1)
    is_gend[0] = 1;
  else
    is_gend[0] = 0;

  if (lrank[1]==procs[1]-1)
    is_gend[1] = 1;
  else
    is_gend[1] = 0;

  if (lrank[2]==procs[2]-1)
    is_gend[2] = 1;
  else
    is_gend[2] = 0;

  for(int d=0;d<3;d++) {
    int vs   = lsize[d]+2*ghosts[d];
    bufsize *= vs;
    gnnodes *= gsize[d];
    volume *= gxmax[d]-gxmin[d];
    vsize = vsize > vs ? vsize : vs;
  }

  #ifdef CACHE_SIZE
  if(vsize > CACHE_SIZE) {
    std::cerr << "CACHE_SIZE needs to at least be " << vsize << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  #else
  ACCH::shared_size = vsize;
  #endif

  ACCH::Create(this, sizeof(GridData));

  U =   (cState*) ACCH::Malloc(bufsize*sizeof(cState));
  U0 =  (cState*) ACCH::Malloc(bufsize*sizeof(cState));
  Res = (cState*) ACCH::Malloc(bufsize*sizeof(cState));

  pres = (double*) ACCH::Malloc(bufsize*sizeof(double));
  temp = (double*) ACCH::Malloc(bufsize*sizeof(double));

  divB = (double*) ACCH::Malloc(bufsize*sizeof(double));
  phi  = (double*) ACCH::Malloc(bufsize*sizeof(double));

  ne    = (double*) ACCH::Malloc(bufsize*sizeof(double));
  rhoi  = (double*) ACCH::Malloc(bufsize*sizeof(double));
  amb   = (double*) ACCH::Malloc(bufsize*sizeof(double));
 
  Qtot = (double*) ACCH::Malloc(bufsize*sizeof(double));
  Jtot = (double*) ACCH::Malloc(bufsize*sizeof(double));
  Stot = (double*) ACCH::Malloc(bufsize*sizeof(double));
  Tau  = (double*) ACCH::Malloc(bufsize*sizeof(double));

  if(Physics.rt_ext[i_ext_cor]>=1) {
   Qthin = (double*) ACCH::Malloc(bufsize*sizeof(double));
  }
    
  if(Physics.params[i_param_spitzer] > 0.0){
    sflx0  = (double*) ACCH::Malloc(bufsize*sizeof(double));
    sflx   = (double*) ACCH::Malloc(bufsize*sizeof(double));
    Rflx   = (double*) ACCH::Malloc(bufsize*sizeof(double));
    BgradT = (double*) ACCH::Malloc(bufsize*sizeof(double));
  }

  if(Physics.params[i_param_eta] > 0.0){
    curlB = (Vector*) ACCH::Malloc(bufsize*sizeof(Vector));
  }

  if(Physics.params[i_param_ambipolar] > 0.0){
    v_amb   = (Vector*) ACCH::Malloc(bufsize*sizeof(Vector));
    v0_amb  = (Vector*) ACCH::Malloc(bufsize*sizeof(Vector));
    R_amb   = (Vector*) ACCH::Malloc(bufsize*sizeof(Vector));
    curlBxB = (Vector*) ACCH::Malloc(bufsize*sizeof(Vector));
  }

  if(Run.diagnostics){
    tvar1 = (double*) ACCH::Malloc(bufsize*sizeof(double));
    tvar2 = (double*) ACCH::Malloc(bufsize*sizeof(double));
    tvar3 = (double*) ACCH::Malloc(bufsize*sizeof(double));
    tvar4 = (double*) ACCH::Malloc(bufsize*sizeof(double));
    tvar5 = (double*) ACCH::Malloc(bufsize*sizeof(double));
    tvar6 = (double*) ACCH::Malloc(bufsize*sizeof(double));
    tvar7 = (double*) ACCH::Malloc(bufsize*sizeof(double));
    tvar8 = (double*) ACCH::Malloc(bufsize*sizeof(double));

    Qres = (double*) ACCH::Malloc(bufsize*sizeof(double));
    Qvis = (double*) ACCH::Malloc(bufsize*sizeof(double));
   
    if(Physics.params[i_param_ambipolar] > 0.0) {
      Qamb = (double*) ACCH::Malloc(bufsize*sizeof(double));
    }
  }
    
  for(int d=0;d<3;d++) {
    dx[d] = (gxmax[d]-gxmin[d])/gsize[d];

    w1[d] = 8./(12.*dx[d]);
    w2[d] =-1./(12.*dx[d]);
    lxmin[d] = gxmin[d]+(  beg[d]-gbeg[d])*dx[d];
    lxmax[d] = gxmin[d]+(1+end[d]-gbeg[d])*dx[d];
    cellvol *= dx[d];
  }

  ACCH::UpdateGPU(this, sizeof(GridData));

}

void GridData::Show() const {
  cout << " ------------ Grid Parameter Settings -------------" << endl;
  cout << "Decomposition = " << procs[0] << 'x' << procs[1] << 'x' << procs[2] << endl
       << "gsize         = " << gsize[0] << ' ' << gsize[1] << ' ' << gsize[2] << endl
       << "lsize         = " << lsize[0] << ' ' << lsize[1] << ' ' << lsize[2] << endl
       << "ghosts        = " << ghosts[0] << ' ' << ghosts[1] << ' ' << ghosts[2] << endl
       << "beg           = " << beg[0] << ' ' << beg[1] << ' ' << beg[2] << endl
       << "end           = " << end[0] << ' ' << end[1] << ' ' << end[2] << endl
       << "lbeg          = " << lbeg[0] << ' ' << lbeg[1] << ' ' << lbeg[2] << endl
       << "lend          = " << lend[0] << ' ' << lend[1] << ' ' << lend[2] << endl
       << "gbeg          = " << gbeg[0] << ' ' << gbeg[1] << ' ' << gbeg[2] << endl
       << "gend          = " << gend[0] << ' ' << gend[1] << ' ' << gend[2] << endl
       << "stride        = " << stride[0] << ' ' << stride[1] << ' ' << stride[2] << endl
       << "gxmin         = " << gxmin[0] << ' ' << gxmin[1] << ' ' << gxmin[2] << endl
       << "gxmax         = " << gxmax[0] << ' ' << gxmax[1] << ' ' << gxmax[2] << endl
       << "dx            = " << dx[0] << ' ' << dx[1] << ' ' << dx[2] << endl
       << "lxmin         = " << lxmin[0] << ' ' << lxmin[1] << ' ' << lxmin[2] << endl
       << "lxmax         = " << lxmax[0] << ' ' << lxmax[1] << ' ' << lxmax[2] << endl;
}
