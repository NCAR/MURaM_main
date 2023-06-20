#ifndef __GRID_INCLUDED__
#define __GRID_INCLUDED__

#include "precision.h"
#include "physics.H"
#include "run.H"

class GridData{
public:
  int SZ;
  int NDIM;
  int bufsize;
  int gnnodes;
  int vsize;
  int v_nvar;

  int lrank[3];
  int leftr[3];
  int rightr[3];
  int rank;
  int gsize[3];
  int lsize[3];
  int beg[3];
  int end[3];
  int gbeg[3];
  int gend[3];
  int lbeg[3];
  int lend[3];
  int is_gbeg[3];
  int is_gend[3];
  int ghosts[3];
  int periods[3];
  int stride[3];
  int pardim[3];
  int procs[3];

  real w1[3];
  real w2[3];
  real gxmin[3];
  real gxmax[3];
  real lxmin[3];
  real lxmax[3];
  real dx[3];
  real cellvol;
  real volume;

  cState* Urestart;
  cState* U;
  cState* U0;
  cState* Res;

  // eos
  double* pres;
  double* temp;
  double* ne;
  double* rhoi;
  double* amb;

  // RT
  double* Qtot;
  double* Jtot;
  double* Stot;
  double* Tau;

  // Thin Losses
  double* Qthin;

  // Resisitive and Viscous Heating
  double* Qres;
  double* Qvis;

  // Heat Conduction
  double* sflx0;
  double* sflx;
  double* Rflx;
  double* BgradT;

  // divB Cleaning
  double* divB;
  double* phi;

  // Current (for explicit diffusivity)
  Vector* curlB;

  // Hyperbolic ambipolar diffusion
  Vector* v_amb;
  Vector* v0_amb;
  Vector* R_amb;
  Vector* curlBxB;
  double* Qamb;

  // temporary storage variables
  double* tvar1;
  double* tvar2; 
  double* tvar3;
  double* tvar4;
  double* tvar5;
  double* tvar6;
  double* tvar7;
  double* tvar8;

  GridData();
  ~GridData();

  void Init(const RunData &Run, const PhysicsData &Physics);
  void Show() const;

  // Remember, use only local indexes !!!
  int node(const int i, const int j, const int k) const {
#ifdef DEBUG
    int offset = i*stride[0]+j*stride[1]+k*stride[2];
    if(offset >= bufsize) {
      cerr << "Accessing grid out of bounds: "
	   << offset << " > " << bufsize
	   << endl;
    }
    return offset;
#else
    return i*stride[0]+j*stride[1]+k*stride[2];
#endif
  }
  real coord(const int ind, const int dir) const {
    return lxmin[dir]+(0.5+(ind-lbeg[dir]))*dx[dir];
  }
};

#define FULL_LOOP(G,i,j,k) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lend[1]+(G).ghosts[1];(j)++) \
  for((i)=(G).lbeg[0]-(G).ghosts[0];(i)<=(G).lend[0]+(G).ghosts[0];(i)++)

#define LOCAL_LOOP(G,i,j,k) \
  for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
  for((j)=(G).lbeg[1];(j)<=(G).lend[1];(j)++) \
  for((i)=(G).lbeg[0];(i)<=(G).lend[0];(i)++)

#define XBEG_GHOST_LOOP(G,i,j,k) \
  for((i)=(G).lbeg[0]-(G).ghosts[0];(i)<=(G).lbeg[0]-1;(i)++) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lend[1]+(G).ghosts[1];(j)++)

#define XEND_GHOST_LOOP(G,i,j,k) \
  for((i)=(G).lend[0]+1;(i)<=(G).lend[0]+(G).ghosts[0];(i)++) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lend[1]+(G).ghosts[1];(j)++)

#define YBEG_GHOST_LOOP(G,i,j,k) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lbeg[1]-1;(j)++) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((i)=(G).lbeg[0]-(G).ghosts[0];(i)<=(G).lend[0]+(G).ghosts[0];(i)++)

#define YEND_GHOST_LOOP(G,i,j,k) \
  for((j)=(G).lend[1]+1;(j)<=(G).lend[1]+(G).ghosts[1];(j)++) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((i)=(G).lbeg[0]-(G).ghosts[0];(i)<=(G).lend[0]+(G).ghosts[0];(i)++)

#define ZBEG_GHOST_LOOP(G,i,j,k) \
  for((k)=(G).lbeg[2]-(G).ghosts[2];(k)<=(G).lbeg[2]-1;(k)++) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lend[1]+(G).ghosts[1];(j)++) \
  for((i)=(G).lbeg[0]-(G).ghosts[0];(i)<=(G).lend[0]+(G).ghosts[0];(i)++)

#define ZEND_GHOST_LOOP(G,i,j,k) \
  for((k)=(G).lend[2]+1;(k)<=(G).lend[2]+(G).ghosts[2];(k)++) \
  for((j)=(G).lbeg[1]-(G).ghosts[1];(j)<=(G).lend[1]+(G).ghosts[1];(j)++) \
  for((i)=(G).lbeg[0]-(G).ghosts[0];(i)<=(G).lend[0]+(G).ghosts[0];(i)++)

#endif
