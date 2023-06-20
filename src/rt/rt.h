#ifndef __RT_RT_H__   // __RT_RT_H__
#define __RT_RT_H__

#include <mpi.h>
#include "grid.H"
#include "run.H"
#include "physics.H"

#define RT_DEFAULT    0
#define RT_SCATTER    1

#define NMU 3

#define UP 0
#define DOWN 1
#define UPWIND 0
#define DOWNWIND 1
#define RIGHT 0
#define LEFT 1
#define FWD  0
#define BWD  1

//             ODF-Kram

#define Nlam 328         // Number of intervals in ODF
#define Nbin 12          //  # bins per interval

#define Rgas 8.31451e7
#define PI 3.14159265e0
#define Cvac 2.99792458e10 // cm/s
#define H 6.626196e-27     // erg*sec
#define k_B 1.380622e-16   // erg/K
#define SIG 5.66961e-05
#define SFLUX 6.28e10
#define TENLOG 2.30258509299405e0

#define dtau_min  1.0E-5
#define dtau_min2  1.0E-10
#define threshold 0.001

#define C   2.99792458E10            // cm/s
#define h   6.6260687652E-27         // erg s
#define kb  1.3806503E-16            // erg/K
#define me  9.1093818872E-28         // g
#define mp  1.67262158E-24           // g
#define mn  1.67492716E-24           // g
#define amu  1.66053892173E-24        // g

#define qe  4.8032042510E-10         // statC
#define Ryd 1.097373156853955E5      // cm-1
#define e_rydberg 13.59844*1.602189E-12  //ergs //Ion.pot.H
#define e0  7.957747154594775607E-02 // 1/(4 pi)
#define R   8.31447215E7             // erg mol-1 K-1
#define ev  1.602189E-12             // erg
#define a0  5.2918E-9               //Bohr radius in cm;

#define km_to_cm 1.e+5
#define cm_to_km 1.e-5
#define m_to_cm 1.e+2
#define cm_to_m 1.e-2
  
#define inv3 0.33333333333333333
#define inv6 0.16666666666666666

class RTS{
protected:
  // For Kappa Tables
  int Nbands, NT, Np;
  int bin_ind;
  int nu_ind;

  // Is there a reference bin
  int N5000;
  // NDIM, 1,2,3
  int NDIM;
  /* fullodf - full odf mode, 0,1.
  If 0) then run opacity binned MURaM setup as usual.
  If 1) then for the scattering case we write out files of J and S for each band.
  And for each 12 bins we weight the intensity used to calculate J's and F's
  based on the ODF. w_sbin above. */
  int fullodf;
  int rttype;
  int scatter;
  int Npp;

  double eps_const;

  // Do i need to output intensity for a slice write?
  int need_I;
  // Use line or optically thin cooling
  int ELTE_on;
// MPI
  int myrank,verbose,lrank[3];
  int cart_sizes[3];
  int leftr[3],rightr[3];
  MPI_Comm  **comm_col;

// dimensions
  int xl,xh,nx,xo;
  int yl,yh,ny,yo;
  int zl,zh,nz,zo;
// driver
  int ***** numits;
  double *I_o, *I_n, *I_n1;
  int ibase[NMU],ixstep[4],iystep[4],izstep[4];
  double *Fx,*Fy,*Fz;
  double a_00[3][NMU],a_01[3][NMU],a_10[3][NMU],a_11[3][NMU];
  double *Tau;
  double  xmu[3][NMU],wmu[NMU];
  double  ds_upw[NMU];//, ds_dnw[NMU] ,dz_upw, dz_dnw;
  double *coeff,*coeff1,*coeff2,*dc;
  //double dc0,dc1,dc2,dc3;
  //int dsize;
// wrapper
  double F_o,dt_rad;
  double * Fr_mean, * gFr_mean;
  /* Column outputs
   * 0: J_col - Angle averaged intensity
   * 1: S_col - Total Source function (=B in LTE)
   * 2: kap_col - total opacity
   * 3: abs_col - absorption opacity (= kap in LTE)
   * 4: sig_col - scattering opacity (= sig in LTE)
   * 5: B_col - Planck function (=S in LTE)
   * 6: tau_col - band dependant tau
   * 7: Qj_col - Q from radiative energy imbalance
   * 8: Qf_col - Q from flux divergence
   */
  int save_col, num_col;
  double avg_col;
  int col_nz, col_nzt, col_offz, col_nvar, col_bnd[4];
  double *** Col_out;
 
  // Multi-dim arrays
  double *Qt,*St,*Jt;
  double ** sbuf;
  double ** rbuf;
  double * Qtemp;

  double *rho,*lgTe,*lgPe,*ne;
  double *kap, *abn, *sig;
  double *B,*J_band;
  int * T_ind, * P_ind;

  // Stuff for Kappa Tables
  // band-P-T tabulated opacities
  double * tab_T;
  double * invT_tab;
  double * tab_p;
  double * invP_tab;
  float * kap_tab;
  float *** sig_tab;
  float *** abn_tab;

  // band-T tabulated Source function
  float * B_tab;

  // For Full ODF
  float * nu_tab;
  float *** acont_pT;
  float *** kcont_pT;

  // Reference wavelength grid
  float * B_5000_tab;
  float * kap_5000_tab;

  // Plane Parallel band-tau5000 tabulated opacities
  double * tau_pp_tab;
  double * invtau_pp_tab;
  float ** kap_pp_tab;
  float ** abn_pp_tab;
  float ** sig_pp_tab;

  //
  double  delta_t1,delta_t2;
  real *x_oldbuf,*y_oldbuf,*z_oldbuf;
  real *x_sbuf,*y_sbuf,*z_sbuf;
  real *x_rbuf,*y_rbuf,*z_rbuf;
  int isgbeg[3],isgend[3],next[3];
  int call_count;
  int fail_count;
  int total_iter;

  double * I_band;

  int *tr_switch;
  void driver(double DZ, double DX, double DY, int band); 
  void interpol(int,int,int,int,int,int,int,int,int,int,double*,double*);
  void IntegrateSetup(int yi_i, int xi_i, int zi_i, int ystep, int xstep, int zstep);
  void integrate(const double c[4]);
  void interpol_and_integrate(const double c[4], double * Ss, int l);
  double error(int,int,int,int,int,double);
  void readbuf(int band,int l,int  DIR,int XDIR,int YDIR);
  void writebuf(int band, int l,int DIR,int XDIR,int YDIR);
  void exchange(int band,int l,int DIR,int XDIR,int YDIR);
  void flux(int l,int DIR,int XDIR,int YDIR);
  void qrad(const double DZ,const double DX,const double DY,const int band);
  void tauscale_qrad(int band, double DX,double DY,double DZ,double* Ss);
  void calc_Qtot_and_Tau(GridData&, const RunData&, const PhysicsData&);
  void get_Tau_and_Iout(GridData&, const RunData&, const PhysicsData&, double DZ, float * B_Iout_tab, float * kap_Iout_tab, double * I_band, int calc_Int);
  void load_bins(char *);
  void save_1D_avg(char*,int,double);

  double * tI_n;
  double *trho, *tkap, *tB;

public:
  RTS(GridData & Grid, RunData & Run, PhysicsData & Physics);
  virtual ~RTS(void);
  virtual double wrapper(int,GridData&,RunData&,const PhysicsData&);
  double tau(int,int,int); 
  double Qtot(int,int,int); 
  double Stot(int,int,int);
  double Jtot(int,int,int);
  double Iout(int,int); 
  double Fout(void){ return F_o; }
  void UpdateIout();
  int grey_rt(void){ return 1; } //Nbands == 1? 1:0;
};


RTS *rt_new(GridData&,RunData&,PhysicsData&);

bool CheckDependency(const int step[4]);
void Transpose_Rho_Kap_B(
  double * rho, double * trho,
  double * kap, double * tkap,
  double * B,   double * tB,
  const int nx, const int ny, const int nz
);
void Transpose_In(double * I_n, double * tI_n, const int nx, const int ny, const int nz);
void Untranspose_In(double * I_n, double * tI_n, const int nx, const int ny, const int nz);
void Transpose_integrate_kernel(
  double * I_n, double * coeff,
  const double c[4], const int off[4],
  const int bounds[3], const int str[3],
  const int strc[3], const int i[3],
  const int step, const int n
);
void Transpose_integrate(
  double * I_n, double * tI_n, double * coeff,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR,
  const int ixstep[4], const int iystep[4], const int izstep[4],
  int x_i, int y_i, int z_i,
  const double c[4]
);
void Transpose_interpol_kernel(
  double * rho, double * kap, double * Ss, double * coeff,
  const double c[4], const int off[4],
  const int bounds[3], const int str[3],
  const int strc[3], const int i[3],
  const int step, const int n,
  const double ds3, const double ds6
);
void Transpose_interpol(
  double * rho, double * trho,
  double * kap, double * tkap,
  double * Ss, double * tSs, double * coeff,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l,
  const int ixstep[4], const int iystep[4], const int izstep[4],
  int x_i, int y_i, int z_i, const double c[4], const double ds_upw[3]
);
void Transpose_readbuf_kernel(
  real * xsb, real * xrb, real * xob,
  real * ysb, real * yrb, real * yob,
  real * zsb, real * zrb, real * zob,
  double * I_n, const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR,
  const int strx, const int stry, const int strz
);
void Transpose_readbuf(
  real * xsb, real * xrb, real * xob,
  real * ysb, real * yrb, real * yob,
  real * zsb, real * zrb, real * zob,
  double * I_n, double * tI_n, const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l, const int band,
  const int ixstep[4], const int iystep[4], const int izstep[4]
);
void Transpose_writebuf_kernel(
  real * xsb, real * ysb, real * zsb, double * I_n,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR,
  const int strx, const int stry, const int strz
);
void Transpose_writebuf(
  real * xsb, real * ysb, real * zsb, double * I_n, double * tI_n,
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l, const int band,
  int ixstep[4], int iystep[4], int izstep[4]
);
void Transpose_flux_kernel(
  double * I_n, double * Fx, double * Fy, double * Fz, double * J_band,
  const double c_J, const double c_x, const double c_y, const double c_z,
  const int nx, const int ny, const int nz
);
void Transpose_flux(
  double * I_n, double * tI_n, double * Fx, double * Fy, double * Fz, double * J_band,
  double wmu[3], double xmu[3][3],
  const int nx, const int ny, const int nz,
  const int XDIR, const int YDIR, const int ZDIR, const int l,
  int ixstep[4], int iystep[4], int izstep[4]
);

#endif                // __RT_RT_H__
