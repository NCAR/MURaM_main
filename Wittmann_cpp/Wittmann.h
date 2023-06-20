#ifndef __IONS_WITTMANN__   // __IONS_WITTMANN__
#define __IONS_WITTMANN__

#include "atom.h"

/* Print more information y/n */
const int debug_verbose = 0;

/* number of atoms and molecules */
int const ncontr=15;
int const nimax=5;
int const nmol=2;

/* Tolerance of electron number iteration, and of temperature/energy
 * serach iteration, as well as max no. of iterations */
const double pe_tol = 1.0e-10;
const double search_tol = 1.0e-10;
const int itermax = 1e7;

/* internal energy offset from MURaM EOS, if required, will need to calculate more accurately
 * Energy calculated to match pressure of OPAL at eps=2.0e12,rho=1.0e-8
 * entropy offset calculated to match at the same point */
const double eps_off =  7.51e10;
const double ss_off =  -2066949262.257303;

/* EOS variables - tables, axis, number of points, spacing, inverse spacing */
float ** ttbl, ** ptbl, ** netbl, ** rhoitbl, ** rhontbl, ** ambtbl, ** stbl, ** rhotbl, ** epstbl;
double * eps_grid, * r_grid, *s_grid, *p_grid;
int Ns,Nr,Neps,Np;
double delta_eps,delta_r,delta_s,delta_p;
double inv_delta_r,inv_delta_eps,inv_delta_s,inv_delta_p;

double * pe_pg10(const atom & at,double,double,double, double);
double saha(double,double,double,double,double);
double acota(double,double,double);
double acotasig(double,double,double);
double * read_abun(int);
double *** read_ie(int);
double * Heps(double tt);
double * invert_eos_newton(double p0,double s0,double rho, double eps,double pmin,double pmax,double smin,double smax,double emin,double emax,double rmin,double rmax,double fct, int maxiter, double tol);
double bilinear(int n1, int n2,double *x,double *y,float **fxy,double xx,double yy);

double p_interp(double ee, double dd);
double s_interp(double ee, double dd);

#endif                // __IONS_WITTMANN__

