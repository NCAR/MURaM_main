#ifndef __PRECISION_INCLUDED__
#define __PRECISION_INCLUDED__

#ifdef SINGLE_PRECISION
typedef float  real;
#define REALTYPE MPI_FLOAT
#else
typedef double real;
#define REALTYPE MPI_DOUBLE
#endif

#endif
