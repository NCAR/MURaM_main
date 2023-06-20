#ifndef __ATOM__   // __ATOM__
#define __ATOM__

#include<string.h>
#define Hperg 4.407e23
#define k 1.3806503e-16
#define Rgas 8.314e7
#define C 2.99792458e10
#define HP 6.6260687652E-27
#define PI 3.14159265e0
#define amu 1.66053892173E-24
#define ev 1.60218e-12
#define inz_H 13.59844
#define inz_Hm1 0.754
#define D0_h2 4.478
#define D0_h2p 2.654
#define mh 1.67262158e-24
#define a0 5.2918e-9
#define Ry 2.1798741e-11

using std::string;

class atom{
 public:
  
      int na;
      char elm[2];
      string *ATM,*ATM2;
      int *Z;

      double *W,*mass; 
      double *abu, *perg;
      double muavg;
      double abutot;

      int * nion;
      int * pf_sw;
      double ** chi;
      double ** uu;
      double ** nlvl;
      double ** G0;
      double ** E0;
      double ** Gjk;
      double ** epsjk;
      double ** m;

      atom(int,int,double*, double ***);
      ~atom();
     void partf(int,double,double);
};

class molecb{
 public:
    int nm; 
    double *Y,*dY;
    molecb(int,double);
   ~molecb();
};

int * sort(double * list, int size);
#endif                // __ATOM__

