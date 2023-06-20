/* ATMDATB is the modified ATMDAT routine so that it outputs not only U1, U2, U3, but also the derivatives of
the logarithms of said energy levels with respect to temperature DU1, DU2, DU3 (except items 39 to 77) plus 56.
the rest are taken as 0 at the moment. SUPPLY ATOMIC PARAMETERS FOR 30 ELEMENTS 
A.D. WITTMANN, GOETTINGEN (1975) */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "atom.h"

int * sort(double * list, int size)
{
  int * ind = new int [size];
  for(int i=0; i<size; i++)
    ind[i] = i;

  for(int i=0; i<size; i++)
  {
    for(int j=size-1; j>i; j--)
    {
      if(list[j]>list[j-1])
      {
        double swap=list[j-1];
        int swapint = ind[j-1];

        list[j-1]=list[j];
        ind[j-1] = ind[j];

        list[j]= swap;
        ind[j] = swapint;
      }
    }
  }
  return ind;
}

using std::string;
/* Names of the elements in capital letters */
string ATM0[]={"H","HE","LI","BE","B","C","N","O","F","NE",
  "NA","MG","AL","SI","P","S","CL","AR","K","CA","SC","TI","V","CR",
  "MN","FE","CO","NI","CU","ZN"};

/*     Names of the elements in small letters  */
string ATM20[]={"h","he","li","be","b","c","n","o","f","ne",
  "na","mg","al","si","p","s","cl","ar","k","ca","sc","ti","v","cr",
  "mn","fe","co","ni","cu","zn"};

double W0[]={1.008,4.003,6.939,9.012,10.811,12.011,14.007,16.,18.998,
  20.183,22.99,24.312,26.982,28.086,30.974,32.064,35.453,39.948,39.102,
  40.08,44.956,47.90,50.942,51.996,54.938,55.847,58.933,58.71,
  63.54,65.37};

int Z0[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
  24,25,26,27,28,29,30};

atom::atom(int na, int nimax, double* eps0, double *** param0){

  ATM = new string [na];        
  ATM2 = new string [na];        

  Z = new int [na];

  W = new double [na];        
  mass = new double [na];
  perg = new double [na];
  abu = new double [na];
  nion = new int [na];
  pf_sw = new int [na];

  uu = new double * [na];
  nlvl = new double * [na];
  chi = new double * [na];
  G0 = new double * [na];
  E0 = new double * [na];
  epsjk = new double * [na];
  Gjk = new double * [na];
  m = new double * [na];

  /* Sort abundances */
  int * na2 = sort(eps0,30);

  muavg=0.0;
  abutot=0.0;
  for (int i=0;i<na;i++){
    int j = na2[i];

    ATM[i]=ATM0[j];
    ATM2[i]=ATM20[j];
    Z[i] = Z0[j];
    W[i] = W0[j];

    nion[i] = (Z[i]<nimax?Z[i]:nimax)+1;

    if (Z[i]>2)
    {
      for (int ni=0;ni<nion[i]-1;ni++)
        if (param0[j][ni][3] ==-1)
          nion[i]=nion[i]<ni+1?nion[i]:ni+1;

      nion[i] = nion[i]>3?nion[i]:3;
      
      /* if we have data use cardona partition functions */
      
      if ((param0[j][1][3]==-1)||(param0[j][0][3]==-1))
        pf_sw[i] = 0;
      else
        pf_sw[i] = 1;
    }
    else
      pf_sw[i] = 1;

    uu[i] = new double [nion[i]];
    nlvl[i] = new double [nion[i]];
    chi[i] = new double [nion[i]];
    G0[i] = new double [nion[i]];
    E0[i] = new double [nion[i]];
    epsjk[i] = new double [nion[i]];
    Gjk[i] = new double [nion[i]];
    m[i] = new double [nion[i]];

    uu[i][0] = 0.0;
    nlvl[i][0] = 0.0;
    chi[i][0] = 0.0;
    G0[i][0] = 0.0;
    E0[i][0] = 0.0;
    epsjk[i][0] = 0.0;
    Gjk[i][0] = 0.0;
    m[i][0] = 0.0;


    for (int ni=1;ni<nion[i];ni++)
    {
      uu[i][ni] = 0.0;
      nlvl[i][ni] = 0.0;
      chi[i][ni] = param0[j][ni-1][0];
      G0[i][ni] = param0[j][ni-1][1];
      E0[i][ni] = param0[j][ni-1][2]*ev/k;
      epsjk[i][ni] = param0[j][ni-1][3]*1.2398e-4*ev/k;
      Gjk[i][ni] = param0[j][ni-1][4];
      m[i][ni] = param0[j][ni-1][5];
    }
    abu[i]=pow(10.,(eps0[i]-12.));
    mass[i] = W[i]*amu;
    abutot+=abu[i];
    muavg+=abu[i]*W[i];
  }

  muavg/=abutot;

  double invmuav = 1.0/(abutot*muavg*amu);
  for (int i=0;i<na;i++)
  {
    perg[i] = abu[i]*invmuav;
  }

  for (int i=0;i<na;i++)
    fprintf(stdout,"atom %d: Z = %d, abundance = %e, perg = %e, nion = %d, pf = %d \n",i,Z[i], abu[i],perg[i], nion[i], pf_sw[i]);

  fprintf(stdout, "-------------------------------\n");
  fprintf(stdout, "Mean mu = %e \n", muavg);
  delete [] na2;
}
atom::~atom(){
  fprintf(stderr,"atom destructor called\n");
  delete [] ATM;
  delete [] ATM2;
  delete [] W;
  delete [] Z;
  delete [] abu;
  delete [] perg;
  delete [] mass;
  delete [] nion;
  delete [] pf_sw;

  for (int i=0;i<na;i++)
  {
    delete [] uu[i];
    delete [] nlvl[i];
    delete [] chi[i];
    delete [] G0[i];
    delete [] E0[i];
    delete [] epsjk[i];
    delete [] Gjk[i];
    delete [] m[i];
  }
  delete [] uu;
  delete [] nlvl;
  delete [] chi;
  delete [] G0;
  delete [] E0;
  delete [] epsjk;
  delete [] Gjk;
  delete [] m;
  fprintf(stderr,"atom destructor complete\n");
}
