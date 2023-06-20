#include<stdio.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include "atom.h"

double saha(double theta,double chi,double u1,double u2,double pe){
  //     SAHA-EGGERT EQUATION  
  double saha_fr=u2*pow(10.,(9.0805126-theta*chi))/(u1*pe*pow(theta,2.5) );
  //       fprintf(stderr,"saha_fr=%e\n",saha_fr);
  return saha_fr;
}
