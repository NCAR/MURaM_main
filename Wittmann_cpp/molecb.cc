#include<stdio.h>
#include<string.h>
#include<math.h>
#include<iostream>

#include "atom.h"

molecb::molecb(int nm,double X){
  double dX;
  // Molecular partition functions
  Y=new double [nm+1];
  dY=new double [nm+1];
  Y[0] = 0.0 ;
  dY[0] = 0.0;
  Y[1]=-11.206998+X*(2.7942767+X*(7.9196803E-2-X*2.4790744E-2));
  Y[2]=-12.533505+X*(4.9251644+X*(-5.6191273E-2+X*3.2687661E-3));

  dX=(-X*X)/5040.;
  dY[1]=dX*(2.7942767+X*(2*7.9196803e-2-X*3*2.4790744e-2));
  dY[2]=dX*(4.9251644+X*(-2*5.6191273e-2-X*3*3.2687661e-3));
}

molecb::~molecb(){
  delete [] Y;
  delete [] dY;
}


