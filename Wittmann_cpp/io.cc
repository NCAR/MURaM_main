#include<stdio.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include <stdlib.h>
#include <pwd.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <ctime>
#include <cstdlib>
#include <cfloat>


double* read_abun(int ncontr){         
  int ii;
  double aa=0.;   
  double *eps0=new double [ncontr];
  FILE *fabu;

  char filename[50];
  sprintf(filename,"abundance.txt");
  if(!(fabu=fopen(filename,"r")))
    fprintf(stderr,"error: %s\n",strerror(errno));

  for(int i=0;i<ncontr;i++){
    //read from file
    fscanf(fabu,"%d",&ii);
    fscanf(fabu,"%lf\n",&aa);
    eps0[i]=aa;
  } 
  
  fclose(fabu);

  return eps0;
}

double *** read_ie(int ncontr){
  double ie=0.0;
  FILE * fion;
  
  double *** param_in = new double ** [ncontr];
  for (int aa = 0;aa<ncontr;aa++)
  {
    int ni = aa+1<10?aa+1:10;
    param_in[aa] = new double * [ni];

    for (int ii =0;ii<ni;ii++)
      param_in[aa][ii] = new double [6];
  }

  char filename[50]; 
  sprintf(filename,"atomic_parameters.txt");
  
  if(!(fion=fopen(filename,"r")))
    fprintf(stderr,"error: %s\n",strerror(errno));

  for (int aa=0;aa<ncontr;aa++)
  {
    int ni = aa+1<10?aa+1:10;
    for (int ii=0;ii < ni;ii++)
    {
      int Z=0, J=0,gg=0,mm=0;
      double ie=0.0,e0=0.0,g0=0.0,ee=0.0;
      fscanf(fion,"%d" ,&Z);
      fscanf(fion,"%d" ,&J);
      fscanf(fion,"%lf",&ie);
      fscanf(fion,"%lf",&g0);
      fscanf(fion,"%lf",&e0);
      fscanf(fion,"%lf ",&ee);
      fscanf(fion,"%d " ,&gg);
      fscanf(fion,"%d \n" ,&mm);
      /* fprintf(stdout, "%d %d %lf %lf %lf %lf %d %d \n",Z,J,ie,g0,e0,ee,gg,mm); */
      param_in[aa][ii][0] = ie;
      param_in[aa][ii][1] = g0;
      param_in[aa][ii][2] = e0;
      param_in[aa][ii][3] = ee;
      param_in[aa][ii][4] = (double) gg;
      param_in[aa][ii][5] = (double) mm;
 
    }
  }
  fclose(fion);

  return param_in;
}
