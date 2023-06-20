#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<errno.h>
#include<iostream>


double* leeabun(int ncontr){         
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
