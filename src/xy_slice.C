#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include <iostream>

using std::cout;
using std::endl;

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);

void xy_slice(const RunData&  Run, const GridData& Grid, 
	      const PhysicsData& Physics) {

  static int ini_flag = 1;

  register int i, j, k, node, ind, nsl, v;

  int ibeg = Grid.lbeg[0], iend = Grid.lend[0];
  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];

  int localsize  = Grid.lsize[0]*Grid.lsize[1];
   
  float* iobuf;

  char filename[128];

  static int nslice; 
  static int* ixpos;
  static int nslvar;

  FILE* fhandle=NULL;

  //MPI_File fhandle_mpi;
  //int offset;
  
  if (ini_flag ==1){
    nslice=Physics.slice[i_sl_xy];
    ixpos = (int*) malloc(nslice*sizeof(int));
    for (i=0;i<nslice;i++){
      ixpos[i] = Physics.xy_lev[i];
    }

    if (Run.rank == 0) {
      if(Run.verbose > 1)cout << "xy_slice: " << nslice << endl;
      for (i=0;i<nslice;i++){
      if(Run.verbose > 1) cout << "xy_slice: " << ixpos[i]<< endl;
      }     
    }

    nslvar = 0;
    for (v=0;v<12;v++){
      if (Physics.xy_var[v] == 1) nslvar+=1;
    }

    ini_flag = 0;
  }

  iobuf = (float*) malloc(nslvar*localsize*sizeof(float));

  for (nsl = 0; nsl<nslice; nsl++){

    if ( (Grid.beg[2] <= ixpos[nsl]+Grid.gbeg[2] ) and 
         (Grid.end[2] >= ixpos[nsl]+Grid.gbeg[2] )){

      for (i=ibeg; i<=iend; i++){
        for (j=jbeg; j<=jend; j++){
          ind  = j-jbeg + (i-ibeg)*Grid.lsize[1];
          k    = Grid.lbeg[2]+ixpos[nsl]+Grid.gbeg[2]-Grid.beg[2];
          node = Grid.node(i,j,k);

	  if (Physics.xy_var[0] == 1){
	    iobuf[ind] = (float) Grid.U[node].d;
	    ind += localsize;
	  }
	  if (Physics.xy_var[1] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.x;
	    ind += localsize;
	  }
	  if (Physics.xy_var[2] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.y; 
	    ind += localsize;
	  }
	  if (Physics.xy_var[3] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.z;
	    ind += localsize;
	  }
	  if (Physics.xy_var[4] == 1){
	    iobuf[ind] = (float) (Grid.U[node].e/Grid.U[node].d);
	    ind += localsize;
	  }
	  if (Physics.xy_var[5] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.x;
	    ind += localsize;
	  }
	  if (Physics.xy_var[6] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.y;  
	    ind += localsize;
	  }
	  if (Physics.xy_var[7] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.z;
	    ind += localsize;
	  }
	  if (Physics.xy_var[8] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].M.sqr());
	    ind += localsize;
	  }
	  if (Physics.xy_var[9] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].B.sqr());
	    ind += localsize;
	  }
	  if (Physics.xy_var[10] == 1){
	    iobuf[ind] = (float) Grid.temp[node];
	    ind += localsize;
	  }
	  if (Physics.xy_var[11] == 1){
	    iobuf[ind] = (float) Grid.pres[node];
	  }
	}
      }

      if(Physics.slice[i_sl_collect] == 0) {
	if(xy_rank == 0) { 
	  sprintf(filename,"%s%s_%04d.%06d",Run.path_2D,"xy_slice",ixpos[nsl],
		  Run.globiter);
	    fhandle=fopen(filename,"w");
	    
	    float header[4];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[1];
	    header[2] = (float) Grid.gsize[0];
	    header[3] = (float) Run.time;
	    fwrite(header,sizeof(float),4,fhandle);
	}
	
	slice_write(Grid,0,iobuf,localsize,nslvar,1,0,fhandle);
	
	if(xy_rank == 0) 
	  fclose(fhandle);
	
      } else {
	if(xy_rank == 0) { 
	  sprintf(filename,"%s_%04d.dat","xy_slice",ixpos[nsl]);
	  fhandle=fopen(filename,"a");
	}
	
	slice_write(Grid,0,iobuf,localsize,nslvar,1,0,fhandle);
	
	if(xy_rank == 0){ 
	  fclose(fhandle);
	  
	  sprintf(filename,"%s_%04d.log","xy_slice",ixpos[nsl]);
	  
      std::fstream fptr;
	  int newfile = 0;
	  fptr.open(filename,std::ios::in);
	  if (!fptr) newfile = 1;
	  fptr.close();
	  
	  fptr.open(filename,std::ios::out|std::ios::app);
	  fptr.precision(10);
	  if (newfile) {      
	    fptr <<  nslvar << ' ' <<  Grid.gsize[0] << ' ' 
		 << Grid.gsize[1] << endl;
	    fptr << Physics.xy_var[0]  << ' ' 
		 << Physics.xy_var[1]  << ' ' 
		 << Physics.xy_var[2]  << ' ' 
		 << Physics.xy_var[3]  << ' ' 
		 << Physics.xy_var[4]  << ' ' 
		 << Physics.xy_var[5]  << ' ' 
		 << Physics.xy_var[6]  << ' ' 
		 << Physics.xy_var[7]  << ' ' 
		 << Physics.xy_var[8]  << ' ' 
		 << Physics.xy_var[9]  << ' ' 
		 << Physics.xy_var[10] << ' ' 
		 << Physics.xy_var[11] << endl;
	  }
	  fptr << Run.globiter << ' ' << Run.time << endl;
	  fptr.close();
	}
      }
    }     
  }

  free(iobuf);
}

