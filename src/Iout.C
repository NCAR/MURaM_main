#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "eos.H"
#include "comm_split.H"
#include "rt/rt.h"

using namespace std;

extern double get_Iout(int,int);

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);
//======================================================================
void Iout(const RunData&  Run, const GridData& Grid, 
	  const PhysicsData& Physics,RTS *rts) {

  register int ind;

  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
  int kbeg = Grid.lbeg[2], kend = Grid.lend[2];

  int localsize  = Grid.lsize[1]*Grid.lsize[2];

  float* icloc;  

  char filename[128];

  FILE* fhandle=NULL;

  //MPI_File fhandle_mpi;
  //int offset;

  if(Grid.is_gend[0]){  
    icloc = (float*) malloc(localsize*sizeof(float));
    rts->UpdateIout();
    for (int k=kbeg; k<=kend; k++){
      for (int j=jbeg; j<=jend; j++){
	    ind = j-jbeg + (k-kbeg)*Grid.lsize[1];       
	    icloc[ind] = (float) rts->Iout(j,k);
      }
    }

    if(Physics.slice[i_sl_collect] == 0) {
      if(yz_rank == 0) {  
        sprintf(filename,"%s%s.%06d",Run.path_2D,"I_out",Run.globiter);
        fhandle=fopen(filename,"w");
        
        float header[4];            
        header[0] = (float) 1;
        header[1] = (float) Grid.gsize[1];
        header[2] = (float) Grid.gsize[2];
        header[3] = (float) Run.time;
        fwrite(header,sizeof(float),4,fhandle);
      }
      
      slice_write(Grid,0,icloc,localsize,1,1,2,fhandle);
      
      if(yz_rank == 0)
	fclose(fhandle);
      
    } else {
      if(yz_rank == 0){
	sprintf(filename,"I_out.dat");
	fhandle=fopen(filename,"a");
      }
      
      slice_write(Grid,0,icloc,localsize,1,1,2,fhandle);
      
      if(yz_rank == 0){
	fclose(fhandle);
	
	fstream fptr;
	int newfile = 0;
	fptr.open("I_out.log",ios::in);
	if (!fptr) newfile = 1;
	fptr.close();
	
	fptr.open("I_out.log",ios::out|ios::app);
	fptr.precision(10);
	if (newfile)       
	  fptr << '1' << ' ' <<  Grid.gsize[1] << ' ' 
	       << Grid.gsize[2] << endl;
	fptr << Run.globiter << ' ' << Run.time << endl;
	fptr.close();
      }
    }
    free(icloc);
  }
}
