#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "rt/rt.h"
#include <iostream>

using std::cout;
using std::endl;

typedef double realtype;

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);

//======================================================================
void tau_slice(const RunData&  Run, const GridData& Grid, 
	      const PhysicsData& Physics,RTS *rts) {

  static int ini_flag = 1;

  const int iroot = 0;

  register int i, j, k, ind, nsl, v, ind1, node1, node2;

  int ibeg = Grid.lbeg[0], iend = Grid.lend[0];
  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
  int kbeg = Grid.lbeg[2], kend = Grid.lend[2];

  int localsize  = Grid.lsize[1]*Grid.lsize[2];
   
  float* iobuf;
  float* iosum;

  char filename[128];

  double q1,q2;

  static int nslice; 
  static double* tau_lev;
  static int nslvar;

  FILE* fhandle=NULL;

  //MPI_File fhandle_mpi;
  //int offset;
 
  if (ini_flag ==1){
    nslice=Physics.slice[i_sl_tau];
    tau_lev = (realtype*) malloc(nslice*sizeof(realtype));
    for (i=0;i<nslice;i++){
      tau_lev[i] = Physics.tau_lev[i];
    }
    if (Run.rank == 0) {
      cout << "tau_slice: " << nslice << endl;
      for (i=0;i<nslice;i++){
	cout << "tau_slice: " << tau_lev[i]<< endl;
      } 
    }   

    nslvar = 0;
    for (v=0;v<14;v++){
      if (Physics.tau_var[v] == 1) nslvar+=1;
    }

    ini_flag = 0;
  }

  iobuf = (float*) malloc(nslvar*localsize*sizeof(float));
  iosum = (float*) malloc(nslvar*localsize*sizeof(float));

  for (nsl = 0; nsl<nslice; nsl++){
 
    for(v=0;v<nslvar*localsize;v++){
      iobuf[v] = 0.0;
      iosum[v] = 0.0;
    }

    for (k=kbeg; k<=kend; k++){
      for (j=jbeg; j<=jend; j++){
	ind  = j-jbeg + (k-kbeg)*Grid.lsize[1];
        for (i=ibeg; i<=iend; i++){
	  node1 = Grid.node(i,j,k);
	  node2 = Grid.node(i+1,j,k);
	  
	  if( (Grid.Tau[node1] >= tau_lev[nsl]) && (tau_lev[nsl] > Grid.Tau[node2]) ){
	    q1 = (tau_lev[nsl]-Grid.Tau[node2])/(Grid.Tau[node1]-Grid.Tau[node2]);
	    //q1 = (log(tau_lev[nsl])-log(Grid.Tau[node2]))/(log(Grid.Tau[node1])-log(Grid.Tau[node2]));
	    q2 = 1.0-q1;
	    
	    ind1 = ind;
	    
	    if (Physics.tau_var[0] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].d*q1+Grid.U[node2].d*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[1] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.x*q1+Grid.U[node2].M.x*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[2] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.y*q1+Grid.U[node2].M.y*q2); 
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[3] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.z*q1+Grid.U[node2].M.z*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[4] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].e/Grid.U[node1].d*q1+Grid.U[node2].e/Grid.U[node2].d*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[5] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.x*q1+Grid.U[node2].B.x*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[6] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.y*q1+Grid.U[node2].B.y*q2);  
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[7] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.z*q1+Grid.U[node2].B.z*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[8] == 1){
	      iobuf[ind1] = (float) (sqrt(Grid.U[node1].M.sqr())*q1+sqrt(Grid.U[node2].M.sqr())*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[9] == 1){
	      iobuf[ind1] = (float) (sqrt(Grid.U[node1].B.sqr())*q1+sqrt(Grid.U[node2].B.sqr())*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[10] == 1){
	      iobuf[ind1] = (float) (Grid.temp[node1]*q1+Grid.temp[node2]*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[11] == 1){
	      iobuf[ind1] = (float) (Grid.pres[node1]*q1+Grid.pres[node2]*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[12] == 1){
	      iobuf[ind1] = (float) rts->Iout(j,k);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[13] == 1){
	      iobuf[ind1] = float(Grid.coord(i,0)*q1+Grid.coord(i+1,0)*q2)/float(Grid.gxmax[0]);
	    }
	  }
	}
      }
    }

    MPI_Reduce(iobuf,iosum,nslvar*localsize,MPI_FLOAT,MPI_SUM,iroot,
		  XCOL_COMM);

    if (xcol_rank == iroot){

      /*
      sprintf(filename,"%s_%.3f.%06d","tau_slice",tau_lev[nsl],
	      Run.globiter);
      MPI_File_open(YZ_COMM,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
		    MPI_INFO_NULL,&fhandle_mpi);
      
      if(yz_rank == 0) {
	float header[4];            
	header[0] = (float) nslvar;
	header[1] = (float) Grid.gsize[1];
	header[2] = (float) Grid.gsize[2];
	header[3] = (float) Run.time;
	MPI_File_write(fhandle_mpi,header,4,MPI_FLOAT,MPI_STATUS_IGNORE);
      }
      
      offset = 4*sizeof(float);
      
      MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,yz_subarray,"native",
			MPI_INFO_NULL);
      MPI_File_write_all(fhandle_mpi,&(iosum[0]),nslvar*localsize,
			 MPI_FLOAT,MPI_STATUS_IGNORE);
      
      MPI_File_close(&fhandle_mpi);   
      */

      if(Physics.slice[i_sl_collect] == 0) {	
	if(yz_rank == 0) {
	  if(tau_lev[nsl] >= 1e-3)
	    sprintf(filename,"%s%s_%.3f.%06d",Run.path_2D,"tau_slice",
		    tau_lev[nsl],Run.globiter);
	  else
	    sprintf(filename,"%s%s_%.6f.%06d",Run.path_2D,"tau_slice",
		    tau_lev[nsl],Run.globiter);
	  fhandle=fopen(filename,"w");
	  
	  float header[4];            
	  header[0] = (float) nslvar;
	  header[1] = (float) Grid.gsize[1];
	  header[2] = (float) Grid.gsize[2];
	  header[3] = (float) Run.time;
	  fwrite(header,sizeof(float),4,fhandle);
	}
	
	slice_write(Grid,0,&(iosum[0]),localsize,nslvar,1,2,fhandle);
	
	if(yz_rank == 0)
	  fclose(fhandle);
	
      } else {
	if(yz_rank == 0){
	  sprintf(filename,"%s_%.3f.dat","tau_slice",tau_lev[nsl]);
	  fhandle=fopen(filename,"a");
	}
	
	slice_write(Grid,0,&(iosum[0]),localsize,nslvar,1,2,fhandle);
	
	if(yz_rank == 0){ 
	  fclose(fhandle);
	  
	  sprintf(filename,"%s_%.3f.log","tau_slice",tau_lev[nsl]);
	  
      std::fstream fptr;
	  int newfile = 0;
	  fptr.open(filename,std::ios::in);
	  if (!fptr) newfile = 1;
	  fptr.close();
	  
	  fptr.open(filename,std::ios::out|std::ios::app);
	  fptr.precision(10);
	  if (newfile) {      
	    fptr <<  nslvar << ' ' <<  Grid.gsize[1] << ' ' 
		 << Grid.gsize[2] << endl;
	    fptr << Physics.tau_var[0]  << ' ' 
		 << Physics.tau_var[1]  << ' ' 
		 << Physics.tau_var[2]  << ' ' 
		 << Physics.tau_var[3]  << ' ' 
		 << Physics.tau_var[4]  << ' ' 
		 << Physics.tau_var[5]  << ' ' 
		 << Physics.tau_var[6]  << ' ' 
		 << Physics.tau_var[7]  << ' ' 
		 << Physics.tau_var[8]  << ' ' 
		 << Physics.tau_var[9]  << ' ' 
		 << Physics.tau_var[10] << ' ' 
		 << Physics.tau_var[11] << ' '
		 << Physics.tau_var[12] << ' '        
		 << Physics.tau_var[13] << endl;
	  }
	  fptr << Run.globiter << ' ' << Run.time << endl;
	  fptr.close();
	}
      }
    }
  }

  free(iobuf);
  free(iosum);
}

