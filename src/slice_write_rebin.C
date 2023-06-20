#include <mpi.h>
#include "grid.H"
#include "run.H"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "comm_split.H"

using std::cout;
using std::endl;
//========================================================================
void slice_write_rebin(const GridData& Grid,const int iroot,float* vloc,
		       const int nloc,const int nvar,const int n0,
		       const int n1,const int sm_x,const int sm_y,
		       FILE* fhandle){
  
  MPI_Comm comm=MPI_COMM_NULL;
  int rank=-1;

  if ( (n0 == 2) and (n1 == 0)){
    comm = XZ_COMM; rank = xz_rank;
  } else if ( (n0 == 1) and (n1 == 0)){
    comm = XY_COMM; rank = xy_rank;
  } else if ( (n0 == 1) and (n1 == 2)){
    comm = YZ_COMM; rank = yz_rank;
  } else {
    cout << "slice_write: mode unknown " << n0 << ' ' << n1 << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  
  int np,i,k,iloc,iglo,v;
  int Gb[2], Ge[2], bounds[4];
  int isz,ksz,ioff,koff;

  int nprocs = Grid.procs[n0]*Grid.procs[n1];
  int slsize = Grid.gsize[n0]*Grid.gsize[n1];
  int slsize_sm = (Grid.gsize[n0]/sm_x)*(Grid.gsize[n1]/sm_y);
  
  int* proc_bounds;
  int* recvcounts;
  int* offsets;

  float* recvbuf=NULL;
  float* iobuf=NULL;
  float* iobuf_sm=NULL;

  int recvbufsize;

  proc_bounds = (int*) malloc(4*nprocs*sizeof(int));
  recvcounts  = (int*) malloc(nprocs*sizeof(int));
  offsets     = (int*) malloc(nprocs*sizeof(int));

  bounds[0]=Grid.beg[n0];
  bounds[1]=Grid.end[n0];
  bounds[2]=Grid.beg[n1];
  bounds[3]=Grid.end[n1];  
  
  MPI_Allgather(bounds,4,MPI_INT,proc_bounds,4,MPI_INT,comm);
    
  for (np=0;np<nprocs;np++){
    Gb[0]=proc_bounds[4*np+0];
    Ge[0]=proc_bounds[4*np+1];
    Gb[1]=proc_bounds[4*np+2];
    Ge[1]=proc_bounds[4*np+3];  
    recvcounts[np]=(Ge[0]-Gb[0]+1)*(Ge[1]-Gb[1]+1)*nvar; 
  }
  offsets[0]=0;
  for (np=1;np<nprocs;np++){  
    offsets[np]=offsets[np-1]+recvcounts[np-1];
  }

  recvbufsize=0;
  for (np=0;np<nprocs;np++){   
    recvbufsize+=recvcounts[np];
  }
    
  if(rank == iroot){
    if (slsize*nvar != recvbufsize){
      cout << "slice_write: nglo != recvbufsize " << rank << ' ' 
	   << slsize*nvar << ' ' << recvbufsize << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    recvbuf  = (float*) malloc(recvbufsize*sizeof(float));
    iobuf    = (float*) malloc(slsize*sizeof(float));
  }
  
  if(nloc*nvar != recvcounts[rank]){
    cout << "slice_write: nloc != recvcounts " << rank << ' ' 
         << nloc << ' ' << recvcounts[rank] << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  
  MPI_Gatherv(vloc,nloc*nvar,MPI_FLOAT,recvbuf,recvcounts,offsets,
              MPI_FLOAT,iroot,comm);
  
  if (rank == iroot){
    for (v=0;v<nvar;v++){
      for(np=0; np<nprocs; np++){
	Gb[0]=proc_bounds[4*np+0];
	Ge[0]=proc_bounds[4*np+1];
	Gb[1]=proc_bounds[4*np+2];
	Ge[1]=proc_bounds[4*np+3];
	
	isz  = Ge[0]-Gb[0]+1;
	ksz  = Ge[1]-Gb[1]+1;
	ioff = Gb[0]-Grid.gbeg[n0];
	koff = Gb[1]-Grid.gbeg[n1];

	for (k=0;k<ksz;k++){ 
	  for (i=0;i<isz;i++){
	    iloc = i + isz*(k+v*ksz) + offsets[np];
	    iglo = (i+ioff) + (k+koff)*Grid.gsize[n0];
	    iobuf[iglo] = recvbuf[iloc];
	  }
	}
      }

      if (sm_x*sm_y > 1) {
	iobuf_sm = (float*) malloc(slsize_sm*sizeof(float));
	for (k=0;k<Grid.gsize[n1]/sm_y;k++){
	  for (i=0;i<Grid.gsize[n0]/sm_x;i++){
	    iloc= i + k*Grid.gsize[n0]/sm_x;
	    iobuf_sm[iloc] = 0.0;
	    int is,ks;
	    for (ks=0;ks<sm_y;ks++){
	      for (is=0;is<sm_x;is++){
		iglo = sm_x*i+is + (sm_y*k+ks)*Grid.gsize[n0];
		iobuf_sm[iloc] += iobuf[iglo]/float(sm_x*sm_y);
	      }
	    }
	  }
	}
	fwrite(iobuf_sm,sizeof(float),slsize_sm,fhandle);
	free(iobuf_sm);
      } else {
	fwrite(iobuf,sizeof(float),slsize,fhandle);
      }
    }

    free(recvbuf);
    free(iobuf);
  }

  free(proc_bounds);
  free(recvcounts);
  free(offsets);

}
