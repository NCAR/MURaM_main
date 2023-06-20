#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include "grid.H"
#include "run.H"

MPI_Comm cart_comm;
MPI_Comm XY_COMM,XZ_COMM,YZ_COMM;
MPI_Comm XCOL_COMM,YCOL_COMM,ZCOL_COMM;
MPI_Datatype xy_subarray,xz_subarray,yz_subarray;
int xy_rank,xz_rank,yz_rank,xcol_rank,ycol_rank,zcol_rank;

void comm_split_finalize(){
  MPI_Comm_free(&XY_COMM);
  MPI_Comm_free(&XZ_COMM);
  MPI_Comm_free(&YZ_COMM);
  MPI_Comm_free(&XCOL_COMM);
  MPI_Comm_free(&YCOL_COMM);
  MPI_Comm_free(&ZCOL_COMM);
  MPI_Type_free(&xy_subarray);
  MPI_Type_free(&xz_subarray);
  MPI_Type_free(&yz_subarray); 
}

void comm_split_init(RunData&  Run, const GridData& Grid) {
  
  // slice subcommunicators
  int lrank[3];
  int x_color,y_color,z_color;
  int ndim=3;  
        
  MPI_Cart_coords(cart_comm,Run.rank,ndim,lrank);
 
  MPI_Comm_split(MPI_COMM_WORLD,lrank[2],Run.rank,&XY_COMM);    
  MPI_Comm_split(MPI_COMM_WORLD,lrank[1],Run.rank,&XZ_COMM);
  MPI_Comm_split(MPI_COMM_WORLD,lrank[0],Run.rank,&YZ_COMM);

  x_color = lrank[1] + lrank[2]*Grid.procs[1];
  y_color = lrank[0] + lrank[2]*Grid.procs[0];
  z_color = lrank[0] + lrank[1]*Grid.procs[0];

  // column subcommunicators

  MPI_Comm_split(MPI_COMM_WORLD,x_color,Run.rank,&XCOL_COMM);
  MPI_Comm_split(MPI_COMM_WORLD,y_color,Run.rank,&YCOL_COMM);
  MPI_Comm_split(MPI_COMM_WORLD,z_color,Run.rank,&ZCOL_COMM);

  MPI_Comm_rank(XY_COMM,&xy_rank);
  MPI_Comm_rank(XZ_COMM,&xz_rank);
  MPI_Comm_rank(YZ_COMM,&yz_rank);
  MPI_Comm_rank(XCOL_COMM,&xcol_rank);
  MPI_Comm_rank(YCOL_COMM,&ycol_rank);
  MPI_Comm_rank(ZCOL_COMM,&zcol_rank);

  Run.zrank = xcol_rank;

  // 2d subarrays for mpi io
  int array_of_sizes[2];
  int array_of_subsizes[2];
  int array_of_starts[2];

  // xy, transpose array
  array_of_sizes[0]=Grid.gsize[1];
  array_of_sizes[1]=Grid.gsize[0];
  array_of_subsizes[0]=Grid.lsize[1];
  array_of_subsizes[1]=Grid.lsize[0];
  array_of_starts[0]=Grid.beg[1]-Grid.gbeg[1];
  array_of_starts[1]=Grid.beg[0]-Grid.gbeg[0];
  
  MPI_Type_create_subarray(2,array_of_sizes,array_of_subsizes,
			   array_of_starts,MPI_ORDER_FORTRAN,
			   MPI_FLOAT,&xy_subarray);
    
  MPI_Type_commit(&xy_subarray);

  // xz, transpose array
  array_of_sizes[0]=Grid.gsize[2];
  array_of_sizes[1]=Grid.gsize[0];
  array_of_subsizes[0]=Grid.lsize[2];
  array_of_subsizes[1]=Grid.lsize[0];
  array_of_starts[0]=Grid.beg[2]-Grid.gbeg[2];
  array_of_starts[1]=Grid.beg[0]-Grid.gbeg[0];
  
  MPI_Type_create_subarray(2,array_of_sizes,array_of_subsizes,
			   array_of_starts,MPI_ORDER_FORTRAN,
			   MPI_FLOAT,&xz_subarray);
    
  MPI_Type_commit(&xz_subarray);

  //yz
  array_of_sizes[0]=Grid.gsize[1];
  array_of_sizes[1]=Grid.gsize[2];
  array_of_subsizes[0]=Grid.lsize[1];
  array_of_subsizes[1]=Grid.lsize[2];
  array_of_starts[0]=Grid.beg[1]-Grid.gbeg[1];
  array_of_starts[1]=Grid.beg[2]-Grid.gbeg[2];
  
  MPI_Type_create_subarray(2,array_of_sizes,array_of_subsizes,
			   array_of_starts,MPI_ORDER_FORTRAN,
			   MPI_FLOAT,&yz_subarray);
    
  MPI_Type_commit(&yz_subarray);
}


