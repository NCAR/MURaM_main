#ifdef MURAM_FFTW

#include <mpi.h>
#include <fftw3-mpi.h>
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include <iostream>

using std::cout;
using std::endl;

/*
  Routine follows following naming conventions regardless of data layout
  z: vertical direction
  x: fast horizontal direction (second fftw dimension)
  y: slow horizontal direction (first fftw dimension, direction of fftw-mpi decomposition)

  Currently all operations are within the top layer of processors, i.e. it is sufficient to call
  the routine only in the boundary layer. 
*/

void x_gather(const GridData&, const int, double*, double*, const int, const int, const int, const MPI_Comm);
void x_scatter_gh(const GridData&, const int, double*, double*, const int, const int, const int, const MPI_Comm);
void y_exchange(const GridData&, double*, const int, const int);

void potential_ext_fftw3(const GridData& Grid, double *bz0, double *b_ext){

  static int ini_flag = 1;

  static fftw_complex *kernel_fft[8];
  static double *bz0_x, *b_ext_x;
  static fftw_complex *data_in, *data_out, *bz_fft;
  static fftw_plan plan_fwd,plan_bwd;

  //double time,stime,fft_time,ttime;
  
  ptrdiff_t N0, N1, alloc_local, local_n0, local_0_start, i, j, k;

  // ==============================================================
  // ******** Change this for different dimensional layout ********
  int hdir_x = 1; // fast dimension
  int hdir_y = 2; // slow dimension, direction of fftw-mpi decomposition

  MPI_Comm x_comm = YCOL_COMM;
  MPI_Comm y_comm = ZCOL_COMM;

  int x_rank = ycol_rank;
  int y_rank = zcol_rank;
  //===============================================================

  int gx = Grid.gsize[hdir_x];
  int gy = Grid.gsize[hdir_y];

  int ly = Grid.lsize[hdir_y];

  int gh_x = Grid.ghosts[hdir_x];
  int gh_y = Grid.ghosts[hdir_y];
  
  int gx_gh = gx+2*gh_x;
  int ly_gh = ly+2*gh_y;

  N0 = gy;
  N1 = gx;

  int buf_sz    = ly*gx;
  int buf_sz_gh = ly_gh*gx_gh;

  if( (Grid.NDIM == 3) and (Grid.dx[hdir_x] != Grid.dx[hdir_y]) ){
    cout << "dx=dz required for potential field bnd: " 
	 << Grid.dx[hdir_x] << ' ' << Grid.dx[hdir_y] << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
    
  double norm,gnorm;
    
  char fnames[8][128];
  sprintf(fnames[0],"%s","PSF-kernel-z-1.dat");
  sprintf(fnames[1],"%s","PSF-kernel-z-2.dat");
  sprintf(fnames[2],"%s","PSF-kernel-x-0.dat");
  sprintf(fnames[3],"%s","PSF-kernel-x-1.dat");
  sprintf(fnames[4],"%s","PSF-kernel-x-2.dat");
  sprintf(fnames[5],"%s","PSF-kernel-y-0.dat");
  sprintf(fnames[6],"%s","PSF-kernel-y-1.dat");
  sprintf(fnames[7],"%s","PSF-kernel-y-2.dat");
  
  if(ini_flag){

    fftw_mpi_init();

    if(x_rank == 0){

      // ======== fftw mpi and data layout ========
      
      alloc_local = fftw_mpi_local_size_2d(N0, N1, y_comm,&local_n0, &local_0_start);

      if(alloc_local != buf_sz){
        cout << "FFTW buffer inconsistent with MURaM decomposition:" << alloc_local << ' ' << buf_sz << endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(local_n0 != ly){
	cout << "FFTW decomposition inconsistent with MURaM decomposition:" << local_n0 << ' ' << ly << endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      
      for(i=0;i<8;i++)
	kernel_fft[i] = fftw_alloc_complex(alloc_local);

      bz0_x   = new double[buf_sz];
      b_ext_x = new double[8*buf_sz_gh];
      for(i=0;i<8*buf_sz_gh;i++) b_ext_x[i] = 0.0;
      
      data_in  = fftw_alloc_complex(buf_sz);
      data_out = fftw_alloc_complex(buf_sz);
      bz_fft   = fftw_alloc_complex(buf_sz);
      
      plan_fwd = fftw_mpi_plan_dft_2d(N0, N1, data_in, data_out, y_comm, FFTW_FORWARD,  FFTW_MEASURE);
      plan_bwd = fftw_mpi_plan_dft_2d(N0, N1, data_in, data_out, y_comm, FFTW_BACKWARD, FFTW_MEASURE);
      
      //  ======== Create 1D subarray for parallel IO ======== 
      MPI_Datatype io_subarray;
      int sz[2],sub_sz[2],strt[2];
      
      sz[0]     = N1;
      sz[1]     = N0;
      sub_sz[0] = N1;
      sub_sz[1] = local_n0;
      strt[0]   = 0;
      strt[1]   = local_0_start;
      MPI_Type_create_subarray(2,sz,sub_sz,strt,MPI_ORDER_FORTRAN,MPI_FLOAT,&io_subarray);
      MPI_Type_commit(&io_subarray);
      
      //  ======== read kernels and compute fft ========
      MPI_File fhandle_mpi;
      float header[3];
      float *iobuf = new float[buf_sz];
      
      for(k=0;k<8;k++){
	MPI_File_open(y_comm,fnames[k],MPI_MODE_RDONLY,MPI_INFO_NULL,&fhandle_mpi);
	
	if(y_rank == 0) {
	  MPI_File_read(fhandle_mpi,header,3,MPI_FLOAT,MPI_STATUS_IGNORE);
	  if( (int(header[0]) != Grid.gsize[hdir_x]) or 
	      (int(header[1]) != Grid.gsize[hdir_y]) or
	      (header[2] != float(Grid.dx[0]/Grid.dx[hdir_x])) ){
	    cout << "inconsistent header: " 
		 << int(header[0]) << " <-> " << Grid.gsize[hdir_x] << " | " 
		 << int(header[1]) << " <-> " << Grid.gsize[hdir_y] << " | "
		 << header[2] << " <-> " << float(Grid.dx[1]/Grid.dx[0]) << endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}
	
	int offset = 3*sizeof(float);
	
	MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,io_subarray,"native",MPI_INFO_NULL);
	MPI_File_read_all(fhandle_mpi,iobuf,buf_sz,MPI_FLOAT,MPI_STATUS_IGNORE);
	
	MPI_File_close(&fhandle_mpi);
	
	norm=0.0;
	for (i=0;i<buf_sz;i++) norm+= iobuf[i];
	MPI_Reduce(&norm,&gnorm,1,MPI_DOUBLE,MPI_SUM,0,y_comm);
	if(y_rank == 0) 
	  cout << "Kernel [" << k << "] norm: " << gnorm << endl;
	
	for (i=0;i<buf_sz;i++){
	  data_in[i][0] = iobuf[i];
	  data_in[i][1] = 0.0;
	}
	
	fftw_execute(plan_fwd);
	
	for (i=0;i<buf_sz;i++){
	  kernel_fft[k][i][0] = data_out[i][0]/double(N0*N1);
	  kernel_fft[k][i][1] = data_out[i][1]/double(N0*N1);
	}
      }
      
      // ======= clean up temporary stuff ========
      
      delete[] iobuf;
      MPI_Type_free(&io_subarray);
      
    }
   
    ini_flag = 0;
  }

  //stime = MPI_Wtime();
  //fft_time = 0.0;
 
  x_gather(Grid,0,bz0,bz0_x,hdir_x,hdir_y,x_rank,x_comm);

  if(x_rank == 0){
    
    // ======== forward fft of bz ========
    for (i=0;i<buf_sz;i++){
      
      data_in[i][0] = bz0_x[i];
      data_in[i][1] = 0.0;
    }

    //time = MPI_Wtime();
    fftw_execute(plan_fwd);
    //fft_time += MPI_Wtime()-time; 

    for (i=0;i<buf_sz;i++){
      bz_fft[i][0] = data_out[i][0];
      bz_fft[i][1] = data_out[i][1];
    }

    //  ======== convolution with kernels ========
    
    for(k=0;k<8;k++){
      for (i=0;i<buf_sz;i++){
	data_in[i][0] = bz_fft[i][0]*kernel_fft[k][i][0]-bz_fft[i][1]*kernel_fft[k][i][1];
	data_in[i][1] = bz_fft[i][0]*kernel_fft[k][i][1]+bz_fft[i][1]*kernel_fft[k][i][0];
      }

      //time = MPI_Wtime();
      fftw_execute(plan_bwd);
      //fft_time += MPI_Wtime()-time; 

      for(j=0;j<ly;j++)
	for(i=0;i<gx;i++)
	  b_ext_x[k*buf_sz_gh+(j+gh_y)*gx_gh+i+gh_x]=data_out[j*gx+i][0];
    }

    // ======= fill periodic boundary values in x-direction ========

    for(k=0;k<8;k++){
      for(j=0;j<ly;j++){
	for(i=0;i<gh_x;i++){
	  b_ext_x[k*buf_sz_gh+(j+gh_y)*gx_gh+i] = b_ext_x[k*buf_sz_gh+(j+gh_y)*gx_gh+i+gx];
	  b_ext_x[k*buf_sz_gh+(j+gh_y)*gx_gh+i+gx+gh_x] = b_ext_x[k*buf_sz_gh+(j+gh_y)*gx_gh+i+gh_x];
	}
      }
    }
    
  }

  // ======= scatter in x-direction including ghostcells ========
  
  x_scatter_gh(Grid,0,b_ext_x,b_ext,hdir_x,hdir_y,x_rank,x_comm);    

  // ======= fill ghost cells in y-direction ======== 
  
  y_exchange(Grid,b_ext,hdir_x,hdir_y);

  //ttime = MPI_Wtime()-stime;

  //if(x_rank*y_rank == 0) cout << "Potential time: " << fft_time << ' ' << ttime << ' ' << fft_time/ttime << endl;
 
}

// =========================================================================================
void x_gather(const GridData& Grid, const int root, double* bz0, double* bz0_x,
	      const int hdir_x, const int hdir_y, const int x_rank, const MPI_Comm x_comm){

  static int ini_flag = 1;
  
  int np,i,j,iloc,iglo;
  int Gb[2], Ge[2], bounds[4];
  const int nprocs = Grid.procs[hdir_x];
  
  static int* proc_bounds;
  static int* recvcounts;
  static int* offsets;
  static int recvbufsize;
 
  double* recvbuf=NULL;

  if (ini_flag == 1){
    proc_bounds = new int[4*nprocs];
    recvcounts  = new int[nprocs];
    offsets     = new int[nprocs];
  
    bounds[0]=Grid.beg[hdir_x];
    bounds[1]=Grid.end[hdir_x];
    bounds[2]=Grid.beg[hdir_y];
    bounds[3]=Grid.end[hdir_y];

    MPI_Allgather(bounds,4,MPI_INT,proc_bounds,4,MPI_INT,x_comm);
    
    for (np=0;np<nprocs;np++){
      Gb[0]=proc_bounds[4*np+0];
      Ge[0]=proc_bounds[4*np+1];
      Gb[1]=proc_bounds[4*np+2];
      Ge[1]=proc_bounds[4*np+3];
      recvcounts[np]=(Ge[0]-Gb[0]+1)*(Ge[1]-Gb[1]+1); 
    }
    offsets[0]=0;
    for (np=1;np<nprocs;np++) offsets[np]=offsets[np-1]+recvcounts[np-1];
    recvbufsize=0;
    for (np=0;np<nprocs;np++) recvbufsize+=recvcounts[np];

    ini_flag = 0;
  }


  if(x_rank == root){
    recvbuf = new double[recvbufsize];
    for(i=0;i<recvbufsize;i++) recvbuf[i] = 0.0;
  }

  MPI_Gatherv(bz0,recvcounts[x_rank],MPI_DOUBLE,recvbuf,recvcounts,offsets,
	      MPI_DOUBLE,root,x_comm);

  if (x_rank == root){   
    for(np=0; np<nprocs; np++){
      Gb[0]=proc_bounds[4*np+0];
      Ge[0]=proc_bounds[4*np+1];
      Gb[1]=proc_bounds[4*np+2];
      Ge[1]=proc_bounds[4*np+3];

      for (j=0; j< Ge[1]-Gb[1]+1;j++){ 
	for (i=0; i< Ge[0]-Gb[0]+1;i++){
          iloc = j*(Ge[0]-Gb[0]+1)+i+offsets[np];
          iglo = j*Grid.gsize[hdir_x]+i+Gb[0]-Grid.gbeg[hdir_x];
          bz0_x[iglo] = recvbuf[iloc];
        }
      }
    }
  }

  if(x_rank == root) delete[] recvbuf;

}
// ============================================================================================
void x_scatter_gh(const GridData& Grid, const int root, double* b_ext_x, double* b_ext,
		  const int hdir_x, const int hdir_y, const int x_rank, const MPI_Comm x_comm){
  
  static int ini_flag = 1;
  int np,i,j,k,iloc,iglo;
  int Gb[2], Ge[2], bounds[4];
  const int nprocs = Grid.procs[hdir_x];
  
  static int* proc_bounds;
  static int* sendcounts;
  static int* offsets;
  static int sendbufsize;

  double* sendbuf=NULL;

  if (ini_flag == 1){
    proc_bounds = new int[4*nprocs];
    sendcounts  = new int[nprocs];
    offsets     = new int[nprocs];
  
    bounds[0]=Grid.beg[hdir_x]-Grid.ghosts[hdir_x];
    bounds[1]=Grid.end[hdir_x]+Grid.ghosts[hdir_x];
    bounds[2]=Grid.beg[hdir_y]-Grid.ghosts[hdir_y];
    bounds[3]=Grid.end[hdir_y]+Grid.ghosts[hdir_y];

    MPI_Allgather(bounds,4,MPI_INT,proc_bounds,4,MPI_INT,x_comm);
    
    for (np=0;np<nprocs;np++){
      Gb[0]=proc_bounds[4*np+0];
      Ge[0]=proc_bounds[4*np+1];
      Gb[1]=proc_bounds[4*np+2];
      Ge[1]=proc_bounds[4*np+3];
      sendcounts[np]=8*(Ge[0]-Gb[0]+1)*(Ge[1]-Gb[1]+1);
    }
    offsets[0]=0;
    for (np=1;np<nprocs;np++) offsets[np]=offsets[np-1]+sendcounts[np-1];
    sendbufsize=0;
    for (np=0;np<nprocs;np++) sendbufsize+=sendcounts[np];

    ini_flag = 0;
  }

   if(x_rank == root){
     sendbuf = new double[sendbufsize];
    for(i=0;i<sendbufsize;i++) sendbuf[i] = 0.0;
  }
  
  if (x_rank == root){
    for(np=0; np<nprocs; np++){
      Gb[0]=proc_bounds[4*np+0];
      Ge[0]=proc_bounds[4*np+1];
      Gb[1]=proc_bounds[4*np+2];
      Ge[1]=proc_bounds[4*np+3];

      for(k=0;k<8;k++){
	for (j=0; j< Ge[1]-Gb[1]+1;j++){ 
	  for (i=0; i< Ge[0]-Gb[0]+1;i++){
	    iloc = j*(Ge[0]-Gb[0]+1)+i+offsets[np];
	    iglo = j*(Grid.gsize[hdir_x]+2*Grid.ghosts[hdir_x])+i+Gb[0]-(Grid.gbeg[hdir_x]-Grid.ghosts[hdir_x]);
	    iloc+= k*(Ge[0]-Gb[0]+1)*(Ge[1]-Gb[1]+1);
	    iglo+= k*(Grid.gend[hdir_x]-Grid.gbeg[hdir_x]+1+2*Grid.ghosts[hdir_x])*(Ge[1]-Gb[1]+1);
	    sendbuf[iloc] = b_ext_x[iglo]; 
	  }
	}
      }
    }
  }

  MPI_Scatterv(sendbuf,sendcounts,offsets,MPI_DOUBLE,b_ext,sendcounts[x_rank],MPI_DOUBLE,root,x_comm); 

  if(x_rank == root) delete[] sendbuf;
}
// ======================================================================================================================================
void y_exchange(const GridData& Grid, double* b_ext, const int hdir_x, const int hdir_y){

  static int ini_flag=1;
 
  register int i,j,k,buf;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int* leftr =Grid.leftr;
  const int* rightr=Grid.rightr;
  
  static int x_sz,y_sz,y_gh,y_beg,y_end,buf_sz;
  static double *sndbuf_l , *recbuf_l, *sndbuf_r, *recbuf_r;

  if(ini_flag){  
    y_gh  = Grid.ghosts[hdir_y];
    y_beg = Grid.ghosts[hdir_y];
    y_end = Grid.lsize[hdir_y]+Grid.ghosts[hdir_y];

    x_sz  = Grid.lsize[hdir_x]+2*Grid.ghosts[hdir_x];
    y_sz  = Grid.lsize[hdir_y]+2*Grid.ghosts[hdir_y];
    
    buf_sz = x_sz*y_gh;
    
    sndbuf_l = new double[8*buf_sz];
    recbuf_l = new double[8*buf_sz];
    sndbuf_r = new double[8*buf_sz];
    recbuf_r = new double[8*buf_sz];

    ini_flag=0;
  }
  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=y_beg;j<y_beg+y_gh;j++){
      for(i=0;i<x_sz;i++){
	sndbuf_l[buf++] = b_ext[k*x_sz*y_sz+j*x_sz+i];
      }  
    }
  }

  MPI_Irecv(recbuf_r,8*buf_sz,MPI_DOUBLE,rightr[hdir_y],0,MPI_COMM_WORLD,rr); 
  MPI_Isend(sndbuf_l,8*buf_sz,MPI_DOUBLE,leftr[hdir_y], 0,MPI_COMM_WORLD,rr+1);

  buf = 0;
  for(k=0;k<8;k++){
    for(j=y_end-y_gh;j<y_end;j++){
      for(i=0;i<x_sz;i++){
	sndbuf_r[buf++] = b_ext[k*x_sz*y_sz+j*x_sz+i];
      }  
    }
  }

  MPI_Irecv(recbuf_l,8*buf_sz,MPI_DOUBLE,leftr[hdir_y], 1,MPI_COMM_WORLD,rl);   
  MPI_Isend(sndbuf_r,8*buf_sz,MPI_DOUBLE,rightr[hdir_y],1,MPI_COMM_WORLD,rl+1);

  MPI_Waitall(2,rr,st);
  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=y_end;j<y_end+y_gh;j++){
      for(i=0;i<x_sz;i++){
	b_ext[k*x_sz*y_sz+j*x_sz+i] = recbuf_r[buf++];
      }  
    }
  }

  MPI_Waitall(2,rl,st);

  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=0;j<y_gh;j++){
      for(i=0;i<x_sz;i++){
	b_ext[k*x_sz*y_sz+j*x_sz+i] = recbuf_l[buf++];
      }  
    }
  } 
}

#endif
