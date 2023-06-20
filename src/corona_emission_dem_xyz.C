#include <mpi.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "io.H"
#include "comm_split.H"
#include <iostream>

using std::cout;
using std::endl;
using std::max;
using std::min;

#define OUTER_LOOP(G,i,j,d1,d2) \
  for((j)=(G)[(d2)][0];(j)<=(G)[(d2)][1];(j)++) \
  for((i)=(G)[(d1)][0];(i)<=(G)[(d1)][1];(i)++)

//extern void slice_write(const GridData&,const int,float*,int,int,const int,
//			const int,FILE*);

extern void slice_write_rebin(const GridData&,const int,float*,const int,const int,
                              const int,const int,const int,const int,FILE*);

inline int imin(int a, int b) { return a < b ? a : b; }
inline int imax(int a, int b) { return a > b ? a : b; }
//======================================================================
void corona_emission_dem_xyz(const RunData&  Run, const GridData& Grid, 
		               const PhysicsData& Physics) {

  static int ini_flag = 1;
  
  double clock1;
  const int iroot = 0;
  
  const double lgTmin = 4.5;
  const double dellgT = 0.1;
  
  const int nout = 4;
  
  const int rebin[3] = {3,1,1};
  
  register int i, j, k, node1, node2, ind, v, v1, v2, ind1, d, d1, d2, d3;
  
  int str, offset;
  
  int stride[3];
  int kbeg, kend, jbeg, jend, ibeg, iend, glsize;
  const int bufsize = Grid.bufsize;

  cState* U  = (cState*) Grid.U;
 
  stride[0] = Grid.stride[0];
  stride[1] = Grid.stride[1];
  stride[2] = Grid.stride[2];
  
  int bounds[3][2];
  
  for(d=0;d<3;d++){
    bounds[d][0] = Grid.lbeg[d];
    bounds[d][1] = Grid.lend[d];
  }
  
  const int loop_order[3][3] = {{ 1, 2, 0 },{ 0, 1, 2 },{ 2, 1, 0 }};
  
  double* los_sum_loc;
  double* los_sum;
  float*  io_buf;

  char filename[128];

  double r14a,r14b,r14,t6a,t6b,va,vb,vv,rfac,dx,p_a,p_b,s,t1,t2,tmin,tmax,lgTmax;

  double lgT0,del0;

  tmax = 0;
  LOCAL_LOOP(Grid,i,j,k){
    tmax = max(tmax,Grid.temp[Grid.node(i,j,k)]);
  }
  
  MPI_Allreduce(&tmax,&lgTmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  
  lgTmax = log10(lgTmax);

  if (lgTmax < lgTmin)
    lgTmax = lgTmin;

  lgT0=lgTmin;
  del0=1.0/dellgT;
  
  int nslvar = (int) ((lgTmax-lgTmin)/dellgT + 1);
   
  double* tlev = new double[nslvar+1];
  
  for(i=0;i<nslvar+1;i++)
    tlev[i] = lgTmin+double(i)*dellgT;

  if(ini_flag){
    if(Run.rank==0){
      cout << "DEM: Use log(T), log(rho)" << endl;
    }
    ini_flag = 0;
  }

  if(Run.rank==0)
    cout << "DEM: lgTmin = " << lgTmin <<  " lgTmax = " << lgTmax << ' ' << nslvar  << endl;

  FILE* fhandle=NULL;

  for(d=0;d<Grid.NDIM;d++){
    d1=loop_order[d][0];
    d2=loop_order[d][1];
    d3=loop_order[d][2];

    int localsize  = Grid.lsize[d2]*Grid.lsize[d3];
 
    los_sum_loc = new double[nout*nslvar*localsize];
    los_sum     = new double[nout*nslvar*localsize];
    io_buf      = new float[nout*nslvar*localsize];

#pragma acc parallel loop copy(los_sum_loc[:nout*nslvar*localsize])
    for(v=0;v<nout*nslvar*localsize;v++){
      los_sum_loc[v] = 0.0;
      los_sum[v]     = 0.0;
      io_buf[v]      = 0.0;
    }

    str = stride[d1];
    dx  = Grid.dx[d1];

    clock1=MPI_Wtime();

    kbeg = bounds[d3][0];
    kend = bounds[d3][1]; 
    jbeg = bounds[d2][0];
    jend = bounds[d2][1];
    ibeg = bounds[d1][0];
    iend = bounds[d1][1];
    glsize = Grid.lsize[d2];

#pragma acc parallel loop collapse(2) gang                          \
  copy(los_sum_loc[:nout*nslvar*localsize])                         \
    present(Grid[:1], Grid.temp[:bufsize], U[:bufsize]            \
      ) copyin(tlev[:nslvar+1], dx, str,stride[:3], localsize,d1,d2,d3)  default(present)\
  private(ind, offset, k, j,ind1,node1, node2, t6a,\
         t6b, r14a, r14b, va, vb, tmin, tmax, \
              ind1, v1, v2,t1, t2, rfac, s, p_a, p_b, r14, vv, i, v)
    for(int k = kbeg; k <= kend;k++){
    for(int j = jbeg; j <= jend;j++){
      offset = j*stride[d2]+k*stride[d3];
      ind = j-jbeg +(k-kbeg)*glsize; 
        #pragma acc loop seq
        for(i=ibeg;i<=iend;i++){

        node1 = offset+i*str;
	node2 = offset+(i+1)*str;

	t6a  = log10(Grid.temp[node1]);
	t6b  = log10(Grid.temp[node2]);
	r14a = log(U[node1].d*1e14);
	r14b = log(U[node2].d*1e14);
	
	if(d == 0)
	{
	va   = U[node1].M.y;
	vb   = U[node2].M.y;
	}
	else if(d == 1)
	{
	va   = U[node1].M.x;
	vb   = U[node2].M.x;
	}
	else
	{
	va   = U[node1].M.z;
	vb   = U[node2].M.z;
	}
	
	ind1 = ind;

	tmin = min(t6a,t6b);
	tmax = max(t6a,t6b);

	v1 = (int) ( (tmin-lgT0)*del0 );
	v2 = (int) ( (tmax-lgT0)*del0 );

	v1 = imax(0,v1);
	v2 = imin(nslvar-1,v2);

	ind1 += v1*localsize;
        #pragma acc loop seq
	for(v=v1;v<=v2;v++){
  
	  t1=max(tmin,tlev[v]);
	  t2=min(tmax,tlev[v+1]);
		 
	  if(t2-t1 > 1e-6){
	    rfac = (t2-t1)/(tmax-tmin);
	    s    = 0.5*(t1+t2);
	    p_a  = (t6b-s)/(t6b-t6a);
	    p_b  = (s-t6a)/(t6b-t6a);
	  } else if (t2-t1 >= 0.0) {
	    rfac = 1.0;
	    p_a  = 0.5;
	    p_b  = 0.5;
	  } else {
	    rfac = 0.0;
	    p_a  = 0.5;
	    p_b  = 0.5;
	  }

	  r14  = exp(p_a*r14a+p_b*r14b);
	  r14  = r14*r14;
	  vv   = p_a*va+p_b*vb;

	  if(r14 > 1e10) rfac = 0.0;
	    
	  los_sum_loc[ind1]                    += rfac*dx;
	  los_sum_loc[ind1+1*nslvar*localsize] += rfac*dx*r14;
	  los_sum_loc[ind1+2*nslvar*localsize] += rfac*dx*r14*vv;
	  los_sum_loc[ind1+3*nslvar*localsize] += rfac*dx*r14*vv*vv;
	  
	  ind1 += localsize;
	}
	
      }
    } //OUTER_LOOP
    }

    if(d1 == 0){
      MPI_Reduce(los_sum_loc,los_sum,nout*nslvar*localsize,MPI_DOUBLE,MPI_SUM,iroot,
		 XCOL_COMM);

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+nslvar*localsize]; 
	if( rfac > 0.0 ){
	  los_sum[i+2*nslvar*localsize] /= rfac; 
	  los_sum[i+3*nslvar*localsize] /= rfac;
	}
      }

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+3*nslvar*localsize]-los_sum[i+2*nslvar*localsize]*los_sum[i+2*nslvar*localsize];
	los_sum[i+3*nslvar*localsize] = sqrt(max(rfac,0.0));
      }

      for(v=0;v<nout*nslvar*localsize;v++)
	io_buf[v] = (float) los_sum[v];
      
      if (xcol_rank == iroot){
	for (v=0;v<nout;v++){
	  if(yz_rank == 0) {
	    if(v == 0)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_fil_x",Run.globiter);
	    if(v == 1)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_dem_x",Run.globiter);
	    if(v == 2)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vlos_x",Run.globiter);
	    if(v == 3)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vrms_x",Run.globiter);
	    
	    fhandle=fopen(filename,"w");

	    float header[6];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    header[3] = (float) Run.time;
	    header[4] = (float) lgTmin;
	    header[5] = (float) dellgT;
	    fwrite(header,sizeof(float),6,fhandle);
	  }
	
	  //slice_write(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,fhandle);
	  slice_write_rebin(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,rebin[d2],rebin[d3],fhandle);
	  
	  if(yz_rank == 0)
	    fclose(fhandle);
	}
      }
    }
    
    if(d1 == 1){    
      MPI_Reduce(los_sum_loc,los_sum,nout*nslvar*localsize,MPI_DOUBLE,MPI_SUM,iroot,
		 YCOL_COMM);
      
      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+nslvar*localsize]; 
	if( rfac > 0.0 ){
	  los_sum[i+2*nslvar*localsize] /= rfac; 
	  los_sum[i+3*nslvar*localsize] /= rfac;
	}
      }

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+3*nslvar*localsize]-los_sum[i+2*nslvar*localsize]*los_sum[i+2*nslvar*localsize];
	los_sum[i+3*nslvar*localsize] = sqrt(max(rfac,0.0));
      }

      for(v=0;v<nout*nslvar*localsize;v++)
	io_buf[v] = (float) los_sum[v];
      
      if (ycol_rank == iroot){
	for (v=0;v<nout;v++){
	  if(xz_rank == 0) {
	    if(v == 0)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_fil_y",Run.globiter);
	    if(v == 1)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_dem_y",Run.globiter);
	    if(v == 2)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vlos_y",Run.globiter);
	    if(v == 3)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vrms_y",Run.globiter);
	    
	    fhandle=fopen(filename,"w");
	    
	    float header[6];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    header[3] = (float) Run.time;
	    header[4] = (float) lgTmin;
	    header[5] = (float) dellgT;
	    fwrite(header,sizeof(float),6,fhandle);
	  }
	
	  //slice_write(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,fhandle);
	  slice_write_rebin(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,rebin[d2],rebin[d3],fhandle);
	  
	  if(xz_rank == 0)
	    fclose(fhandle);
	}
      }
    }
      
    if(d1 == 2){    
      MPI_Reduce(los_sum_loc,los_sum,nout*nslvar*localsize,MPI_DOUBLE,MPI_SUM,iroot,
		 ZCOL_COMM);
    
      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+nslvar*localsize]; 
	if( rfac > 0.0 ){
	  los_sum[i+2*nslvar*localsize] /= rfac; 
	  los_sum[i+3*nslvar*localsize] /= rfac;
	}
      }

      for(i=0;i<nslvar*localsize;i++){
	rfac = los_sum[i+3*nslvar*localsize]-los_sum[i+2*nslvar*localsize]*los_sum[i+2*nslvar*localsize];
	los_sum[i+3*nslvar*localsize] = sqrt(max(rfac,0.0));
      }

      for(v=0;v<nout*nslvar*localsize;v++)
	io_buf[v] = (float) los_sum[v];

      if (zcol_rank == iroot){
	for (v=0;v<nout;v++){
	  if(xy_rank == 0) {
	    if(v == 0)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_fil_z",Run.globiter);
	    if(v == 1)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_dem_z",Run.globiter);
	    if(v == 2)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vlos_z",Run.globiter);
	    if(v == 3)
	      sprintf(filename,"%s%s.%06d",Run.path_2D,"corona_emission_adj_vrms_z",Run.globiter);
	
	    fhandle=fopen(filename,"w");

	    float header[6];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[d2]/rebin[d2];
	    header[2] = (float) Grid.gsize[d3]/rebin[d3];
	    header[3] = (float) Run.time;
	    header[4] = (float) lgTmin;
	    header[5] = (float) dellgT;
	    fwrite(header,sizeof(float),6,fhandle);
	  }
	
	  //slice_write(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,fhandle);
	  slice_write_rebin(Grid,0,&(io_buf[v*nslvar*localsize]),localsize,nslvar,d2,d3,rebin[d2],rebin[d3],fhandle);
	  
	  if(xy_rank == 0)
	    fclose(fhandle);
	}
      }
    }
   
    delete[] los_sum_loc;
    delete[] los_sum;
    delete[] io_buf;
  }

  delete[] tlev;

}

