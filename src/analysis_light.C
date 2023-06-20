#include <mpi.h>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "analysis.H"
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "rt/rt.h"

using std::cout;
using std::endl;
using std::min;
using std::max;

void AnalyzeSolution(const RunData& Run,const GridData& Grid,const PhysicsData& Physics,RTS* rts) {

  register int i,j,k,node,node2;

  static const cState* U = Grid.U;

  double gsize = double(Grid.gsize[0])*double(Grid.gsize[1])*double(Grid.gsize[2]);
  double hsize = double(Grid.gsize[1])*double(Grid.gsize[2]);

  double vmaxl=0.0;
  double bmaxl=0.0;
  double rminl=1.e99;
  double rmaxl=0.0;
  double vrmsl=0.0;
  double brmsl=0.0;
  double abyml=0.0;
  
  double rho_min,rho_max,vmax,bmax,vrms,brms,abym,q1,q2;
 
  LOCAL_LOOP(Grid,i,j,k) {
    node  = Grid.node(i,j,k);
    node2 = Grid.node(i+1,j,k);
    
    rminl = min(rminl,U[node].d);
    rmaxl = max(rmaxl,U[node].d);
    vmaxl = max(vmaxl,sqrt(U[node].Velocity().sqr()));
    bmaxl = max(bmaxl,sqrt(U[node].B.sqr()));

    vrmsl+= U[node].Velocity().sqr();
    brmsl+= U[node].B.sqr();

    // interpolate Bx on tau=1 surface
    if( (Grid.Tau[node] >= 1.0) && (Grid.Tau[node2] < 1.0) ){
      q1 = (1.0-Grid.Tau[node2])/(Grid.Tau[node]-Grid.Tau[node2]);
      q2 = 1.0-q1;

      abyml += fabs(q1*U[node].B.x+q2*U[node2].B.x);
    }
  }

  MPI_Reduce(&rminl,&rho_min ,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&rmaxl,&rho_max ,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&vmaxl,&vmax    ,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&bmaxl,&bmax    ,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&vrmsl,&vrms    ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&brmsl,&brms    ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&abyml,&abym    ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if(Run.rank==0){ // MPI_Reduce results are only meaningful on rank 0!
    vrms  = sqrt(vrms/gsize);
    brms  = sqrt(brms/gsize);
    abym  = abym/hsize;

    // v in km/s
    vrms *= 1.e-5;
    vmax *= 1.e-5;

    // B in Gauss
    bmax *= 3.5449;
    brms *= 3.5449;
    abym *= 3.5449;

    std::ofstream fptr(Run.anlfile,std::ios::out|std::ios::app);
    fptr.precision(10);
    fptr << Run.globiter << ' '
	 << Run.time << ' '
	 << Run.dt << " | "
	 << vmax << ' '
	 << vrms << " | "
	 << bmax << ' '
	 << brms << ' '
         << abym << " | "
	 << rho_min << ' '
	 << rho_max << endl;
    
    fptr.close();
  }
  
}
