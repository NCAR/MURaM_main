#include <fstream>
#include <mpi.h>
#include "analysis.H"
#include "precision.h"
#include "physics.H"
#include "grid.H"
#include "run.H"

using namespace std;

inline real max(real a, real b) { return a > b ? a : b; }
inline real min(real a, real b) { return a < b ? a : b; }
inline real sqr(real a) { return a*a; }

void AnalyzeSolution(const RunData& Run,
		     const GridData& Grid,
		     const PhysicsData& Physics) {
  register int i,j,k,d,node;

  static const int x=0,y=1,z=2;

  static const int ndim = Grid.NDIM;
  static const real gridvol = Grid.volume;
  static const real cellvol = Grid.cellvol;
  static const int* stride  = Grid.stride;
  static const real* dx     = Grid.dx;
  static const cState* U0   = Grid.U0;
  static const cState* U    = Grid.U;

  static const real ycarea  = dx[0]*dx[1]; // Cell area in y direction

  static const real kvisc   = 1./Physics.Re;  // Kinematic viscosity
  static const real mvisc   = 1./Physics.Rem; // Magnetic viscosity

  static real   rho,p,J2,VoB,divV;
  static Vector W,J,V,B;  // Vorticity, current, fields
  static Vector dV[3];    // Velocity derivative
  static Vector dB[3];    // B-field derivative
  static Vector aB;       // |B| estimate

  real KEG,KE = 0.; // Kinetic energy
  real MEG,ME = 0.; // Magnetic energy
  real IEG,IE = 0.; // Internal energy

  real KEFG,KEF = 0.; // Kinetic energy flux
  real MEFG,MEF = 0.; // Magnetic energy flux
  real IEFG,IEF = 0.; // Internal energy flux

  real PWG,PW = 0.; // Pressure work
  real LFG,LF = 0.; // Lorentz force work
  real VHG,VH = 0.; // Viscous heating
  real OHG,OH = 0.; // Ohmic heating

  real KHG,KH = 0.; // Cross helicity
  real VxBG,VxB = 0.; // |VxB|
  real JxBG,JxB = 0.; // |JxB|
  real Wmax,W2max = 0.; // Max vorticity
  real WxG,Wx = 0.; // x component of vorticity
  real WyG,Wy = 0.; // y component of vorticity
  real WzG,Wz = 0.; // z component of vorticity
  real Bmax,B2max = 0.; // Max magnetic field
  real Jmax,J2max = 0.; // Max current
  real JxG,Jx = 0.; // x component of vorticity
  real JyG,Jy = 0.; // y component of vorticity
  real JzG,Jz = 0.; // z component of vorticity
  real divBG,divB = 0.; // |Div B|
  real VsquareG, Vsquare = 0.; // |V|^2
  real V2l=0.;
  real Tmean=0; // horizontal mean of Temp
  real TflG,Tfl = 0.; // |T-Tmean|^2
  real pmaxG,pmax = 0.;
  real pminG,pmin=1.e10;
  real rhominG,rhomin = 1.e10;
  real V2maxG,V2max = 0.;

  const int ibeg = Grid.lbeg[0];
  const int iend = Grid.lend[0];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];
  const int kbeg = Grid.kbeg[2];
  const int kend = Grid.kend[2];
  const int bufsize = Grid.bufsize;

  if( Run.order==2 || Run.order==5 ) {

#pragma acc parallel loop collapse(3) \
 private(node, KEG, MEG, dV[:3], dB[:3], aB[:3], p, V, B, \
  divV, W, J, J2, VoB, d) \
 reduction(+:IE) reduction(+:KE) reduction(max:B2max) \
 reduction(+:PW) reduction(+:Wx) reduction(+:Wy) \
 reduction(+:Wz) reduction(max:W2max) reduction(+:Jx) \
 reduction(+:Jy) reduction(+:Jz) reduction(max:J2max) \
 reduction(+:OH) reduction(+:VxB) reduction(+:JxB) \
 reduction(+:LF) reduction(+:VH) reduction(+:divB) \
 reduction(+:KH) reduction(+:KEF) reduction(+:MEF) \
 reduction(+:IEF) reduction(+:ME) \
 present(Grid[:1], U[:bufsize], U0[:bufsize])
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++)
    for(i=ibeg; i<=iend; i++) {
      node = Grid.node(i,j,k);

      // Volumetric averages
      KEG = U0[node].KineticEnergy();
      MEG = U0[node].MagneticEnergy();
      IE += U0[node].e-KEG-MEG;
      KE += KEG; ME += MEG;

      // Compute derivatives
      for(d=0;d<ndim;d++) {
	dV[d] = (U[node+stride[d]].Velocity()-
		 U[node-stride[d]].Velocity())/(2.*dx[d]);
	dB[d] = (U[node+stride[d]].B-
		 U[node-stride[d]].B)/(2.*dx[d]);
	aB[d] = (fabs(U[node+stride[d]].B[d])+
		 fabs(U[node-stride[d]].B[d]))/(2.*dx[d]);
      }

      p = U[node].Pressure();
      V = U[node].Velocity();
      B = U[node].B;

      // Compute max B^2
      B2max = max(B2max,B.sqr());

      // Compute pressure work
      divV = dV[x].x+dV[y].y+dV[z].z;
      PW += p*divV;

      // Compute vorticity
      W = Vector(dV[y].z-dV[z].y,dV[z].x-dV[x].z,dV[x].y-dV[y].x);
      Wx += fabs(W.x);
      Wy += fabs(W.y);
      Wz += fabs(W.z);
      W2max = max(W2max,W.sqr()); // max vorticity^2

      // Compute current
      J = Vector(dB[y].z-dB[z].y,dB[z].x-dB[x].z,dB[x].y-dB[y].x);

      J2 = J.sqr();

      Jx += fabs(J.x);
      Jy += fabs(J.y);
      Jz += fabs(J.z);
      J2max = max(J2max,J2/aB.sqr()); // max current^2

      // Compute Ohmic heating
      OH += mvisc*J2;

      // Compute Lorentz force work
      VxB += (V^B).sqr();
      W    = (J^B);
      JxB += W.sqr();
      LF  += V*W;

      // Compute viscous heating
      VH += kvisc*(2.*(sqr(dV[x].x)+sqr(dV[y].y)+sqr(dV[z].z))+
		   sqr(dV[y].x+dV[x].y)+sqr(dV[z].x+dV[x].z)+
		   sqr(dV[z].y+dV[y].z)-(2./3.)*divV*divV);

      // Compute |divB|
      divB += fabs(dB[x].x+dB[y].y+dB[z].z)/(fabs(aB.x)+fabs(aB.y)+fabs(aB.z));

      // Compute cross helicity
      VoB = V*B;
      KH += VoB;

      // Compute boundary fluxes
      if( Grid.is_gbeg[1] && k == Grid.lbeg[1] ) {
	KEG  = 0.5*U[node].d*V.sqr();
	KEF += V.y*(KEG+p)-kvisc*(V.x*dV[x].y+V.z*dV[z].y-
				  (2./3.)*V.y*(dV[x].x+dV[z].z));

	MEG  = 0.5*B.sqr();
	MEF += 2.*V.y*MEG-B.y*VoB+mvisc*(B.x*dB[x].y+B.z*dB[z].y);

	IEF += V.y*(U[node].e-KEG-MEG);
      }

      if( Grid.is_gend[0] && k == Grid.lend[0] ) {
	KEG  = 0.5*U[node].d*V.sqr();
	KEF -= V.y*(KEG+p)-kvisc*(V.x*dV[x].y+V.z*dV[z].y-
				  (2./3.)*V.y*(dV[x].x+dV[z].z));

	MEG  = 0.5*B.sqr();
	MEF -= 2.*V.y*MEG-B.y*VoB+mvisc*(B.x*dB[x].y+B.z*dB[z].y);

	IEF -= V.y*(U[node].e-KEG-MEG);
      }
    }

    MPI_Reduce(&KE   ,&KEG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&ME   ,&MEG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&IE   ,&IEG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&KEF  ,&KEFG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&MEF  ,&MEFG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&IEF  ,&IEFG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&PW   ,&PWG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&LF   ,&LFG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&VH   ,&VHG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&OH   ,&OHG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&KH   ,&KHG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&VxB  ,&VxBG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&JxB  ,&JxBG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Wx   ,&WxG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Wy   ,&WyG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Wz   ,&WzG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&W2max,&Wmax,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&Jx   ,&JxG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Jy   ,&JyG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Jz   ,&JzG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&J2max,&Jmax,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&B2max,&Bmax,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&divB,&divBG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);

    if( Run.rank==0 ) {
      ofstream fptr(Run.anlfile,ios::out|ios::app);
      fptr.precision(10);
      fptr << Run.globiter << ' '
	   << Run.time << ' '
	   << Run.dt << ' '
	   << cellvol*KEG/gridvol << ' '
	   << cellvol*MEG/gridvol << ' '
	   << cellvol*IEG/gridvol << ' '
	   << ycarea*KEFG/gridvol << ' '
	   << ycarea*MEFG/gridvol << ' '
	   << ycarea*IEFG/gridvol << ' '
	   << cellvol*PWG/gridvol << ' '
	   << cellvol*LFG/gridvol << ' '
	   << cellvol*VHG/gridvol << ' '
	   << cellvol*OHG/gridvol << ' '
	   << cellvol*KHG/gridvol << ' '
	   << sqrt(cellvol*VxBG/gridvol) << ' '
	   << sqrt(cellvol*JxBG/gridvol) << ' '
	   << cellvol*WxG/gridvol << ' '
	   << cellvol*WyG/gridvol << ' '
	   << cellvol*WzG/gridvol << ' '
	   << sqrt(Wmax) << ' '
	   << cellvol*JxG/gridvol << ' '
	   << cellvol*JyG/gridvol << ' '
	   << cellvol*JzG/gridvol << ' '
	   << sqrt(Jmax) << ' '
	   << sqrt(Bmax) << ' '
	   << cellvol*divBG/gridvol << endl;
      fptr.close();
    }
  }

  /* This part needs a lot of work */
  else if( Run.order==4 ) {

#pragma acc parallel loop collapse(3) \
 private(node, KEG, MEG, dV[:3], dB[:3], aB[:3], p, V, B, \
  divV, W, J, J2, VoB, d, rho, V2l) \
 reduction(+:IE) reduction(+:KE) reduction(max:B2max) \
 reduction(+:PW) reduction(+:Wx) reduction(+:Wy) \
 reduction(+:Wz) reduction(max:W2max) reduction(+:Jx) \
 reduction(+:Jy) reduction(+:Jz) reduction(max:J2max) \
 reduction(+:OH) reduction(+:VxB) reduction(+:JxB) \
 reduction(+:LF) reduction(+:VH) reduction(+:divB) \
 reduction(+:KH) reduction(+:KEF) reduction(+:MEF) \
 reduction(+:IEF) reduction(+:ME) reduction(max:V2max) \
 reductino(+:Vsquare) reduction(max:pmax), reduction(min:rhomin) \
 reduction(min:pmin) \
 present(Grid[:1], U[:bufsize], U0[:bufsize])
    for(k=kbeg; k<=kend; k++)
    for(j=jbeg; j<=jend; j++)
    for(i=ibeg; i<=iend; i++) {
      node = Grid.node(i,j,k);

      // Volumetric averages
      KEG = U0[node].KineticEnergy();
      MEG = U0[node].MagneticEnergy();
      IE += U0[node].e-KEG-MEG;
      KE += KEG; ME += MEG;

      /*
      KEG = U0[node].KineticEnergy();
      MEG = U0[node].MagneticEnergy();
      IE += 0.75*(U0[node].e-KEG-MEG);
      KE += 0.75*KEG;
      ME += 0.75*MEG;
      */

      // Compute derivatives
      for(d=0;d<ndim;d++) {
	/*
	KEG = U0[node+stride[d]].KineticEnergy();
	MEG = U0[node+stride[d]].MagneticEnergy();
	IE += (U0[node+stride[d]].e-KEG-MEG)/24.;
	KE += KEG/24.;
	ME += MEG/24.;

	KEG = U0[node-stride[d]].KineticEnergy();
	MEG = U0[node-stride[d]].MagneticEnergy();
	IE += (U0[node-stride[d]].e-KEG-MEG)/24.;
	KE += KEG/24.;
	ME += MEG/24.;
	*/

	dV[d] = (-U[node+2*stride[d]].Velocity()
		 +8*U[node+stride[d]].Velocity()
		 -8*U[node-stride[d]].Velocity()
		 +U[node-2*stride[d]].Velocity())/(12*dx[d]);
	dB[d] = (-U[node+2*stride[d]].B
		 +8*U[node+stride[d]].B
		 -8*U[node-stride[d]].B
		 +U[node-2*stride[d]].B)/(12*dx[d]);
	aB[d] = (fabs(U[node+stride[d]].B[d])+
		 fabs(U[node-stride[d]].B[d]))/(2.*dx[d]);
      }

      rho = U[node].Density();
      p = U[node].Pressure();
      V = U[node].Velocity();
      B = U[node].B;

      // Compute max B^2
      B2max = max(B2max,B.sqr());

      // Compute pressure work
      divV = dV[x].x+dV[y].y+dV[z].z;
      PW += p*divV;

      // Compute vorticity
      W = Vector(dV[y].z-dV[z].y,dV[z].x-dV[x].z,dV[x].y-dV[y].x);
      Wx += fabs(W.x);
      Wy += fabs(W.y);
      Wz += fabs(W.z);
      W2max = max(W2max,W.sqr()); // max vorticity^2

      // Compute current
      J = Vector(dB[y].z-dB[z].y,dB[z].x-dB[x].z,dB[x].y-dB[y].x);
      J2 = J.sqr();

      Jx += fabs(J.x);
      Jy += fabs(J.y);
      Jz += fabs(J.z);
      // J2max = max(J2max,J2/aB.sqr()); // max current^2

      // Compute Ohmic heating
      OH += mvisc*J2;

      // Compute Lorentz force work
      VxB += (V^B).sqr();
      W    = (J^B);
      JxB += W.sqr();
      LF  += V*W;

      // Compute viscous heating
      VH += kvisc*(2.*(sqr(dV[x].x)+sqr(dV[y].y)+sqr(dV[z].z))+
		   sqr(dV[y].x+dV[x].y)+sqr(dV[z].x+dV[x].z)+
		   sqr(dV[z].y+dV[y].z)-(2./3.)*divV*divV);

      // Compute |divB|
      //  divB += fabs(dB[x].x+dB[y].y+dB[z].z)/(fabs(aB.x)+fabs(aB.y)+fabs(aB.z));

      // Compute cross helicity
      VoB = V*B;
      KH += VoB;

      // Compute Vsquare
      V2l = V.sqr();
      V2max = max(V2max,V2l);      
      Vsquare +=V2l;
      pmax = max(pmax,p);
      rhomin = min(rhomin,rho);
      pmin = min(pmin,p);    
  


      
// Compute boundary fluxes
      if( Grid.is_gbeg[1] && k == Grid.lbeg[1] ) {
	KEG  = 0.5*U[node].d*V.sqr();
	KEF += V.y*(KEG+p)-kvisc*(V.x*dV[x].y+V.z*dV[z].y-
				  (2./3.)*V.y*(dV[x].x+dV[z].z));

	MEG  = 0.5*B.sqr();
	MEF += 2.*V.y*MEG-B.y*VoB+mvisc*(B.x*dB[x].y+B.z*dB[z].y);

	IEF += V.y*(U[node].e-KEG-MEG);
      }

      if( Grid.is_gend[0] && k == Grid.lend[0] ) {
	KEG  = 0.5*U[node].d*V.sqr();
	KEF -= V.y*(KEG+p)-kvisc*(V.x*dV[x].y+V.z*dV[z].y-
				  (2./3.)*V.y*(dV[x].x+dV[z].z));

	MEG  = 0.5*B.sqr();
	MEF -= 2.*V.y*MEG-B.y*VoB+mvisc*(B.x*dB[x].y+B.z*dB[z].y);

	IEF -= V.y*(U[node].e-KEG-MEG);
      }
    }

    MPI_Reduce(&KE   ,&KEG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&ME   ,&MEG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&IE   ,&IEG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&KEF  ,&KEFG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&MEF  ,&MEFG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&IEF  ,&IEFG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&PW   ,&PWG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&LF   ,&LFG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&VH   ,&VHG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&OH   ,&OHG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&KH   ,&KHG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&VxB  ,&VxBG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&JxB  ,&JxBG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Wx   ,&WxG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Wy   ,&WyG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Wz   ,&WzG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&W2max,&Wmax,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&Jx   ,&JxG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Jy   ,&JyG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Jz   ,&JzG ,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&J2max,&Jmax,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&B2max,&Bmax,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&divB,&divBG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Vsquare,&VsquareG,1,REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&V2max  ,&V2maxG  ,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&pmax   ,&pmaxG   ,1,REALTYPE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&rhomin ,&rhominG ,1,REALTYPE,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(&pmin   ,&pminG   ,1,REALTYPE,MPI_MIN,0,MPI_COMM_WORLD);


    if( Run.rank==0 ) {
      ofstream fptr(Run.anlfile,ios::out|ios::app);
      fptr.precision(10);
      fptr << Run.globiter << ' '
	   << Run.time << ' '
	   << Run.dt << ' '
	   << cellvol*KEG/gridvol << ' '
	   << cellvol*MEG/gridvol << ' '
	   << cellvol*IEG/gridvol << ' '
	   << ycarea*KEFG/gridvol << ' '
	   << ycarea*MEFG/gridvol << ' '
	   << ycarea*IEFG/gridvol << ' '
	   << cellvol*PWG/gridvol << ' '
	   << cellvol*LFG/gridvol << ' '
	   << cellvol*VHG/gridvol << ' '
	   << cellvol*OHG/gridvol << ' '
	   << cellvol*KHG/gridvol << ' '
	   << sqrt(cellvol*VxBG/gridvol) << ' '
	   << sqrt(cellvol*JxBG/gridvol) << ' '
	   << cellvol*WxG/gridvol << ' '
	   << cellvol*WyG/gridvol << ' '
	   << cellvol*WzG/gridvol << ' '
	   << sqrt(Wmax) << ' '
	   << cellvol*JxG/gridvol << ' '
	   << cellvol*JyG/gridvol << ' '
	   << cellvol*JzG/gridvol << ' '
	   << sqrt(Jmax) << ' '
	   << sqrt(Bmax) << ' '
	   << cellvol*divBG/gridvol << ' '
           << cellvol*VsquareG/gridvol << ' '
	   << V2maxG  <<' '
	   << pmaxG   <<' '
	   << rhominG <<' '
	   << pminG   << endl;


     
      fptr.close();
    }
  }
}
