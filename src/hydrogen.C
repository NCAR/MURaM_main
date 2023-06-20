#include <mpi.h>
#include <iostream>
#include "grid.H"
#include "hydrogen_eos.H"

/* Hydrogen ionisation routine, following the prescription of Leenaarts 2011
Using the updated H2 rates from Gudiksen 
The Gist of this routine is:
At the beginning of an interation, while calling const_to_prim we 
additionally call hydrogen_eos::solve.
This will:
1) set up the rate equations.
2) Fill the jacobian.
3) Invert the jacobian.
4) Perform a newton-raphson update.
5) Iterate 1-4 until convergence criterion or max iterations are reached.
*/

const double NH = 4.407e23;
const double cs = 2.99792458e10;
const double pics = 8.0*pi/cs/cs;

hydrogen::hydrogen() {
  // First we need the 6 hydrogen levels, nH2, ne and T, Maybe this will go in a structure.
  n1 = NULL;
  n2 = NULL;
  n3 = NULL;
  n4 = NULL;
  n5 = NULL;
  n6 = NULL;
  nH = NULL;
  ne = NULL;
  T = NULL;
}
hydrogen::~hydrogen(){
  // Delete
  delete [] n1;
  delete [] n2;
  delete [] n3;
  delete [] n4;
  delete [] n5;
  delete [] n6;
  delete [] nH2;
  delete [] ne;
  delete [] T;
}
hydrogen::hydrogen_init() {
  // Allocate`
  n1 = new double [buffsize]();
  n2 = new double [buffsize]();
  n3 = new double [buffsize]();
  n4 = new double [buffsize]();
  n5 = new double [buffsize]();
  n6 = new double [buffsize]();
  nH2 = new double [buffsize]();
  ne = new double [buffsize]();
  T = new double [buffsize]();

  // Calculate frequently used constants

}

hydrogen::update(max_iter,hyd_tol){

//define stuff I need and read from memory

int hiter=0;
double error = 1.0e10
while((h_iter<=max_iter)||(h_err>=h_tol)){
// Calculate Rates as in Sollum
calc_rates();

// Fill jacobian with rate equations, chemical equilibrium conservation of hydrogen nuclei etc
fill_jacobian();

// Calculate Matrix Inversion
invert_jacobian();

//H = J*H_old

h_iter = h_iter+1;
// Calculate error
h_err = error();
}
}

hydrogen::calculate_rates(){

  double var[9][v_length], res[9][v_length];

  // First Loop, calculate Rate equations
  OUTER_LOOP(bounds,j,k,d2,d3){
    offset = j*stride[d2]+k*stride[d3];
    for(i=0;i<sz;i++){
      node = offset+i*str;

      // We need density for ntot and Temperature for the rates

      rho[i] = Grid.U[node].d;
      tem[i] = H.T[node];
     
      ne[i] = H.ne[node];
      ntot[i] =  rho[node]*NH;
      
      double sqrttem = power(tem[i],0.5);
      double invtem = 1.0/tem[i];
      double invne = 1.0/ne[i];
     
      // Calculate radiative and collisional rates as in Leenarts et al 2011

      for (int u=0;u<5);
        // bound bound rates
        for (int l=0;l<5);

          // First Calculate collisional rates
          // In the future, precompute glgu
          double glgu = g[l]/g[u];
          Cul[l][u][i] = glgu*ne[i]*Cexc[i]*sqrttem;
          double edif = -(xi[u]-xi[l])*intem/k;
          double expedif = exp(edif);
          Clu[l][u][i] = ne[i]*Cexc[i]*sqrttem*expedif;

          dCuldne[l][u][i] = Cul[l][u][i]*invne;
          dCludne[l][u][i] = Clu[l][u][i]*invne;
 
          dCuldT[l][u][i] = glgu*ne[i]*sqrttem[i]*expedif*(dCexcdT[i] + Cexc[i]*0.5*invtem);
          dCludT[l][u][i] = ne[i]*sqrttem[i]*expedif*(dCexcdT[i] + Cexc[i]*(edif*invtem+0.5*invtem));

          // Next radiative rates
       
          Rul[l][u][i] = glgu*rconst[l][u]*Rul_expTrad[i];
          Rlu[l][u][i] = rconst[l][u]*Rlu_expTrad[i];

          dRuldt[l][u][i] = glgu*drdtconst[l][u]*dRludT_expTrad[i]
          dRludt[l][u][i] = drdtconst[l][u]*dRuldT_expTrad[i];

       }
       // Bound free rates
       //
       // collisional

       double nin6 = ne[i]*0.5*glgu[l][u]*pow(2*pi*m_e*k*T[i]/h/h,1.5);
       double edif = -(xi[6]-xi[l])*invtem/k;
       double expedif = exp(edif);
       double glgu = g[l]/g[5];
       Clu[l][5][i] = ne[i]*Cion[i]*sqrttem*expedif; 
       Cul[l][5][i] = nin6*Clu[l][5][i];

       dCuldne[l][5][i] = 2.0*Cul[l][5][i]*invne;
       dCludne[l][5][i] = Clu[l][5][i]*invne;
       
       dCludT[l][5][i] = ne[i]*sqrttem[i]*expedif*(dCexcdT[i] + Cexc[i]*(edif*invtem+0.5*invtem));
       dCuldT[l][5][i] = (dCludT[l][5][i]- Cul[l][5][i]*invtem*(1.5+edif))*nin6;

       //Radiative
       Rul[l][5][i]=pics*alpha_0*pow(freq0,3)*infsum_A51(nh*freq0/k/Trad)
       Rlu[l][5][i]=pics*alpha_0*pow(freq0,3)*nin6*infsum_A54(nh*freq0/k/Trad)

       dRuldne[l][5][i]=Rul[l][5][i]/ne[i];
       dRludne[l][5][i]=0.0;
       
       dRuldt[l][5][i]=pics*alpha_0*pow(freq0,3)*infsum_A57(nh*freq0/k/Trad)
       dRludt[l][5][i]=pics*alpha_0*pow(freq0,3)*invtem;
     }
}

hydrogen::fill_jacobian(){

}

hydrogen::invert_jacobian(){

}

hyrogen::inf_sum_A51(){

}

hyrogen::inf_sum_A54(){

}

hyrogen::inf_sum_A57(){

}

hydrogen::error(){

}

