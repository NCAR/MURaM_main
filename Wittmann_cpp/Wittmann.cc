#include<stdio.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include <stdlib.h>
#include <pwd.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <ctime>
#include <cstdlib>
#include <cfloat>

#include "Wittmann.h"
#include "atom.h"

const double eos_gamma = 1.65;
const double eos_mu = 0.62;
const double eos_mu_n = 1.3;
const double c_temp= (eos_gamma-1.)*eos_mu/Rgas;
const double c_pres= (eos_gamma-1.);

float** array_2d_contiguous(int nx,int ny){
  float **array = new float*[nx];

  array[0] = new float[nx*ny];

  for(int i=0;i<nx;i++)
    array[i] = array[0]+i*ny;

  return array;
}

using namespace std;
  int main(){
    /* Calls Wittman routines, provided from Borrero, then down the chain:
    Robert Cameron (cameron@mps.mpg.de)
    L. S. Anusha (bhasari@mps.mpg.de)
    D.P. Version to make 3D EOS for Hion MURaM (przybylski@mps.mpg.de) */
        
    /* Minimum and Maximum T */
    double T_min, T_max;
    double pg,T,E, pe;
    double E_min, E_max,pe_a,pe_b;
    double E_old, E_a, E_b;
    double T_old, T_a, T_b;

    /* density grid */
    Nr = 3000;
    r_grid = new double [Nr];
    double rmin=log(1.0e-20);
    double rmax=log(1.0e-4);
    
    /* energy grid */
    Neps = 3000;
    eps_grid = new double [Neps];
    double emin=log(1.0e12);
    double emax=log(2.0e15);
    
    /* inverse grids */
    Np = 400;
    p_grid = new double [Np];
    /* Inverse grids only used for bottom boundary, so low pressures not needed */
    double pmin = 10.0;
    double pmax = 22.0;

    Ns = 300;
    s_grid = new double [Ns];
    double smin = -1.5e9;
    double smax = -5.0e8;

    double n_i[ncontr];
    
    ttbl = array_2d_contiguous(Nr,Neps);
    ptbl = array_2d_contiguous(Nr,Neps);
    netbl = array_2d_contiguous(Nr,Neps);
    rhontbl = array_2d_contiguous(Nr,Neps);
    rhoitbl = array_2d_contiguous(Nr,Neps);
    ambtbl = array_2d_contiguous(Nr,Neps);
    stbl = array_2d_contiguous(Nr,Neps);
    rhotbl = array_2d_contiguous(Ns,Np);
    epstbl = array_2d_contiguous(Ns,Np);

    /* create density grid */

    r_grid[0]=rmin; 
    delta_r=(rmax-rmin)/(Nr-1);
    inv_delta_r = 1.0/delta_r;

    for(int ir=1;ir<Nr-1;ir++)
      r_grid[ir]=r_grid[ir-1]+delta_r;
    
    r_grid[Nr-1]=rmax;
    
    /* create energy grid */

    eps_grid[0]=emin;
    delta_eps=(emax-emin)/(Neps-1);
    inv_delta_eps = 1.0/delta_eps;

    for(int ieps=1;ieps<Neps-1;ieps++)
      eps_grid[ieps]=eps_grid[ieps-1] + delta_eps;

    eps_grid[Neps-1]=emax;
    
    fprintf(stdout, "Read in abundances \n");
    double *eps0=read_abun(30);

    fprintf(stdout, "Read in ionisation energies \n");
    double ***param0=read_ie(30);
    
    /* create atom object */
    atom at(ncontr,nimax,eps0,param0);

    for (int aa=0;aa<30;aa++)
    {
      int ni=aa+1<10?aa+1:10;
      for(int ii=0;ii<ni;ii++)
        delete [] param0[aa][ii];
      delete [] param0[aa];
    }

    delete [] param0;
    delete [] eps0;

    int non_converge_count=0;
    int t_out_of_range_count=0;
    double T_search=0.0,E_search=0.0;

    fprintf(stdout, "Begin calculation of EOS tables\n");

    for(int ir=0;ir<Nr;ir++)
    {
      double rh=exp(r_grid[ir]);
      double xntot=0.0;

      /* computer Ntot for all included elements */
      for (int i=0;i<ncontr;i++)
      {
        n_i[i]=rh*at.perg[i];
        xntot+=n_i[i];
      }
      double xnhtot = n_i[0];
      double xnohtot = xntot-xnhtot;

      for(int ieps=0;ieps<Neps;ieps++)
      {
        double rhoi = 0.0,rhon=0.0;

        double ep = exp(eps_grid[ieps]);
      
        fprintf(stdout,"ir=%d rho=%e | ieps=%d eps=%e \n",ir,rh,ieps,ep);

        double factor = 0.2;
        int converged = 0;
        int minit=0;
        while ((!converged)&&(minit<1000))
        { 
          if (ieps > 1)
            T_min = exp(ttbl[ir][ieps-1]) - (1+factor)*fabs(exp(ttbl[ir][ieps-1])-exp(ttbl[ir][ieps-2]));
          else if (ir > 1)
            T_min = exp(ttbl[ir-1][ieps]) - (1+factor)*fabs(exp(ttbl[ir-1][ieps])-exp(ttbl[ir-2][ieps]));
          else if (ieps == 1)
            T_min=0.95*T_search;
          else if (ir == 1)
            T_min = exp(ttbl[ir-1][ieps])*0.95;
          else
            T_min = 1.0e3;

          T_min = max(1.0e3,T_min);

          T=T_min;
          T_a=T_min;
          if (debug_verbose > 0)
            fprintf(stdout,"------------  Tmin, Temp=%e -------------------------\n",T);
          


          /* initial guess for pe */
          if (ieps > 1)
            pe = exp(netbl[ir][ieps-1]) - (1+factor)*fabs(exp(netbl[ir][ieps-1])-exp(netbl[ir][ieps-2]));
          else if (ir > 1)
            pe = exp(netbl[ir-1][ieps]) - (1+factor)*fabs(exp(netbl[ir-1][ieps])-exp(netbl[ir-2][ieps]));
          else 
            pe = xntot;
          
          if (pe!=pe)
            pe = xntot;

          double NN = xntot + pe;

          pe*=k*T;

          double err_pe = 1.0e8;
          int pe_iter=0;
          int max_pe_iter = 1e6;

          double *pres = new double [7];

          while ((err_pe > pe_tol)&&(pe_iter<max_pe_iter))
          {
            for(int i=0;i<ncontr;i++)
              at.partf(i,T,NN);

            double * pres1=pe_pg10(at, T,pe,xnhtot*k*T,rh);

            err_pe = fabs(pe-pres1[0])/pe;
            pe = pres1[0];
            pe_iter+=1;

            NN = xnohtot+(pres1[1]+pres1[2]+pres1[3]+pres1[4]+pres1[5])*xnhtot+pe/(k*T);

            memcpy(pres,pres1,7*sizeof(double));
            delete [] pres1;
          }

          double xnhtot1=0.0;
          for (int i=1;i<4;i++)
            xnhtot1+=pres[i];
          for (int i=4;i<6;i++)
            xnhtot1+=pres[i]*2;
          xnhtot1*=xnhtot;

          pe_a = pe;

          pg = pres[6];

          if (debug_verbose >0)
            fprintf(stdout,"Hn_1=%e,Hp_1=%e,Hm_1=%e H2n_1=%e H2p_1=%e ne=%e| sumH/xnhtot=%lf \n",
                pres[1]*xnhtot,pres[2]*xnhtot,pres[3]*xnhtot,pres[4]*xnhtot,pres[5]*xnhtot,pres[0]/(k*T),xnhtot1/xnhtot);
            
          /* molecular H2, H2+ and H- energies */
          double * Hmol_eps = Heps(T);
          double sumi = Hmol_eps[0]*pres[1] + Hmol_eps[1]*pres[2] + Hmol_eps[2]*pres[3] + Hmol_eps[3]*pres[4] + Hmol_eps[4]*pres[5];
          delete [] Hmol_eps;

          /* Excitation energies */
          for(int i=1;i<ncontr;i++)
          {
            double chisum = 0.0;
            for (int j=1;j<at.nion[i];j++)
            {
              chisum += at.chi[i][j];
              sumi+=at.nlvl[i][j]*chisum*ev;
            }
          }

          sumi*=xnhtot;
          
          double ski = 1.5*pg;

          E_min=ski + sumi;
         
          /* Diffusion Approximation Radiation
           * if (pg > 1.0e3)
           * E_min+=8.0/15.0*PI*(PI*k*T)*pow(PI*k*T/(C*HP),3); */

          E_a=ep-E_min/rh;

          if (E_a < 0)
            factor*=2;
          else if ( E_a!=E_a)
            factor*=0.9;
          else
            converged=1;

          minit+=1;

          delete [] pres;
        }
        
        factor = 0.5;
        converged = 0;
        int maxit = 0;
        while ((!converged)&&(maxit<1000))
        {
          if (ieps > 1)
            T_max = exp(ttbl[ir][ieps-1]) + (1.0+factor)*(exp(ttbl[ir][ieps-1])-exp(ttbl[ir][ieps-2]));
          else if (ir > 1)
            T_max = exp(ttbl[ir-1][ieps]) + (1.0+factor)*(exp(ttbl[ir-1][ieps])-exp(ttbl[ir-2][ieps]));
          else if (ieps == 1)
            T_max=1.5*T_search;
          else if (ir == 1)
            T_max = exp(ttbl[ir-1][ieps])*1.25;
          else
            T_max = 1.0e7;

          T_max = min(1.0e7,T_max);

          T=T_max;
          T_b=T_max;
          if (debug_verbose > 0)
            fprintf(stdout,"------------  Tmax, Temp=%e -------------------------\n",T);


          /* initial guess for pe */
          if (ieps > 1)
            pe = exp(netbl[ir][ieps-1]) + (1+factor)*(exp(netbl[ir][ieps-1])-exp(netbl[ir][ieps-2]));
          else if (ir > 1)
            pe = exp(netbl[ir-1][ieps]) + (1+factor)*(exp(netbl[ir-1][ieps])-exp(netbl[ir-2][ieps]));
          else 
            pe = xntot;

          if (pe!=pe)
            pe = xntot;

          double NN = xntot + pe;

          pe*=k*T;


          double err_pe = 1.0e8;
          int pe_iter=0;
          int max_pe_iter = 1e6;
          
          double *pres = new double [7];

          while ((err_pe > pe_tol)&&(pe_iter<max_pe_iter))
          {
            for(int i=0;i<ncontr;i++)
              at.partf(i,T,NN);

            double * pres1=pe_pg10(at, T,pe,xnhtot*k*T,rh);

            err_pe = fabs(pe-pres1[0])/pe;
            pe = pres1[0];
            pe_iter+=1;

            NN = xnohtot+(pres1[1]+pres1[2]+pres1[3]+pres1[4]+pres1[5])*xnhtot+pe/(k*T);

            memcpy(pres,pres1,7*sizeof(double));
            delete [] pres1;
          }
            
          double xnhtot1=0.0;
          for (int i=1;i<4;i++)
            xnhtot1+=pres[i];
          for (int i=4;i<6;i++)
            xnhtot1+=pres[i]*2;
          xnhtot1*=xnhtot;

          pe_b = pe;
          
          pg = pres[6];

          if (debug_verbose > 0)
            fprintf(stdout,"Hn_1=%e,Hp_1=%e,Hm_1=%e H2n_1=%e H2p_1=%e ne=%e| sumH/xnhtot=%lf \n",
                pres[1]*xnhtot,pres[2]*xnhtot,pres[3]*xnhtot,pres[4]*xnhtot,pres[5]*xnhtot,pres[0]/(k*T),xnhtot1/xnhtot);
          
          /* molecular H2, H2+ and H- energies */
          double * Hmol_eps = Heps(T);
          double sumi = Hmol_eps[0]*pres[1] + Hmol_eps[1]*pres[2] + Hmol_eps[2]*pres[3] + Hmol_eps[3]*pres[4] + Hmol_eps[4]*pres[5];
          delete [] Hmol_eps;

          /* Excitation energies */
          for(int i=1;i<ncontr;i++)
          {
            double chisum = 0.0;
            for (int j=1;j<at.nion[i];j++)
            {
              chisum += at.chi[i][j];
              sumi+=at.nlvl[i][j]*chisum*ev;
            }
          }
          sumi*=xnhtot;
          
          double ski = 1.5*pg;
          
          E_max=ski + sumi;

          /* Diffusion Approximation Radiation
           * if (pg > 1.0e3)
           * E_max+= 8.0/15.0*PI*(PI*k*exp(T))*pow(PI*k*T/(C*HP),3); */

          E_b=ep-E_max/rh;

          if (E_b > 0)
            factor*=2;
          else if (E_b!=E_b)
            factor*=0.9;
          else 
            converged=1;

          maxit+=1;

          delete [] pres;
        }

        T = -E_a*(T_b-T_a)/(E_b-E_a)+T_a;

        fprintf(stdout,"Initial guess Temp=%e | ",T);
        fprintf(stdout,"eps_min=%e eps_max =%e | T_min=%e T_max=%e | pe_min=%e pe_max=%e\n",E_min/rh,E_max/rh,T_min,T_max,pe_a,pe_b);

        /* Initialize */
        E=E_min; 

        if (debug_verbose > 0)
          fprintf(stdout,"------------  Begin iterations -------------------------\n");
        if((E_a > 0. && E_b > 0.)||(E_a < 0. && E_b < 0.))
        {
          fprintf(stdout,"E_a and E_b out of bounds\n");
          fprintf(stdout,"Use ideal gas for T, pg and PV=NRT for electron number\n");
          t_out_of_range_count+=1;
          T = c_temp*ep;
          T_search = T;
          pg = c_pres*ep/rh;

          /* pV=nRT */
          pe = pg-rh/(eos_mu_n*mh)*T*k; 
        }
        else 
        {
          /* Temperature Search */

          double rel_err=1.e+19;
          int iterTn=0;
          int itermax=100000;
          double * pres = new double[7];
          double xnhtot1 = 0.0;

          while(rel_err > search_tol && iterTn < itermax)
          {
            iterTn=iterTn+1;


            /* initial guess for pe */
            pe = -E_a*(pe_b-pe_a)/(E_b-E_a)+pe_a;
          
            if (pe!=pe)
              pe = xntot;

        
            double NN = xntot + pe/(k*T);

            double err_pe = 1.0e8;
            int pe_iter = 0;
            int max_pe_iter = 1e3;
            
            while ((err_pe > pe_tol)&&(pe_iter<max_pe_iter))
            {
              for(int i=0;i<ncontr;i++)
                at.partf(i,T,NN);

              double *pres1=pe_pg10(at, T,pe,xnhtot*k*T,rh);

              err_pe = fabs(pe-pres1[0])/pe;
              pe = pres1[0];
              pe_iter+=1;

              NN = xnohtot+(pres1[1]+pres1[2]+pres1[3]+pres1[4]+pres1[5])*xnhtot+pe/(k*T);

              memcpy(pres,pres1,7*sizeof(double));
              delete [] pres1;
            }

            xnhtot1=0.0;
            for (int i=1;i<4;i++)
              xnhtot1+=pres[i];
            for (int i=4;i<6;i++)
              xnhtot1+=pres[i]*2;
            xnhtot1*=xnhtot;

            pg = pres[6];

            /* Ion density */
            rhoi=(pres[2] + pres[3] + 2.0*pres[5])*at.W[0];
            rhon=(pres[1]+2.0*pres[4])*at.W[0];

            /* molecular H2, H2+ and H- energies */
            double * Hmol_eps = Heps(T);
            double sumi = Hmol_eps[0]*pres[1] + Hmol_eps[1]*pres[2] + Hmol_eps[2]*pres[3] + Hmol_eps[3]*pres[4] + Hmol_eps[4]*pres[5];
            delete [] Hmol_eps;

            /* Excitation energies */
            for(int i=1;i<ncontr;i++)
            {
              double chisum = 0.0;
                for (int j=1;j<at.nion[i];j++)
                {
                  chisum += at.chi[i][j];
                  sumi+=at.nlvl[i][j]*chisum*ev;
                  rhoi+=at.nlvl[i][j]*at.W[i];
                }
              rhon+=(at.nlvl[i][0])*at.W[i];
            }
            sumi*=xnhtot;
            
            rhoi*=amu*xnhtot;
            rhon*=amu*xnhtot;

            double ski=1.5*pg;

            E_old=E; 
            E=ski+sumi;

            /* Add Diffusion approximation Radiation

            double E_rad = 8.0/15.0*PI*(PI*k*T)*pow(PI*k*T/(C*HP),3);
            double p_rad = E_rad/3.0;

            if (pg > 1.0e3)
            {
              pg+= p_rad;
              E+=E_rad;
            }
            */

            /* Difference in energy */
            E_search=ep-E/rh;          

            T_old = T;
            int sw = (int) (E_search*E_a < 0.0);
            T_b = sw*T + (1-sw)*T_b;
            E_b = sw*E_search + (1-sw)*E_b; 
            pe_b = sw*pe + (1-sw)*pe_b;

            int sw2 = (int) (E_search*E_b < 0.0);
            T_a = sw2*T + (1-sw2)*T_a;
            E_a = sw2*E_search + (1-sw2)*E_a;
            pe_a = sw2*pe + (1-sw2)*pe_a;

            //T_search = 0.5*(T_a+T_b);
            T_search = -E_a*(T_b-T_a)/(E_b-E_a)+T_a;
            rel_err=fabs((T_search-T_old)/T_old);
            T=T_search;

            if (debug_verbose>0) 
            {
              fprintf(stdout,"sumi=%e ski=%e pg=%e pe=%e\n",sumi,ski,pg,pe);
              fprintf(stdout,"T_a=%e T_search=%e T_b=%e | E_a=%e E_search=%e E_b=%e\n",T_a,T_search,T_b,E_a,E_search,E_b);
              fprintf(stdout,"T_old=%e T_new=%e | eps_old=%e eps_new=%e\n",T_old,T,E_old/rh,E/rh);
              fprintf(stdout,"iterTn=%d rel_err=%e\n",iterTn,rel_err);
              fprintf(stdout,"-------------------------------------------------\n");
            }
          }
          if(iterTn==itermax)
          {
            fprintf(stdout,"Bisection method did not find a root\n");
            non_converge_count+=1;
          }
          fprintf(stdout,"Converged to solution | ");
          fprintf(stdout,"Hn_1=%e,Hp_1=%e,Hm_1=%e H2n_1=%e H2p_1=%e ne=%e| sumH/xnhtot=%lf \n",
              pres[1]*xnhtot,pres[2]*xnhtot,pres[3]*xnhtot,pres[4]*xnhtot,pres[5]*xnhtot,pres[0]/(k*T),xnhtot1/xnhtot);
          fprintf(stdout,"T=%e pg=%e rhoi=%e rhon=%e ne=%e\n",T,pg,rhoi, rhon,pe/k/T);
          fprintf(stdout," \n");

          delete [] pres;
        }
        ttbl[ir][ieps]= (float) log(T);
        ptbl[ir][ieps]= (float) log(pg);
        rhoitbl[ir][ieps]=(float) log(rhoi);
        rhontbl[ir][ieps] =(float) log(rhon);
        ambtbl[ir][ieps]=(float) 0.0;
        netbl[ir][ieps]=(float) log(pe/k/T);
      }
    }
      
    for (int ir=0;ir<Nr;ir++)
    {
      double ss1=0.0;
      stbl[ir][0] = 0.0;

      double hon2 = delta_eps/2.0;
      double hon3 = delta_eps/3.0;

      stbl[ir][1] = hon2*(exp(eps_grid[0])/exp(ttbl[ir][0])+exp(eps_grid[1])/exp(ttbl[ir][1]));

      ss1=hon3*(exp(eps_grid[0])/exp(ttbl[ir][0])+4.0*exp(eps_grid[1])/exp(ttbl[ir][1]));

      for (int ieps=2;ieps<Neps;ieps+=2)
      {
        double ss2=ss1+hon3*exp(eps_grid[ieps])/exp(ttbl[ir][ieps]);
        stbl[ir][ieps]=ss2;
        ss2+=hon2*(exp(eps_grid[ieps])/exp(ttbl[ir][ieps])+
            exp(eps_grid[ieps+1])/exp(ttbl[ir][ieps+1]));
        stbl[ir][ieps+1]=ss2;
        ss1+=hon3*(2.0*exp(eps_grid[ieps])/exp(ttbl[ir][ieps])+4.0*exp(eps_grid[ieps+1])/exp(ttbl[ir][ieps+1]));
      }
    }

    double ss1=0.0;
    for (int ieps=0;ieps<Neps;ieps++)
    {
      stbl[0][0]+=ss1;
    }

    double hon2 = delta_r/2.0;
    double hon3 = delta_r/3.0;

    double ss2 = -hon2*(exp(ptbl[0][0])/(exp(ttbl[0][0])*exp(r_grid[0]))+exp(ptbl[1][0])/(exp(ttbl[1][0])*exp(r_grid[1])));

    for (int ieps=0;ieps<Neps;ieps++)
    {
      stbl[1][ieps]+=ss2;
    }

    ss1-=hon3*(exp(ptbl[0][0])/(exp(ttbl[0][0])*exp(r_grid[0]))+4.0*exp(ptbl[1][0])/(exp(ttbl[1][0])*exp(r_grid[1])));

    for (int ir=2;ir<Nr;ir+=2)
    {
      ss2=ss1-hon3*exp(ptbl[ir][0])/(exp(ttbl[ir][0])*exp(r_grid[ir]));
      double ss3=ss2-hon2*(exp(ptbl[ir][0])/(exp(ttbl[ir][0])*exp(r_grid[ir]))+exp(ptbl[ir+1][0])/(exp(ttbl[ir+1][0])*exp(r_grid[ir+1])));
      ss1-=hon3*(2.0*exp(ptbl[ir][0])/(exp(ttbl[ir][0])*exp(r_grid[ir]))+4.0*exp(ptbl[ir+1][0])/(exp(ttbl[ir+1][0])*exp(r_grid[ir+1])));
      for (int ieps=0;ieps<Neps;ieps++)
      {
        stbl[ir][ieps]+=ss2;
        stbl[ir+1][ieps]+=ss3;
      }
    }

    fprintf(stdout,"non converge count = %i | t out of range count %i \n",non_converge_count, t_out_of_range_count);
    fprintf(stdout, "setting up inverse tables, smin=%e, smax=%e, pmin=%e,pmax=%e",smin,smax,pmin,pmax);


    p_grid[0]=pmin;
    delta_p=(pmax-pmin)/(Np-1);
    inv_delta_p=1.0/delta_p;

    for(int ip=1;ip<Np-1;ip++)
      p_grid[ip]=p_grid[ip-1]+delta_p;

    p_grid[Np-1]=pmax;

    s_grid[0]=smin;
    delta_s=(smax-smin)/(Ns-1);
    inv_delta_s = 1.0/delta_s;

    for(int is=1;is<Ns-1;is++)
      s_grid[is]=s_grid[is-1]+delta_s;

    s_grid[Ns-1]=smax;
    
    /* Invert for rho and eps tables */
    /* NB Inverse tables are only currently used for the lower boundary conditions.
     * As such we set the min/max rho/eps based on what is expected between a range of depths from -1 Mm
     * the tables currently do not work to a high enough density to do the the 20 Mm deep boxes. Around
     * 3.0e-3 would be needed */
    

    double ee0 = 1.0e12;
    double ee1 = 1.0e14;
    double ss0 = smin;
    ss1 = smax;
    double rr0 = 1.0e-10;
    double rr1 = 1.0e-4;
    double pp0 = exp(pmin);
    double pp1 = exp(pmax);

    fprintf(stdout, "Starting Inversion \n ");
    fprintf(stdout, " -------------------------------------- \n ");
      

    double max_err_p = 0.0;
    double max_err_s = 0.0;
 
    /* Quick First iteration to get approximate solution */
    for(int is=0;is<Ns;is++)
    {
      int ip=0;

      double p0 = exp(p_grid[ip]);
      double s0 = s_grid[is];

      double rh = log(1e-8);
      double ep = log(1.0e13);

      double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1,1.01,1000,0.001);

      rhotbl[is][ip] = output[0];
      epstbl[is][ip] = output[1];
      max_err_p = max(max_err_p,output[2]);
      max_err_s = max(max_err_s,output[3]);
      delete [] output;
    }
    fprintf(stdout, "First iteration, ip= 0, max_err_p %e, max_err_s %e \n",max_err_p,max_err_s);
    for(int ip=1;ip<Np;ip++)
    {
      max_err_p = 0.0;
      max_err_s = 0.0;
      double p0 = exp(p_grid[ip]);
      for(int is=0;is<Ns;is++)
      {
        double s0 = s_grid[is];
        
        double rh = rhotbl[is][ip-1];
        double ep = epstbl[is][ip-1];
      
        double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1,1.01,1000,0.001);

        rhotbl[is][ip] = output[0];
        epstbl[is][ip] = output[1];
        max_err_p = max(max_err_p,output[2]);
        max_err_s = max(max_err_s,output[3]);
        delete [] output;
      }
    fprintf(stdout, "First iteration, ip= %d, max_err_p %e, max_err_s %e \n",ip,max_err_p,max_err_s);
    }
    fprintf(stdout, "First Iteration Complete Starting Inversion Iteration 2\n ");
    fprintf(stdout, " -------------------------------------- \n ");
    /* Second iteration - Average and smooth in p to remove stripes:
     * Should probably iterate edges as well, but it is the edge, so likely doesn't matter */
    for(int ip=1;ip<Np-1;ip++)
    {
      max_err_p = 0.0;
      max_err_s = 0.0;
      double p0 = exp(p_grid[ip]);
      for(int is=0;is<Ns;is++)
      {
        double s0 = s_grid[is];

        double rh = (rhotbl[is][ip+1]+rhotbl[is][ip-1])/2.0;
        double ep = (epstbl[is][ip+1]+epstbl[is][ip-1])/2.0;
        
        double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1,1.01,10000,0.0001);

        rhotbl[is][ip] = output[0];
        epstbl[is][ip] = output[1];
        max_err_p = max(max_err_p,output[2]);
        max_err_s = max(max_err_s,output[3]);
        delete [] output;
      }
      fprintf(stdout, "Second iteration, ip= %d, max_err_p %e, max_err_s %e \n",ip,max_err_p,max_err_s);
    }


    for(int is=0;is<Ns;is++)
    {
      int ip = 0;
      double s0 = s_grid[is];
      double p0 = exp(p_grid[ip]);

      double rh = rhotbl[is][ip+1];
      double ep = epstbl[is][ip+1];
      double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1,1.01,10000,0.0001);

      rhotbl[is][ip] = output[0];
      epstbl[is][ip] = output[1];
      max_err_p = max(max_err_p,output[2]);
      max_err_s = max(max_err_s,output[3]);
      delete [] output;
    }

    for(int is=0;is<Ns;is++)
    {
      int ip = Np-1;
      double s0 = s_grid[is];
      double p0 = exp(p_grid[ip]);

      double rh = rhotbl[is][ip-1];
      double ep = epstbl[is][ip-1];

      double * output = invert_eos_newton(p0,s0,rh,ep,pp0,pp1,ss0,ss1,ee0,ee1,rr0,rr1,1.01,10000,0.0001);

      rhotbl[is][ip] = output[0];
      epstbl[is][ip] = output[1];
      max_err_p = max(max_err_p,output[2]);
      max_err_s = max(max_err_s,output[3]);
      delete [] output;
    }

    FILE *outfp;
    float header[14];

    header[0] = (float) Neps;
    header[1] = (float) Nr;
    header[2] = (float) Np;
    header[3] = (float) Ns;
    header[4] = (float) eps_grid[0];
    header[5] = (float) eps_grid[Neps-1];
    header[6] = (float) r_grid[0];
    header[7] = (float) r_grid[Nr-1];
    header[8] = (float) p_grid[0];
    header[9] = (float) p_grid[Np-1];
    header[10] = (float) s_grid[0];
    header[11] = (float) s_grid[Ns-1];
    header[12] = (float) eps_off;
    header[13] = (float) ss_off;

    outfp = fopen("./eostest.dat","w");
    fwrite(header,sizeof(float),14,outfp);
    fwrite(&ptbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&ttbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&stbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&netbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&rhoitbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&ambtbl[0][0],sizeof(float),Neps*Nr,outfp);
    fwrite(&epstbl[0][0],sizeof(float),Np*Ns,outfp);
    fwrite(&rhotbl[0][0],sizeof(float),Np*Ns,outfp);
    fclose(outfp);

    fprintf(stdout,"Deallocating tables\n");


    delete [] ptbl[0];
    delete [] ttbl[0];
    delete [] stbl[0];
    delete [] netbl[0];
    delete [] rhoitbl[0];
    delete [] rhontbl[0];
    delete [] ambtbl[0];

    delete [] ptbl;
    delete [] ttbl;
    delete [] stbl;
    delete [] netbl;
    delete [] rhoitbl;
    delete [] rhontbl;
    delete [] ambtbl;

    delete [] epstbl[0];
    delete [] rhotbl[0];

    delete [] epstbl;
    delete [] rhotbl;
    fprintf(stdout,"Deallocation complete\n");

  }

/* PE_PG10 (included in EQUISUBMU) evaluates the electonic pressure from the
 * the gaseous pressure and an estimate of the electronic pressure
 * calculates the electronic pressure from the pg and from an estimate of the pe1 */

double* pe_pg10(const atom &at, double t,double pe,double phtot, double rr)
{
  double g1,g2,g3,g4,g5;
  double f1,f2,f3,f4,f5,f6;

  if(t<500.){
    fprintf(stdout,"pe_pg10: temperature < 500 K; temperature = 500 K\n");
    t=500.;
  }
  double theta=5040./t;

  molecb mol(nmol,theta);
  double * cmol=new double [nmol]();
  double * dcmol=new double [nmol]();

  for (int i=0;i<nmol;i++){
    cmol[i]=mol.Y[i+1];
    dcmol[i]=mol.dY[i+1];
  }

  /* to set min vlaue for pe */
  if((pe<=1.0d-42)||(pe!=pe))
  {
    pe=1.0d-42;
    g4=0.;
    g5=0.;
  } 
  else
  {
    g4 = pe*pow(10.,cmol[0]);
    g5 = pe*pow(10.,cmol[1]);
  }

  double a,b,c,d,e;

  g1=0.;

  /* for all other elements except hydrogen, assume LTE and call saha */
  /* Determine how much of each is ionised */
  double epssum = 0.0;
  for(int i=1;i<ncontr;i++)
  {
    a=saha(theta,at.chi[i][1],at.uu[i][0],at.uu[i][1],pe);
    b=saha(theta,at.chi[i][2],at.uu[i][1],at.uu[i][2],pe);
    c=at.abu[i]/(1.+a*(1.+b));
    epssum+=at.abu[i];

    at.nlvl[i][0] = c;
    at.nlvl[i][1] = c*a;
    at.nlvl[i][2] = c*a*b;

    g1+=at.nlvl[i][1]+2.0*at.nlvl[i][2];
  }

  /* p(h+)/p(h) */
  g2=saha(theta,at.chi[0][1],at.uu[0][0],at.uu[0][1],pe);

  /* p(h)/p(h-) */
  g3 = saha(theta,inz_Hm1,1.0,at.uu[0][0],pe);
  
  if (g3 > 0.0)
    g3=1.0/g3;
  else
    g3 = 0.0;

  if (g2 < 1.0e5)
  {
    /* step 1: LTE f1 and f2:later replace by NE */
    a=1.+g2+g3;
    b=2.*(1.+g2/g5*g4);
    c=g5;
    d=g2-g3;
    e=g2/g5*g4;

    double c1=c*b*b+a*d*b-e*a*a;
    double c2=2.0*a*e-d*b+a*b*g1;
    double c3=-(e+b*g1);

    f1=0.5*c2/c1;
    double sgn=c1>0.0?1.00:-1.00;
    f1=-f1+sgn*sqrtl(f1*f1-c3/c1);
    f2=g2*f1;

    g2 = f2/f1;
    e=g2/g5*g4;

    f3=g3*f1;
    f5=(1.0-a*f1)/b;
    f4=e*f5;
    f6=f2-f3+f4+g1;
  } 
  else
  {
    /* If almost totally ionised then ignore H2 */
    f1 = 1.0/(1.0+g2+g3);
    f2 = g2/(1.0+g2+g3);
    f3 = g3/(1.0+g2+g3);
    f4 = 0.0;
    f5 = 0.0;
    f6 = g1 + f2-f3;
  }

  double * pres = new double [7];

  /* Return pe, ph/ph', php/ph', phm/ph', ph2/ph',ph2p/ph', pg/ph' */
  pres[0] = phtot*f6;
  pres[1] = f1;
  pres[2] = f2;
  pres[3] = f3;
  pres[4] = f5;
  pres[5] = f4;
  pres[6] = phtot*(epssum+f1+f2+f3+f4+f5+f6);

  delete [] cmol;
  delete [] dcmol;

  return pres;
}

double acota(double x,double x0,double x1){
  double xx=x;
  if(x<x0)xx=x0;
  if(x>x1)xx=x1;
  return xx;
}

double acotasig(double x,double x0,double x1){
  double xx=x;
  if(x<0.){
    x=-x;
    xx=acota(x,x0,x1);
    xx=-xx;
  }
  else{
    xx=acota(x,x0,x1);
  }
  return xx;
}

double * Heps(double tt)
{
  
  double* E = new double [5];

  double theta = 5040.39/tt;

  /* Vardya Polynomial fits to molecular energies */
  double dE_H2 = 0.0; /*2.6757-1.4772*theta+0.60602*theta*theta-0.12427*theta*theta*theta
    +0.0097503*theta*theta*theta*theta;
  dE_H2*=k*tt;*/

  double dE_H2p = 0.0; /* 2.9216-2.0036*theta+1.7231*theta*theta-0.82685*theta*theta*theta
    +0.15253*theta*theta*theta*theta;
  dE_H2p*=k*tt; */

  /* Hn */
  E[0] = (0.5*D0_h2)*ev;
  /* Hp */
  E[1] = (0.5*D0_h2+inz_H)*ev;
  /* Hm */
  E[2] = (0.5*D0_h2-inz_Hm1)*ev;
  /* H2 */
  E[3] = dE_H2;
  /* H2p */
  E[4] = (dE_H2p+D0_h2-D0_h2p+inz_H)*ev;

  return E;
}

double * invert_eos_newton(double p0,double s0,double rho, double eps,double pmin,double pmax,double smin,double smax,double emin,double emax,double rmin,double rmax,double fct, int maxiter, double tol)
{

  double err_p = 1.0e10;
  double err_s = 1.0e10;
  int it=0;
   
  double rho1 = exp(rho);
  double eps1 = exp(eps);
    
  double p1 = p_interp(eps1,rho1);
  double s1 = s_interp(eps1,rho1);

  while(((err_p > tol)&&(err_s > tol))&&(it<maxiter))
  {
    double p_r = (p_interp(eps1,fct*rho1)-p1)/((fct-1)* rho1);
    double p_e = (p_interp(fct*eps1,rho1)-p1)/((fct-1)* eps1);
    double s_r = (s_interp(eps1,fct*rho1)-s1)/((fct-1)* rho1);
    double s_e = (s_interp(fct*eps1,rho1)-s1)/((fct-1)* eps1);

    double det = p_r*s_e-p_e*s_r;

    double r1 = ( s_e*(p0-p1)-p_e*(s0-s1))/det;
    double e1 = (-s_r*(p0-p1)+p_r*(s0-s1))/det;

    rho1 = min(max(rmin,rho1+r1),rmax);
    eps1 = min(max(emin,eps1+e1),emax);

    p1 = p_interp(eps1,rho1);
    s1 = s_interp(eps1,rho1);

    err_p = fabs(p1-p0)/p0;
    err_s = fabs(s1-s0)/s0;

    it+=1;
  }

  double * output = new double [4];
  output[0] = log(rho1);
  output[1] = log(eps1);
  output[2] = err_p;
  output[3] = err_s;

  return output;
}

double bilinear(int n1, int n2,double *x,double *y,float **fxy,double xx,double yy)
{

  double c00,c01,c10,c11;
  double ff=0.;

  for (int i=0; i<n1-1; i++)
  {
    for (int j=0; j<n2-1; j++)
    {

      if(((xx>=x[i] && xx<x[i+1])||(xx>x[i] && xx<=x[i+1]))&&
          ((yy>=y[j] && yy<y[j+1])||(yy>y[j] && yy<=y[j+1])))
      {
        c00=(x[i+1]-xx)*(y[j+1]-yy)/((x[i+1]-x[i])*(y[j+1]-y[j]));
        c01=(xx-x[i])*(y[j+1]-yy)/((x[i+1]-x[i])*(y[j+1]-y[j]));
        c10=(x[i+1]-xx)*(yy-y[j])/((x[i+1]-x[i])*(y[j+1]-y[j]));
        c11=(xx-x[i])*(yy-y[j])/((x[i+1]-x[i])*(y[j+1]-y[j]));
        ff=c00*fxy[i][j]+c01*fxy[i+1][j]+c10*fxy[i][j+1]+c11*fxy[i+1][j+1];
      }
      else
      {
        fprintf(stdout,"Out of range %e %e %e %e %e %e\n",xx,x[0],x[n1-1],yy,y[0],y[n2-1]);
      }
    }
  }
  return ff;
}

double s_interp(double ee, double dd) {

  int i,j;
  double logr, ss,ee1;

   ee1 = log(ee);
   logr = log(dd);

  i = (int) ((ee1-eps_grid[0])*inv_delta_eps );
  if (i < 0)         i=0;
  if (i > Neps - 2) i=Neps-2;

  j = (int) ( (logr-r_grid[0])*inv_delta_r );
  if (j < 0)        j=0;
  if (j > Nr - 2) j=Nr - 2;

  ss =  ((logr-r_grid[j])   * (ee1-eps_grid[i])   * stbl[j+1][i+1]
       +(logr-r_grid[j])   * (eps_grid[i+1]-ee1) * stbl[j+1][i]
       +(r_grid[j+1]-logr) * (ee1-eps_grid[i])   * stbl[j][i+1]
       +(r_grid[j+1]-logr) * (eps_grid[i+1]-ee1) * stbl[j][i])*inv_delta_eps*inv_delta_r;

  return ss = ss;
}

double p_interp(double ee, double dd) {

  int i,j;
  double logr, logp,ee1;

  ee1 = log(ee);
  logr = log(dd);

  i = (int) ((ee1-eps_grid[0])*inv_delta_eps );
  if (i < 0)         i=0;
  if (i > Neps - 2) i=Neps-2;

  j = (int) ( (logr-r_grid[0])*inv_delta_r );
  if (j < 0)        j=0;
  if (j > Nr - 2) j=Nr - 2;

  logp =  ((logr-r_grid[j])   * (ee1-eps_grid[i])   * ptbl[j+1][i+1]
       +(logr-r_grid[j])   * (eps_grid[i+1]-ee1) * ptbl[j+1][i]
       +(r_grid[j+1]-logr) * (ee1-eps_grid[i])   * ptbl[j][i+1]
       +(r_grid[j+1]-logr) * (eps_grid[i+1]-ee1) * ptbl[j][i])*inv_delta_eps*inv_delta_r;

  logp = exp(logp);

  return logp;
}
