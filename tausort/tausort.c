#include<stdio.h>
#include<stdlib.h> 
#include<string.h>
#include<math.h>
#include "./global_tau.h"

/*            tausort ordnet die ODF-Stufen in Baender  und         */
/*               erstellt Tabellen fuer Bandopazitaeten und         */
/*                  Planck-Funktion B                               */

inline double min(double a, double b) { return a < b ? a : b; }
inline double max(double a, double b) { return a > b ? a : b; }
        

int main(){

  if ((FULLODF)&&(GREY))
    fprintf(stdout," grey can not be combined with full ODF\n");
        
  input();
  printf("%i    ",1);
 
  fprintf(stdout,"Input Finished. Beginning initialisation.\n");
  initialize();
  fprintf(stdout,"Initialisation complete.\n");

  if (FULLODF){
    fprintf(stdout,"Full ODF mode, calculating bands.\n");
    sort_full();
    fprintf(stdout,"band calculation complete.\n");
  } else {
    fprintf(stdout,"Sorting opacities.\n");
    sort();
    fprintf(stdout,"Sorting complete.\n");

    meanop();
    fprintf(stdout,"binned opacities calculated\n");
  }
  fprintf(stdout,"beginning output\n");
  output();
  fprintf(stdout,"output complete\n");
}


/*    *********************************** */
/*    *            input                * */ 
/*    *********************************** */

void input(void){
 FILE *fp, *fp2;
 
 fp = fopen("G2_1D.dat","r");
 for (int i=0; i < KMAX ; i++)
   {
     fscanf(fp, "%le ", &zax[i]);
     fscanf(fp, "%le ", &rho[i]);
     fscanf(fp, "%le ", &p[i]);
     fscanf(fp, "%lf ", &T[i]);
   }
 
 if(Nlam == 328)
   fp = fopen("./p00big2_asplund.bdf","r");
 else
   fp = fopen("./p00lit2_new.asc","r");
 
 for (int lam=0; lam<=Nlam-1; lam++){
   fscanf(fp, "%le", &wbeg[lam]);
   fscanf(fp, "%le", &wend[lam]);
   for (int p=0; p<=Np-1; p++){
     for (int t=0; t<=NT-1; t++){
       for (int bin=0; bin<=Nbin-1; bin++){
         fscanf(fp, "%hd", &ODF[lam][t][p][bin]);
       }
     }
   }
 }

 if(Nlam == 328)
   fp = fopen("./asplund_abs_cont.dat","r");
 else
   fp = fopen("./acont_pT_lit.dat","r");
 
 if(Nlam == 328)
   fp2 = fopen("./asplund_sca_cont.dat","r");
 else
   fp2 = fopen("./acont_pT_lit.dat","r");
 
 for (int lam=0; lam<=Nlam-1; lam++)
   for (int t=0; t<=NT-1; t++)
     for (int p=0; p<=Np-1; p++){
       double temp1,temp2;
       fscanf(fp, "%le", &temp1);
       kcont_pT[lam][t][p]=(float) temp1;
       fscanf(fp2, "%le", &temp2);
       scont_pT[lam][t][p]=(float) temp2;
       acont_pT[lam][t][p] = kcont_pT[lam][t][p]+scont_pT[lam][t][p];
 }
}

/*  **************************************** */
/*  *                 initialize           * */
/*  **************************************** */

/* Define the frequency intervals
   determine kappa (nu,z) and tau(nu,z) */

  void initialize(void){
  int i,j,k,l,m,nu,bin;

  if (GREY) for (bin=0;bin<Nbands;bin++) 
  {
    Level_tau_bot[bin] = -99.0-bin;
    Level_tau_top[bin] = 99.0-bin;
    Level_nu_bot[bin] = 0;
    Level_nu_top[bin] = 999;
  }

  for (nu=0; nu<=Nnu-1; nu++)
  {

    dnu[nu]= Cvac * (1./wbeg[nu] - 1./wend[nu]) * 1.e07  ;
    numid[nu]= .5e0 * Cvac * (1./wbeg[nu] + 1./wend[nu])
      * 1.e07; 
    nuout[nu] = (float) numid[nu];

    /* in Hz */
    for (j=0; j<=NT-1; j++)
    {
      B[nu][j] = 2.e0*H/(pow(Cvac,2)) * pow(numid[nu],3) * 
        1.e0/(exp(H*numid[nu]/(k_B * exp(TENLOG*tab_T[j])) ) - 1.e0) ; 

      dBdT[nu][j] = 2.e0*pow(H,2)/(k_B*pow(Cvac,2))
        *pow(numid[nu],4)/pow(exp(TENLOG*tab_T[j]),2)
        *1.e0/( exp(H*numid[nu]/(k_B * exp(TENLOG*tab_T[j])) )
        +exp(-1.e0*H*numid[nu]/(k_B * exp(TENLOG*tab_T[j])) )
        -2.e0 );
     }  
   }

  for (k=0; k<=KMAX-1; k++)
  {
    double lgT= log10(T[k]);
    double lgp= log10(p[k]);
  
    for (l=0; l<=NT-2; l++)
    {
      if ((lgT >= tab_T[l]) && (lgT < tab_T[l+1]))
      {
        for (m=1; m<=Np-2; m++)
        {
          if ((lgp >= tab_p[m]) && (lgp < tab_p[m+1]))
          {
            for (nu=0; nu<=Nnu-1; nu++)
            {
              kcont[k][nu]= 1.0/( tab_T[l+1]-tab_T[l])*1.0/(tab_p[m+1]-tab_p[m])
              *( (tab_T[l+1]-lgT)* (tab_p[m+1]-lgp) * log10(kcont_pT[nu][l][m])  
              + (lgT- tab_T[l])*( tab_p[m+1]-lgp) * log10(kcont_pT[nu][l+1][m])
              + (tab_T[l+1]-lgT)*(lgp-tab_p[m]) * log10(kcont_pT[nu][l][m+1])
              +  (lgT- tab_T[l])*(lgp-tab_p[m]) * log10(kcont_pT[nu][l+1][m+1])); 
              kcont[k][nu] = exp(TENLOG*kcont[k][nu]);
              
              scont[k][nu]= 1.0/( tab_T[l+1]-tab_T[l])*1.0/(tab_p[m+1]-tab_p[m])
              *( (tab_T[l+1]-lgT)* (tab_p[m+1]-lgp) * log10(scont_pT[nu][l][m])  
              + (lgT- tab_T[l])*( tab_p[m+1]-lgp) * log10(scont_pT[nu][l+1][m])
              + (tab_T[l+1]-lgT)*(lgp-tab_p[m]) * log10(scont_pT[nu][l][m+1])
              +  (lgT- tab_T[l])*(lgp-tab_p[m]) * log10(scont_pT[nu][l+1][m+1])); 
              scont[k][nu] = exp(TENLOG*scont[k][nu]);

              acont[k][nu] = scont[k][nu]+kcont[k][nu];
             
              for (bin=0; bin<=Nbin-1; bin++){
                line_opacity[nu][bin][k] = 1./( tab_T[l+1]-tab_T[l]) * 1./ (tab_p[m+1]- tab_p[m])
                *( (tab_T[l+1]-lgT)* (tab_p[m+1]-lgp) * ODF[nu][l][m][bin]  
                + (lgT- tab_T[l])*( tab_p[m+1]-lgp) * ODF[nu][l+1][m][bin]
                + (tab_T[l+1]-lgT)*(lgp-tab_p[m]) * ODF[nu][l][m+1][bin]
                +  (lgT- tab_T[l])*(lgp-tab_p[m]) * ODF[nu][l+1][m+1][bin] ); 

                line_opacity[nu][bin][k] = exp(TENLOG*1.e-3* line_opacity[nu][bin][k]);
                kap[nu][bin][k] = line_opacity[nu][bin][k] + acont[k][nu];
              }
            }
          }
        }
      }
    }
  }

  dzi[0][1]=0.e0;
  dzi[KMAX-1][0]=0.e0;

  for (i=1 ; i <= KMAX-1 ; i++)
    dzi[i][1] = zax[i-1]-zax[i]; 
  for (i=0 ; i <= KMAX-2 ; i++)
    dzi[i][0] = zax[i]-zax[i+1]; 

  // Interpolation of kappa to calculate Tau
  for (nu=0; nu<=Nnu-1; nu++){
    for (bin=0; bin<=Nbin-1; bin++){
      tau_nu[nu][bin][0]=1.0e-10 ;
      for (i=1; i<=KMAX-1; i++){
        tau_nu[nu][bin][i] =  tau_nu[nu][bin][i-1] + dzi[i][1] *
        ((1.e0/3.e0)*kap[nu][bin][i]*rho[i]
        +(1.e0/6.e0)*kap[nu][bin][i-1]*rho[i]
        +(1.e0/6.e0)*kap[nu][bin][i]*rho[i-1]
        +(1.e0/3.e0)*kap[nu][bin][i-1]*rho[i-1]);          
      } 
    }
  }

  if(Nlam == 328)
    i_5000 = 169;
  else
    i_5000 = 432;

  tau5000[0] = 1.0e-20;
  Kind = KMAX-1;
  for (i=1; i<=KMAX-1; i++){
     tau5000[i] = tau5000[i-1]+dzi[i][1] * 
       ((1.e0/3.e0)*acont[i][i_5000]*rho[i]
        +(1.e0/6.e0)*acont[i-1][i_5000]*rho[i]
        +(1.e0/6.e0)*acont[i][i_5000]*rho[i-1]
        +(1.e0/3.e0)*acont[i-1][i_5000]*rho[i-1]);

     if (log10(tau5000[i]) >= 5.0)
       Kind = (i < Kind) ? i : Kind;
  }

 for (i=0; i<=KMAX-1; i++)
   printf("%i  %E  %E \n", i,  zax[i]*1.0e-8,log10(tau5000[i]));

}
/*  **************************************** */
/*  *                 sort                 * */
/*  **************************************** */

void sort(void){
  int i,j,k,l,nu,bin;
  struct member *pLL; 
  // initialize pFirst to be 0

  for (i=0; i<=Nbands-1; i++)
    pFirst[i] = 0;

  // Loop over frequency points and ODF sub-bins:
 
  for (nu=0; nu<=Nnu-1; nu++){
    for (bin=0; bin<=Nbin-1; bin++){ 
      Index[nu][bin] = 0 ;
      i=0;
      // Find the point at which tau = 1.0 at this frequency
      while (tau_nu[nu][bin][i] < 1.e0)
        i++; 
      
      // Remember Index, Height and tau5000 of tau_nu = 1.0
      i_at_tau_eq_1[nu][bin] = i;
      z_at_tau_eq_1[nu][bin] = zax[i];
      tref_at_tau_eq_1[nu][bin] = tau5000[i];

      // for each band if tau5000 at tau_nu = 1 is less than the 
      // set bin levels then put it in the bin
      for (j=Nbands-1; j >= 0; j--){
        int tau_range = (log10(tref_at_tau_eq_1[nu][bin])>Level_tau_bot[j])&&
          (log10(tref_at_tau_eq_1[nu][bin])<=Level_tau_top[j]);
        int nu_range = (nu>=Level_nu_bot[j])&&
                    (nu<Level_nu_top[j]);
        if ((tau_range)&&(nu_range)){ 
          Index[nu][bin] = j; 
          break; 
        }
      }

      i = Index[nu][bin];

      pLL = (struct member *) malloc(sizeof(struct member));
      pLL->nuix = nu;
      pLL->binix = bin;
      pLL->pNext = 0;

      // If the lowest wavelength of this band has not been set
      // set both PFirst and PLast to PLL
      if (pFirst[i] == 0){
        pFirst[i] = pLL;
        pLast[i] = pFirst[i];
      }
      // Set the pointer of the last wavelength point to this one
      pLast[i]->pNext = pLL;
      // now set P_last to the new one and repeat
      pLast[i] = pLL;
    }
  }
}

/*  **************************************** */
/*  *              meanop                  * */
/*  **************************************** */

// Calculate Band averaged opacities on p,T grid

void meanop(void){
  struct member *pLL;
  double tau_i; 
  double sum;
  // Loop over each temperature point and frequency (wavelength) band
  for (int t=0; t<=NT-1; t++){
    for (int band=0; band<=Nbands-1; band++){
      B_band[band][t] = 0.e0;
      dBdT_band[band][t] = 0.e0;    
      // calculate the total blackbody source function for each bin
    
      for (pLL=pFirst[band]; pLL; pLL=pLL->pNext){
        int nu = pLL->nuix ;
        int bin= pLL->binix;
        B_band[band][t] = B_band[band][t] + B[nu][t]*dnu[nu]*wt[bin];
        dBdT_band[band][t] = dBdT_band[band][t] + dBdT[nu][t]*dnu[nu]*wt[bin];
      }
    }
    // Calculate blackbody source function at 5000A
    B_5000[t] = log(B[i_5000][t]*dnu[i_5000]);
  }

  // Loop over each temperature, pressure point and frequency (wavelength) band
  for (int t=0; t<=NT-1; t++){
    for (int p=0; p<=Np-1; p++){
      for (int band=0; band<=Nbands-1; band++){
        kap_pl[band][t][p] = 0.e0;
        kap_ro[band][t][p] = 0.e0;

        // Sum Planck and Rosseland intensities in each bin.
        for (pLL=pFirst[band]; pLL; pLL=pLL->pNext){
          int nu = pLL->nuix ;
          int bin= pLL->binix;

          kap_pl[band][t][p] = kap_pl[band][t][p] + B[nu][t]*dnu[nu]*wt[bin]*
          (exp(TENLOG*1.e-3*ODF[nu][t][p][bin]) + acont_pT[nu][t][p]); 
          kap_ro[band][t][p] = kap_ro[band][t][p] + dBdT[nu][t]*dnu[nu]*wt[bin]*
          1.e0/(exp(TENLOG*1.e-3*ODF[nu][t][p][bin]) + acont_pT[nu][t][p]);
        }

        // And normalise
        kap_pl[band][t][p] = kap_pl[band][t][p]/B_band[band][t];
        kap_ro[band][t][p] = dBdT_band[band][t]/kap_ro[band][t][p];

        // tau value used to switch smoothly between Rosseland and Planck
        tau_i = kap_ro[band][t][p]*exp(TENLOG*tab_p[p])/2.74e4;
        // If Ross or Grey only save Rosseland mean
        if(ROSS||GREY)
          kap_mean[band][t][p] = log(kap_ro[band][t][p]);
        else
          kap_mean[band][t][p] = log(kap_pl[band][t][p] * pow(2.e0, -1.e0*tau_i/TAU1_2[band])
          + kap_ro[band][t][p]*(1.e0 - pow(2.e0, -1.e0*tau_i/TAU1_2[band])));

      }
      // Tau5000 bin for reference output
      kap_5000[t][p] = log(acont_pT[i_5000][t][p]);
    }
    for (int band=0; band<=Nbands-1; band++)
      B_band[band][t] = log(B_band[band][t]);
  }
}

/*  **************************************** */
/*  *              sort_full                  * */
/*  **************************************** */

void sort_full(void){
  // Loop over each temperature point and frequency (wavelength) band
  for (int j=0; j<=NT-1; j++)
    for (int nu=0; nu<=Nnu-1; nu++)
      for (int bin=0;bin<=Nbin-1;bin++) {
        B_band[nu*Nbin+bin][j] = log(B[nu][j]*dnu[nu]*wt[bin]);
      }
  
  // Calculate blackbody source function at 5000A
  for (int j=0;j<=NT-1;j++)
    B_5000[j] = log(B[i_5000][j]*dnu[i_5000]);

  // Loop over each temperature, pressure point and frequency (wavelength) band
  for (int j=0; j<=NT-1; j++){    
    for (int k=0; k<=Np-1; k++){
      for(int nu=0;nu<=Nnu-1;nu++){
        for (int bin=0; bin<=Nbin-1; bin++){
          kap_mean[nu*Nbin+bin][j][k] =  log(pow(10,(float) 1.e-3*ODF[nu][j][k][bin])
                + acont_pT[nu][j][k]);
        }
      }
      // Tau5000 bin for reference output
      kap_5000[j][k] = log(acont_pT[i_5000][j][k]);
    }
  }
} 

void output(void){
  int i,j,k,l,nu,bin,start,len;
  FILE *outfp;
  char fname[100];
  struct member *pLL;

  int Nbands_out;
  int header[8];

  int tau5000bin = 0;
  int full_odf = 0;
  int back_heating = 0;
  int scatter_on = 0;
  int pp_axis = 0;

  strcpy(fname,"kappa_");

  if(GREY){
    strcat(fname,"grey");
    Nbands_out = 1;
  } else if (!FULLODF){
    Nbands_out = Nbands;
    char tempstr[10];
    sprintf(tempstr,"%i_%s",Nbands_out,"band");
    strcat(fname,tempstr);
  } else if (FULLODF){
    strcat(fname,"fullodf");
    Nbands_out=Nbin*Nnu;
    full_odf = 1;
  }
  
  strcat(fname,".dat");

  outfp = fopen(fname,"w");

  if (BIN5000)
    tau5000bin = 1;

  header[0] = tau5000bin;  // tau5000 bin? (0) no (1) yes
  header[1] = NT;
  header[2] = Np;
  header[3] = Nbands_out;
  header[4] = pp_axis;
  header[5] = full_odf; // FULL ODF? (0) no (1) yes
  header[6] = scatter_on; // Scattering? (0) no (1) yes
  header[7] = back_heating; //  coronal back heating bin? (0) no (1) yes
  
  fprintf(stdout,"Bins are tau5000 %i, NT %i, NT %i, Nbands_out %i, full ODF %i, scatter %i, back_heating %i \n",tau5000bin,NT,Np,Nbands_out,full_odf,scatter_on,back_heating);

  fwrite(header,sizeof(int),8,outfp);

  fwrite(tab_T,sizeof(double),NT,outfp);
  fwrite(tab_p,sizeof(double),Np,outfp);

  if(BIN5000){
    fwrite(kap_5000,sizeof(float),NT*Np,outfp);
    fwrite(B_5000,sizeof(float),NT,outfp);
  }

  if (GREY){
    fwrite(kap_mean[0],sizeof(float),Nbands_out*NT*Np,outfp);
    fwrite(B_band[0],sizeof(float),Nbands_out*NT,outfp);
  } else {
    fwrite(kap_mean,sizeof(float),Nbands_out*NT*Np,outfp);
    fwrite(B_band,sizeof(float),Nbands_out*NT,outfp);
  }

  if (FULLODF)
    fwrite(nuout,sizeof(float),Nnu,outfp);

  fclose(outfp);
}
