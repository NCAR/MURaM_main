#ifndef __RUN_INCLUDED__
#define __RUN_INCLUDED__

#include "precision.h"

class RunData{
public:
  int rank;
  int zrank;
  int autostep;
  int rttype;
  char anlfile[128];
  char resfile[128];
  char backfile[128];
  char comment[128];
  int iteration;
  int globiter;
  int maxiter;
  int anlfreq;
  int resfreq;
  int backfreq;
  real outcad;
  int slicefreq;
  real dt;
  real time;
  real Tmax;
  real CFL;
  real CFL_tvd;
  real maxWtime;
  int verbose;
  char path_3D[128];
  char path_2D[128];  
  char eos_name[128];
  char kap_name[128];

  int diagnostics;
  int need_diagnostics;

  int DEM;
  int HAVG;
  int RT_HAVG;

  int diag_output[11];
  int eos_output[14];

  RunData();
  ~RunData() {}

  void Init(const int np) {
    rank = np;
    zrank = np;
    autostep = dt==0. ? 1 : 0;
  }

  void Advance() {
    iteration++;
    globiter++;
    time += dt;
  }

  int IsValid() const {
    return globiter < maxiter && time < Tmax;
  }

  double dt_rem() const { 
    if (outcad==0) return 0.0;
    double num = time/outcad;
    int ceil = (int) num;
    if (num == (double) ceil) return 0.0;
    return (ceil+1)*outcad-time;
  }

  int NeedsOutput() const { return !(globiter%resfreq);  }
  int NeedsBackup() const { return !(globiter%backfreq); }
  int NeedsSlice()  const { return (!(globiter%slicefreq) || !(globiter%resfreq)); }

  void Synchronize(real);
  void Show() const;

};

#endif
