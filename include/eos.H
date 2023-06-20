#include "eos_precision.h"

extern int N_eps,N_lr,N_lp3,N_s3;

extern double del_eps,del_lr,del_lp3,del_s3;
extern double inv_del_eps,inv_del_lr,inv_del_lp3,inv_del_s3;
extern double eps_off, ss_off;
extern double *xeps;
extern double *xlr;
extern eos_real **p_eostab;
extern eos_real **T_eostab;
extern eos_real **s_eostab;
extern eos_real **ne_eostab;
extern eos_real **rhoi_eostab;
extern eos_real **amb_eostab;

extern double *xlp3;
extern double *xs3;
extern eos_real **d3_eostab;
extern eos_real **eps3_eostab;

class GridData;
class RunData;

extern double T_interp(double ee, double dd);
#pragma acc routine seq
extern double p_interp(double ee, double dd);
#pragma acc routine seq
extern double s_interp(double ee, double dd);
extern double ne_interp(double ee, double dd);
extern double rhoi_interp(double ee, double dd);
extern double amb_interp(double ee, double dd);
#pragma acc routine seq
extern double d3_interp(double pp, double ss);
#pragma acc routine seq
extern double eps3_interp(double pp, double ss);

extern void eos_init(GridData&, RunData& , const PhysicsData&);
extern void ConsToPrim(GridData&, const PhysicsData&, const RunData&);
