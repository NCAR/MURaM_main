class RunData;
class GridData;
class PhysicsData;

extern double MHD_Residual(const RunData&,GridData&, const PhysicsData&); 
extern void TVDlimit(const RunData&,GridData&, const PhysicsData&,const double);
