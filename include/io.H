class RunData;
class GridData;
class PhysicsData;
class RTS;

extern int  HasBackupFile(const char*);
extern void OutputSolution(const RunData&, const GridData&, const PhysicsData&);
extern void BackupSolution(const RunData&, const GridData&, const PhysicsData&);
extern void RestoreSolution(RunData&, GridData&, const PhysicsData&);
extern void OutputDerivatives(const RunData&, const GridData&, const PhysicsData&);
extern void RestoreDerivatives(RunData&, GridData&, const PhysicsData&);
extern void IO_Init(const GridData&);
extern void IO_Finalize();

extern void Iout(const RunData&, const GridData&, const PhysicsData&, RTS*);
extern void tau_slice(const RunData&, const GridData&, const PhysicsData&, RTS*);
extern void xy_slice(const RunData&, const GridData&, const PhysicsData&);
extern void yz_slice(const RunData&, const GridData&, const PhysicsData&, RTS*);
extern void xz_slice(const RunData&, const GridData&, const PhysicsData&);

extern void diag_output(const RunData&,const GridData&,const PhysicsData&, RTS*);
extern void eos_output(const RunData&,const GridData&,const PhysicsData&, RTS*);
extern void corona_emission_dem_xyz(const RunData&  Run, const GridData& Grid, 
		               const PhysicsData& Physics);
