class GridData;
class PhysicsData;
class RunData;

const int maxstage = 4;

extern void TCheck(const RunData&, GridData&, const PhysicsData&);
extern void Source_Integrate_Tcheck(const RunData&, GridData&, 
                                    const PhysicsData&, const int);
extern void SaveCons(GridData& Grid,const PhysicsData& Physics);
