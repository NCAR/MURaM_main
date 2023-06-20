#ifndef __SOLVER_INCLUDED__
#define __SOLVER_INCLUDED__

class RunData;
class GridData;
class PhysicsData;
class RTS;

extern void ComputeSolution(RunData&,GridData&,const PhysicsData&,RTS*);

#endif
