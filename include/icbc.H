#ifndef __ICBC_INCLUDED__
#define __ICBC_INCLUDED__

class GridData;
class PhysicsData;
class RunData;
class RTS;

extern void SetInitialConditions(GridData&, const PhysicsData&);
extern void SetBoundaryConditions(const RunData&, GridData&,
				  const PhysicsData&, 
				  const int, const int,RTS*);

#endif
