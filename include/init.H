#ifndef __INIT_INCLUDED__
#define __INIT_INCLUDED__
class RunData;
class GridData;
class PhysicsData;
class RTS;

enum {newrun=0, restart=1};

extern int Initialize(RunData&,GridData&,PhysicsData&,RTS*&);

#endif
