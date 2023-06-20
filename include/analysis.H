#ifndef __ANALYSIS_INCLUDED__
#define __ANALYSIS_INCLUDED__

class RunData;
class GridData;
class PhysicsData;
class RTS;

extern void AnalyzeSolution(const RunData&,
			    const GridData&,
			    const PhysicsData&,RTS*);
extern void AnalyzeSolution_VP(const RunData&,
			    const GridData&,
			    const PhysicsData&,RTS*);

#endif
