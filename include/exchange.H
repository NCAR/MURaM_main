class GridData;
class PhysicsData;

extern void exchange_grid(GridData&, const PhysicsData&, const int);
extern void exchange_B(GridData&);
extern void exchange_single(const GridData&,double*);

extern void exchange_grid_acc(GridData&, const PhysicsData&, const int);
extern void exchange_B_acc(GridData&);
extern void exchange_single_acc(const GridData&,double*);


