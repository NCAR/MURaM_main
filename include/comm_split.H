class GridData;
class RunData;

extern MPI_Comm cart_comm;
extern MPI_Comm XY_COMM,XZ_COMM,YZ_COMM;
extern MPI_Comm XCOL_COMM,YCOL_COMM,ZCOL_COMM;
extern MPI_Datatype xy_subarray,xz_subarray,yz_subarray;
extern int xy_rank,xz_rank,yz_rank,xcol_rank,ycol_rank,zcol_rank;

extern void comm_split_init(RunData&, const GridData&);
extern void comm_split_finalize();
