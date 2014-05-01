/* Prototypes */

#include<stdlib.h>
#include<mpi.h>

/* bc_data.c */
BC_Data *NewBCData(Grid_Data *grid_data);
void InitBCData(int *types, double *vals, Grid_Data *grid_data,
                BC_Data *bc_data);
void FreeBCData(BC_Data *bc_data);
void LoadBCDataCB(double **psi, double **psi_ghosted,
                  Grid_Data *grid_data);

/* data_vector.c */
Data_Vector *NewDataVector(Grid_Data *grid_data);
void FreeDataVector(Data_Vector *data_vector);

/* directions.c */
static int Compare(const void *i_v, const void *j_v);
void InitDirections(Grid_Data *grid_data, int num_directions_per_octant);
void InitOmegas(Directions *directions, int *omega_map,
                int *omega_map_inv, int num_directions);
void FreeDirections(Directions *directions);

/* driver.c */
void SweepDriver(User_Data *user_data);

/* grid.c */
Grid_Data *GenGrid(int npx, int npy, int npz,
                   int num_directions,
                   double xmin, double xmax, int nx_g,
                   double ymin, double ymax, int ny_g,
                   double zmin, double zmax, int nz_g,
                   MPI_Comm comm);
void FreeGrid(Grid_Data *grid_data);

/* input_variables.c */
Input_Variables *ReadInput(FILE *in_file);
void PrintInputVariables(Input_Variables *input_variables, FILE *out_file);
void FreeInputVariables(Input_Variables *input_variables);

/* main.c */
int main(int argc, char *argv[]);

/* sweep_kernel.c */
Boltzmann_Solver
*NewBoltzmannSolver(Grid_Data *grid_data, MPI_Comm comm);
int BoltzmannSolverSolve(double **rhs, double **ans, double **tempv,
                         User_Data *user_data);
int FreeBoltzmannSolver(Grid_Data *grid_data,
                        Boltzmann_Solver *boltzmann_solver);

/* sigma_tot.c */
Sigma_Tot *NewSigmaTot(double param0);
void FreeSigmaTot(Sigma_Tot *sigma_tot);
void EvalSigmaTot(Sigma_Tot *sigma_tot, Data_Vector *vector);

/* sweep.c */
void SweepDD(int d, Grid_Data *grid_data, double *volume,
             double *sigt, double *sors, double *psi,
             double *i_plane_psi, double *j_plane_psi,
             double *k_plane_psi,
             double *psi_lf, double *psi_fr, double *psi_bo);

/* sweep_solver.c */
Sweep_Solver_Data *NewSweepSolverData(Grid_Data *grid_data, MPI_Comm comm);
void FreeSweepSolverData (Sweep_Solver_Data *sweep_solver_data);
int SweepSolverSolve(User_Data *user_data, double **rhs, double **ans,
                     double **tempv);
void CreateBufferInfoDD(User_Data *user_data);

/* sweep_solver_kba.c */
int SweepSolverSolveDDKBA(User_Data *user_data, double **rhs, double **ans,
			  double **tempv);
void CreateBufferInfoDDKBA(User_Data *user_data);

/* unit_roundoff.c */
double UnitRoundoff(void);

/* user_data.c */
User_Data *AllocUserData(MPI_Comm comm,
                         Input_Variables *input_vars);
int InitUserData(MPI_Comm comm, User_Data *user_data,
                 Input_Variables *input_vars);
void FreeUserData(User_Data *user_data);
