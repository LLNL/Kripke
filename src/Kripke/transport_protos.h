/* Prototypes */

#include<stdlib.h>
#include<mpi.h>
#include<vector>

struct User_Data;
struct Boltzmann_Solver;

/* directions.c */
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

void EvalSigmaTot(User_Data *user_data, std::vector<double> &vector);

/* sweep.c */
void SweepDD(int d, Grid_Data *grid_data, double *volume,
             std::vector<double> &sigt, double *sors, double *psi,
             double *i_plane_psi, double *j_plane_psi,
             double *k_plane_psi,
             double *psi_lf, double *psi_fr, double *psi_bo);

/* sweep_solver.c */
int SweepSolverSolve(User_Data *user_data, double **rhs, double **ans,
                     double **tempv);
void CreateBufferInfoDD(User_Data *user_data);

