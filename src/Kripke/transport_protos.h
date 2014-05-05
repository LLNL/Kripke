/* Prototypes */

#include<stdlib.h>
#include<mpi.h>
#include<vector>

#include<Kripke/directions.h>

struct User_Data;
struct Boltzmann_Solver;
struct Grid_Data;

/* directions.c */
void InitDirections(User_Data *grid_data, int num_directions_per_octant);

/* driver.c */
void SweepDriver(User_Data *user_data);


void EvalSigmaTot(User_Data *user_data, std::vector<double> &vector);

/* sweep.c */
void SweepDD(int d, Grid_Data *grid_data, std::vector<Directions> &directions,
            std::vector<double> const &volume,
             std::vector<double> &sigt, double *sors, double *psi,
             double *i_plane_psi, double *j_plane_psi,
             double *k_plane_psi,
             double *psi_lf, double *psi_fr, double *psi_bo);

/* sweep_solver.c */
int SweepSolverSolve(User_Data *user_data, double **rhs, double **ans,
                     double **tempv);
void CreateBufferInfoDD(User_Data *user_data);

