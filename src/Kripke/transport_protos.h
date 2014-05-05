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

/* sweep_solver.c */
int SweepSolverSolve(User_Data *user_data);
void CreateBufferInfoDD(User_Data *user_data);

