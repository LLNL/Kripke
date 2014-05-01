/*--------------------------------------------------------------------------
 * Utility functions for the Boltzmann_Solver structure.
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"

/*--------------------------------------------------------------------------
 * NewBoltzmannSolver: Creates a new Boltzmann_Solver
 *                      structure and allocates memory for its data.
 *--------------------------------------------------------------------------*/

Boltzmann_Solver
*NewBoltzmannSolver(Grid_Data *grid_data, MPI_Comm comm)
/*--------------------------------------------------------------------------
 * grid_data      : Grid information structure
 * comm           : MPI communicator to use for preconditioning (r_group)
 *--------------------------------------------------------------------------*/
{
  Boltzmann_Solver   *boltzmann_solver;
  Sweep_Solver_Data  *sweep_solver_data;

  int *nzones = grid_data->nzones;
  int num_zones = nzones[0]*nzones[1]*nzones[2];

  NEW(boltzmann_solver, 1, Boltzmann_Solver *);
  sweep_solver_data = NewSweepSolverData(grid_data, comm);
  /*---------------------------------------------------------------------
   * Load pointers into boltzmann_solver structure
   *---------------------------------------------------------------------*/
  (boltzmann_solver->comm)               = comm;
  (boltzmann_solver->sweep_solver_data) = sweep_solver_data;

  return(boltzmann_solver);
}


/*----------------------------------------------------------------------
 * BoltzmannSolverSolve
 *----------------------------------------------------------------------*/

int BoltzmannSolverSolve(double **rhs, double **ans, double **tempv,
                         User_Data *user_data)
{
  Grid_Data         *grid_data = user_data->grid_data;
  Boltzmann_Solver *boltzmann_solver = user_data->boltzmann_solver;

  int *nzones              = grid_data->nzones;
  int num_directions = grid_data->num_directions;
  int d;

  /*---------------------------------------------------------------------
   * Call SweepSolverSolve to solve the H Psi = R linear system
   * for all groups and directions.
   *---------------------------------------------------------------------*/
  SweepSolverSolve(user_data, rhs, ans, tempv);

  return(0);
}

/*----------------------------------------------------------------------
 * FreeBoltzmannSolver
 *----------------------------------------------------------------------*/

int FreeBoltzmannSolver (Grid_Data *grid_data,
                         Boltzmann_Solver *boltzmann_solver)
{

  FreeSweepSolverData(boltzmann_solver->sweep_solver_data);
  FREE(boltzmann_solver);

  return(0);
}
