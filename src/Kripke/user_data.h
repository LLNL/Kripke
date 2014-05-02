/*--------------------------------------------------------------------------
 * Header file for the User_Data structure
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_USER_DATA_H__
#define KRIPKE_USER_DATA_H__

#include <Kripke/grid.h>
#include <Kripke/input_variables.h>

#include <vector>
#include <mpi.h>


/*--------------------------------------------------------------------------
 * Define the User_Data structure.
 * sigma_tot             : Structure holding info for sigma_tot computations
 * grid_data             : Structure holding all grid information
 * boltzmann_solver      : Structure holding solver data
 * bc_data               : Structure holding all boundary condition data
 * tmp_sigma_tot         : Temporary space for sigma_tot values
 * tmp_source            : Temporary space for source values
 * timing_index          : Array of timing indices used for timing
 * num_timings           : The actual number of timings asked for
 * nlevels_kba           : If > 0, then forces the DD sweep in 3D to use
 *                         KBA for sweeps with nlevels_kba levels.
 * ncalls                : The number of calls to the sweep driver routine.
 *--------------------------------------------------------------------------*/

struct User_Data {
  User_Data(MPI_Comm comm, Input_Variables *input_vars);
  ~User_Data();

  int                  *timing_index;
  int num_timings;

  int ncalls;

  double source_value;

  int bc_types[6];
  double bc_values[6];

  double sigma_tot;

  Grid_Data            *grid_data;
};

#endif
