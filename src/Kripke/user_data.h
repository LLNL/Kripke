/*--------------------------------------------------------------------------
 * Header file for the User_Data structure
 *--------------------------------------------------------------------------*/

#ifndef included_user_data
#define included_user_data

#include "sigma_tot.h"
#include "grid.h"
#include "bc_data.h"

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

typedef struct {
  int                  *timing_index;
  int num_timings;

  int nlevels_kba;
  int ncalls;

  double source_value;

  BC_Data              *bc_data;

  Sigma_Tot            *sigma_tot;

  Grid_Data            *grid_data;

  Boltzmann_Solver    *boltzmann_solver;

  Data_Vector          *tmp_sigma_tot;
  Data_Vector          *tmp_source;

  double               *zonal_tmp;
  double              **psi_i_plane;
  double              **psi_j_plane;
  double              **psi_k_plane;

} User_Data;

#endif
