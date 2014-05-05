/*--------------------------------------------------------------------------
 * Header file for the User_Data structure
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_USER_DATA_H__
#define KRIPKE_USER_DATA_H__

#include <Kripke/grid.h>
#include <Kripke/input_variables.h>
#include <Kripke/timing.h>

#include <vector>
#include <mpi.h>


class Kernel;

struct User_Data {
  User_Data(Input_Variables *input_vars);
  ~User_Data();

  Timing timing;

  int ncalls;

  double source_value;

  int bc_types[6];                          // boundary condition type
  double bc_values[6];                      // Boundary condition value

  std::vector<double> sigma_tot;            // Cross section data

  Grid_Data *grid_data;                     // Spatial grids and variables
    
  size_t global_num_zones;                  // Total zones across all grids
  int num_group_sets;                       // Number of group-sets
  int num_groups_per_set;                   // How many groups in each set
  int num_direction_sets;                   // Number of direction-sets
  int num_directions_per_set;               // Number of directions per dir set
  
  std::vector<Directions> directions;       // Direction data
  std::vector<int> octant_map;              // Direction origination octant
  
  Kernel *kernel;
};

#endif
