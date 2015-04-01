/*--------------------------------------------------------------------------
 * Header file for the User_Data structure
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_USER_DATA_H__
#define KRIPKE_USER_DATA_H__

#include <Kripke/Grid.h>
#include <Kripke/Input_Variables.h>
#include <Kripke/Timing.h>

#include <vector>
#include <mpi.h>


class Kernel;
class Comm;

struct User_Data {
  User_Data(Input_Variables *input_vars);
  ~User_Data();

  Timing timing;

  int niter;

  double source_value;

  std::vector<double> sigma_tot;            // Cross section data

  Grid_Data *grid_data;                     // Spatial grids and variables

  size_t global_num_zones;                  // Total zones across all grids
  int num_group_sets;                       // Number of group-sets
  int num_groups_per_set;                   // How many groups in each set
  int num_direction_sets;                   // Number of direction-sets
  int num_directions_per_set;               // Number of directions per dir set
  int num_zone_sets;

  std::vector<Directions> directions;       // Direction data

  Kernel *kernel;
  Comm *comm;

  void randomizeData();
  void copy(User_Data const &b);
  bool compare(User_Data const &b, double tol, bool verbose);
};

#endif
