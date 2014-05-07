/*--------------------------------------------------------------------------
 * Utility functions for the User_Data structure.
 * Note that user input data is only used in this file.
 *--------------------------------------------------------------------------*/

#include <Kripke/Comm.h>
#include <Kripke/User_Data.h>
#include <Kripke.h>
#include <stdio.h>

/*--------------------------------------------------------------------------
 * AllocUserData : Creates a new User_Data structure and
 *                       allocates memory for its data.
 *--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------
 * comm            : MPI Communicator for use in grid manipulations.
 * input_vars      : Structure with input data.
 *--------------------------------------------------------------------------*/

User_Data::User_Data(Input_Variables *input_vars)
{


  /* Set the processor grid dimensions */
  int R = (input_vars->npx)*(input_vars->npy)*(input_vars->npz);;
  create_R_grid(R);
  MPI_Comm R_group = GetRGroup();

  // Create the spatial grid
  grid_data = new Grid_Data(input_vars, &directions[0]);

  // Create base quadrature set
  InitDirections(this, input_vars->num_dirsets_per_octant * input_vars->num_dirs_per_dirset);

  num_direction_sets = 8*input_vars->num_dirsets_per_octant;
  num_directions_per_set = input_vars->num_dirs_per_dirset;
  num_group_sets = input_vars->num_groupsets;
  num_groups_per_set = input_vars->num_groups_per_groupset;

  // Initialize Group and Direction Set Structures
  grid_data->gd_sets.resize(input_vars->num_groupsets);
  int group0 = 0;
  for(int gs = 0;gs < grid_data->gd_sets.size();++ gs){
    grid_data->gd_sets[gs].resize(8*input_vars->num_dirsets_per_octant);
    int dir0 = 0;
    for(int ds = 0;ds < grid_data->gd_sets[gs].size();++ ds){
      Group_Dir_Set &gdset = grid_data->gd_sets[gs][ds];
      gdset.num_groups = input_vars->num_groups_per_groupset;
      gdset.num_directions = input_vars->num_dirs_per_dirset;

      gdset.group0 = group0;

      gdset.direction0 = dir0;
      gdset.directions = &directions[dir0];

      dir0 += input_vars->num_dirs_per_dirset;
    }

    group0 += input_vars->num_groups_per_groupset;
  }

  /* Set ncalls */
  niter = input_vars->niter;

  // setup cross-sections
  sigma_tot.resize(num_group_sets*num_groups_per_set, 0.0);

  // Compute number of zones
  global_num_zones = (size_t)input_vars->nx 
                   * (size_t)input_vars->ny 
                   * (size_t)input_vars->nz;

  // create the kernel object based on nesting
  kernel = createKernel(NEST_GDZ, 3);

  // Allocate data
  kernel->allocateStorage(this);

  /* Create buffer info for sweeping if using Diamond-Difference */
  CreateBufferInfo(this);
}


User_Data::~User_Data()
{
  /* Free buffers used in sweeping */
  RBufFree();
  delete grid_data;
  delete kernel;
}
