/*--------------------------------------------------------------------------
 * Utility functions for the User_Data structure.
 * Note that user input data is only used in this file.
 *--------------------------------------------------------------------------*/

#include <Kripke/comm.h>
#include <Kripke/user_data.h>
#include <Kripke/transport_protos.h>
#include <stdio.h>


static void InitBCData(int *types, double *vals, Grid_Data *grid_data,
                User_Data *user_data)
/*--------------------------------------------------------------------------
 * types     : Array of integers indicating the boundary condition type
 *             on each face for each energy group.
 *                            Type = 0  => Dirichlet
 *                            Type = 1  => Specular
 * vals      : Array of doubles with boundary data
 * grid_data : Grid information
 * bc_data   : The BC_Data structure of boundary condition
 *             information to be returned
 *--------------------------------------------------------------------------*/
{
  int  *nzones       = grid_data->nzones;
  int i1, i2;
  int ndir;
  int face;
  int index;
  int bc_type;

  if(((types[0] == 1) && (types[1] == 1)) ||
     ((types[2] == 1) && (types[3] == 1)) ||
     ((types[4] == 1) && (types[5] == 1))){
    /* Cannot have opposing reflecting boundaries */
    if(GetRrank() == 0){
      printf("\nERROR: Opposing reflecting boundaries are not allowed.\n\n");
      error_exit(1);
    }
  }

  for(ndir = 0; ndir < 3; ndir++){
    for(face = 0; face < 2; face++){
      if( (grid_data->mynbr[ndir][face]) == -1){
        user_data->bc_types[ndir*2 + face] = types[ndir*2 + face];
        bc_type = types[ndir*2 + face];

        if( (bc_type == 0) || (bc_type == 1) ){
          /* Dirichlet or Neumann condition */
          user_data->bc_values[ndir*2 + face] = vals[ndir*2 + face];
        }
        else {
          /* Illegal type */
          error_exit(1);
        }
      }
    }
  }

}


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
  num_group_sets = input_vars->num_groupsets;


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

      group0 += input_vars->num_groups_per_groupset;
      dir0 += input_vars->num_dirs_per_dirset;
    }
  }

  /* Set ncalls */
  ncalls = input_vars->ncalls;

  // setup cross-sections
  sigma_tot = input_vars->sigma_total_value;

  // Compute number of zones
  global_num_zones = (size_t)input_vars->nx 
                   * (size_t)input_vars->ny 
                   * (size_t)input_vars->nz;

  // Initialize Boundary Conditions
  InitBCData(input_vars->bndry_types, input_vars->bndry_values, grid_data, this);

  // create the kernel object based on nesting
  kernel = createKernel(NEST_GDZ, 3);

  // Allocate data
  kernel->allocateStorage(grid_data);

  /* Create buffer info for sweeping if using Diamond-Difference */
  CreateBufferInfoDD(this);
}


User_Data::~User_Data()
{
  /* Free buffers used in sweeping */
  RBufFree();
  delete grid_data;
  delete kernel;
}
