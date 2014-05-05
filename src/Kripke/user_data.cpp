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
  int dirmin, dirmax;
  int n1, n2;
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
    if( ((ndir+1)%3) > ((ndir+2)%3) ){
      dirmin = (ndir+2)%3;
      dirmax = (ndir+1)%3;
    }
    else {
      dirmin = (ndir+1)%3;
      dirmax = (ndir+2)%3;
    }
    n1 = nzones[dirmin];
    n2 = nzones[dirmax];

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

User_Data::User_Data(MPI_Comm comm, Input_Variables *input_vars)
{
  /* Set the processor grid dimensions */
  int R = (input_vars->npx)*(input_vars->npy)*(input_vars->npz);;
  create_R_grid(R);
  MPI_Comm R_group = GetRGroup();
    
  // Create the spatial grid  
  num_directions = 8*input_vars->num_directions_per_octant;
  num_groups = 1;
  grid_data = new Grid_Data(input_vars, num_directions, num_groups,  R_group);

  // Create base quadrature set
  InitDirections(this, input_vars->num_directions_per_octant);

  /* Set ncalls */
  ncalls = input_vars->ncalls;

  // setup cross-sections
  sigma_tot = input_vars->sigma_total_value;

  // Compute number of zones
  global_num_zones = (size_t)input_vars->nx 
                   * (size_t)input_vars->ny 
                   * (size_t)input_vars->nz;


  /* Create buffer info for sweeping if using Diamond-Difference */
  CreateBufferInfoDD(this);

  // Initialize Boundary Conditions
  InitBCData(input_vars->bndry_types, input_vars->bndry_values, grid_data, this);
}

/*--------------------------------------------------------------------------
 * FreeUserData : Frees a User_Data structure and all memory
 *                      associated with it.
 *--------------------------------------------------------------------------*/

//void FreeUserData(User_Data *user_data)
User_Data::~User_Data()
/*--------------------------------------------------------------------------
 * user_data : The User_Data structure to be deallocated.
 *--------------------------------------------------------------------------*/
{
  /* Free buffers used in sweeping */
  RBufFree();
  delete grid_data;
}
