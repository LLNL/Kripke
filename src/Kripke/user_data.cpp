/*--------------------------------------------------------------------------
 * Utility functions for the User_Data structure.
 * Note that user input data is only used in this file.
 *--------------------------------------------------------------------------*/

#include <Kripke/user_data.h>
#include <Kripke/transport_headers.h>


void InitBCData2(int *types, double *vals, Grid_Data *grid_data,
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
  MPI_Comm R_group;
  double  ***psi_zonal;
  double    *phi_in, *phi_out;

  int num_directions;
  int num_zones;
  int i, R;
  int       *nzones;
  int nx, ny, nz, npx, npy, npz, nx_l, ny_l, nz_l;


  /* Initialize timing data */
  num_timings = 6;
  CNEW(timing_index, num_timings, int *);
  timing_index[SWEEP] = InitializeTiming("Sweeping");

  /* Set the processor grid dimensions */
  R = (input_vars->npx)*(input_vars->npy)*(input_vars->npz);;
  create_R_grid(R);
  R_group = GetRGroup();

  grid_data =
    GenGrid(input_vars->npx, input_vars->npy, input_vars->npz,
            input_vars->num_directions_per_octant,
            input_vars->xmin, input_vars->xmax,
            input_vars->nx, input_vars->ymin, input_vars->ymax,
            input_vars->ny, input_vars->zmin, input_vars->zmax,
            input_vars->nz, R_group);
  num_directions = grid_data->num_directions;

  /* Set ncalls */
  ncalls = input_vars->ncalls;

  sigma_tot = input_vars->sigma_total_value;
  source_value = input_vars->source_value;

  tmp_source = NewDataVector(grid_data);


  nzones = grid_data->nzones;
  nx_l = nzones[0];
  ny_l = nzones[1];
  nz_l = nzones[2];
  num_zones = nzones[0]*nzones[1]*nzones[2];
  /* Zero out tmp_source->data as zero is always used if no volume
     sources are present. tmp_source should not be used for workspace */
  for(i = 0; i < num_zones; i++){
    tmp_source->data[i] = 0.0;
  }

  int i_plane_zones = ny_l * nz_l;
  int j_plane_zones = nx_l * nz_l;
  int k_plane_zones = nx_l * ny_l;
  psi_i_plane.resize(num_directions);
  psi_j_plane.resize(num_directions);
  psi_k_plane.resize(num_directions);
  for(int d=0; d<num_directions; d++){
    psi_i_plane[d].resize(i_plane_zones);
    psi_j_plane[d].resize(j_plane_zones);
    psi_k_plane[d].resize(k_plane_zones);
  }

  tmp_sigma_tot.resize(num_zones);
  zonal_tmp.resize(num_zones);

  /*Only accelerate the 0th moment with DSA. Can eventually make
    num_dsa_moments an input parameter. */
  boltzmann_solver =
    NewBoltzmannSolver(grid_data, R_group);

  /* Create buffer info for sweeping if using Diamond-Difference */
  CreateBufferInfoDD(this);


  double *values = input_vars->bndry_values;
  int *types = input_vars->bndry_types;
  InitBCData2(types, values, grid_data, this);

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
  int *nzones = grid_data->nzones;
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int num_directions = grid_data->num_directions;
  int d;

  FreeDataVector(tmp_source);
  FreeBoltzmannSolver(grid_data, boltzmann_solver);

  /* Free buffers used in sweeping */
  RBufFree();
  FreeGrid(grid_data);
}
