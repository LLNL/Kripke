/*--------------------------------------------------------------------------
 * Utility functions for the User_Data structure.
 * Note that user input data is only used in this file.
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"

/*--------------------------------------------------------------------------
 * AllocUserData : Creates a new User_Data structure and
 *                       allocates memory for its data.
 *--------------------------------------------------------------------------*/

User_Data *AllocUserData(MPI_Comm comm, Input_Variables *input_vars)
/*--------------------------------------------------------------------------
 * comm            : MPI Communicator for use in grid manipulations.
 * input_vars      : Structure with input data.
 *--------------------------------------------------------------------------*/
{
  FILE *errfp;
  MPI_Comm R_group;
  User_Data *user_data;
  double  ***psi_zonal;
  double    *zonal_tmp;
  double    *phi_in, *phi_out;

  int num_directions;
  int num_zones;
  int i, d, R;
  int       *nzones;
  int nx, ny, nz, npx, npy, npz, nx_l, ny_l, nz_l;

  errfp = stdout;

  /* Allocate user_data structure */
  NEW(user_data, 1, User_Data *);

  /* Initialize timing data */
  user_data->num_timings = 6;
  CNEW(user_data->timing_index, user_data->num_timings, int *);
  user_data->timing_index[SWEEP] = InitializeTiming("Sweeping");

  /* Set the processor grid dimensions */
  R = (input_vars->npx)*(input_vars->npy)*(input_vars->npz);;
  create_R_grid(R);
  R_group = GetRGroup();

  user_data->grid_data =
    GenGrid(input_vars->npx, input_vars->npy, input_vars->npz,
            input_vars->num_directions_per_octant,
            input_vars->xmin, input_vars->xmax,
            input_vars->nx, input_vars->ymin, input_vars->ymax,
            input_vars->ny, input_vars->zmin, input_vars->zmax,
            input_vars->nz, R_group);
  num_directions = user_data->grid_data->num_directions;

  /* Set nlevels_kba */
  user_data->nlevels_kba = input_vars->nlevels_kba;
  if((user_data->nlevels_kba > 0) && (input_vars->npz > 1)){
    if(GetRrank() == 0){
      printf("ERROR: When nlevels_kba > 0, npz must equal 1.\n");
    }
    error_exit(1);
  }

  /* Set ncalls */
  user_data->ncalls = input_vars->ncalls;

  user_data->sigma_tot = NewSigmaTot(input_vars->sigma_total_value);
  user_data->source_value = input_vars->source_value;

  user_data->tmp_source = NewDataVector(user_data->grid_data);
  user_data->tmp_sigma_tot = NewDataVector(user_data->grid_data);

  nzones = user_data->grid_data->nzones;
  nx_l = nzones[0];
  ny_l = nzones[1];
  nz_l = nzones[2];
  num_zones = nzones[0]*nzones[1]*nzones[2];
  /* Zero out tmp_source->data as zero is always used if no volume
     sources are present. tmp_source should not be used for workspace */
  for(i = 0; i < num_zones; i++){
    user_data->tmp_source->data[i] = 0.0;
  }

  double **psi_i_plane, **psi_j_plane, **psi_k_plane;
  int i_plane_zones = ny_l * nz_l;
  int j_plane_zones = nx_l * nz_l;
  int k_plane_zones = nx_l * ny_l;
  NEW( psi_i_plane, num_directions, double ** );
  NEW( psi_j_plane, num_directions, double ** );
  NEW( psi_k_plane, num_directions, double ** );
  for(d=0; d<num_directions; d++){
    NEW( psi_i_plane[d], i_plane_zones, double * );
    NEW( psi_j_plane[d], j_plane_zones, double * );
    NEW( psi_k_plane[d], k_plane_zones, double * );
  }
  user_data->psi_i_plane = psi_i_plane;
  user_data->psi_j_plane = psi_j_plane;
  user_data->psi_k_plane = psi_k_plane;

  NEW( zonal_tmp, num_zones, double * );
  user_data->zonal_tmp = zonal_tmp;

  /*Only accelerate the 0th moment with DSA. Can eventually make
    num_dsa_moments an input parameter. */
  user_data->boltzmann_solver =
    NewBoltzmannSolver(user_data->grid_data, R_group);

  /* Create buffer info for sweeping if using Diamond-Difference */
  if(user_data->nlevels_kba > 0){
    CreateBufferInfoDDKBA(user_data);
  }
  else {
    CreateBufferInfoDD(user_data);
  }

  return(user_data);
}

/*--------------------------------------------------------------------------
 * InitUserData : Initializes user data associated with time, and
 *                      boundary and initial conditions.
 *--------------------------------------------------------------------------*/

int InitUserData(MPI_Comm comm, User_Data *user_data,
                 Input_Variables *input_vars)
/*--------------------------------------------------------------------------
 * user_data : The User_Data for which initializations will be done.
 *--------------------------------------------------------------------------*/
{
  double *values;

  int *types;
  int j;
  int ierr = 0;

  NEW(types, 6, int *);
  NEW(values, 6, double *);
  for(j=0; j<6; j++){
    types[j] = input_vars->bndry_types[j];
    values[j] = input_vars->bndry_values[j];
  }
  if(((types[0] == 1) && (types[1] == 1)) ||
     ((types[2] == 1) && (types[3] == 1)) ||
     ((types[4] == 1) && (types[5] == 1))){
    /* Cannot have opposing reflecting boundaries */
    if(GetRrank() == 0){
      printf("\nERROR: Opposing reflecting boundaries are not allowed.\n\n");
      error_exit(1);
    }
  }

  user_data->bc_data = NewBCData(user_data->grid_data);
  InitBCData(types, values, user_data->grid_data, user_data->bc_data);

  FREE(types);
  FREE(values);

  return(ierr);
}

/*--------------------------------------------------------------------------
 * FreeUserData : Frees a User_Data structure and all memory
 *                      associated with it.
 *--------------------------------------------------------------------------*/

void FreeUserData(User_Data *user_data)
/*--------------------------------------------------------------------------
 * user_data : The User_Data structure to be deallocated.
 *--------------------------------------------------------------------------*/
{
  int *nzones = user_data->grid_data->nzones;
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int num_directions = user_data->grid_data->num_directions;
  int d;

  FreeDataVector(user_data->tmp_source);
  FreeDataVector(user_data->tmp_sigma_tot);
  FreeBoltzmannSolver(user_data->grid_data, user_data->boltzmann_solver);
  FreeBCData(user_data->bc_data);
  for(d=0; d<num_directions; d++){
    FREE(user_data->psi_i_plane[d]);
    FREE(user_data->psi_j_plane[d]);
    FREE(user_data->psi_k_plane[d]);
  }
  FREE(user_data->psi_i_plane);
  FREE(user_data->psi_j_plane);
  FREE(user_data->psi_k_plane);
  FREE(user_data->zonal_tmp);
  FreeSigmaTot(user_data->sigma_tot);
  /* Free buffers used in sweeping */
  RBufFree();
  FreeGrid(user_data->grid_data);
  FREE(user_data);
}
