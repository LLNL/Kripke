/*--------------------------------------------------------------------------
 * Utility functions for the User_Data structure.
 * Note that user input data is only used in this file.
 *--------------------------------------------------------------------------*/

#include <Kripke/Comm.h>
#include <Kripke/User_Data.h>
#include <Kripke.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

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
  /* Check size of PQR_group is the same as MPI_COMM_WORLD */
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(R != size){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
      printf("ERROR: Incorrect number of MPI tasks. Need %d MPI tasks.", R);
    }
    error_exit(1);
  }

  // Create the spatial grid
  grid_data = new Grid_Data(input_vars, &directions[0]);

  // Create base quadrature set
  InitDirections(this, input_vars->num_dirsets_per_octant * input_vars->num_dirs_per_dirset);
  //printf("Total directions=%d\n", (int)directions.size());

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
  kernel = createKernel(input_vars->nesting, 3);

  // Allocate data
  kernel->allocateStorage(this);

  /* Create buffer info for sweeping if using Diamond-Difference */
  CreateBufferInfo(this);
}


User_Data::~User_Data()
{
  /* Free buffers used in sweeping */
  delete comm;
  delete grid_data;
  delete kernel;
}

/*
 * Timing timing;

  int niter;

  double source_value;

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
  Comm *comm;
 */
void User_Data::randomizeData(){
  for(int i = 0;i < sigma_tot.size();++i){
    sigma_tot[i] = drand48();
  }

  for(int i = 0;i < directions.size();++i){
    directions[i].xcos = drand48();
    directions[i].ycos = drand48();
    directions[i].zcos = drand48();
  }

  grid_data->randomizeData();
}
/**
 * Copies DATA without changing nesting
 */
void User_Data::copy(User_Data const &b){
  sigma_tot = b.sigma_tot;
  directions = b.directions;
  grid_data->copy(*b.grid_data);
}

bool User_Data::compare(User_Data const &b, double tol, bool verbose){
  bool is_diff = false;

  is_diff |= compareVector("sigma_tot", sigma_tot, b.sigma_tot, tol, verbose);

  for(int i = 0;i < directions.size();++i){
    std::stringstream dirname;
    dirname << "directions[" << i << "]";

    is_diff |= compareScalar(dirname.str()+".xcos",
        directions[i].xcos, b.directions[i].xcos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".ycos",
        directions[i].ycos, b.directions[i].ycos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".zcos",
        directions[i].zcos, b.directions[i].zcos, tol, verbose);
  }

  is_diff |= grid_data->compare(*b.grid_data, tol, verbose);

  return is_diff;
}
