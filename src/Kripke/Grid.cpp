#include <Kripke/Grid.h>

#include <Kripke/Input_Variables.h>
#include <Kripke/Layout.h>
#include <Kripke/SubTVec.h>
#include <cmath>
#include <sstream>


/**
 * Grid_Data constructor.
 * Currently, the spatial grid is calculated so that cells are a uniform
 * length = (xmax - xmin) / nx
 * in each spatial direction.
 *
*/
Grid_Data::Grid_Data(Input_Variables *input_vars)
{
  Layout *layout = createLayout(input_vars);

  /* Set the processor grid dimensions */
  int R = (input_vars->npx)*(input_vars->npy)*(input_vars->npz);

  /* Check size of PQR_group is the same as MPI_COMM_WORLD */
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(R != size){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
      printf("ERROR: Incorrect number of MPI tasks. Need %d MPI tasks.", R);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }


  int npx = input_vars->npx;
  int npy = input_vars->npy;
  int npz = input_vars->npz;

  /* Compute the local coordinates in the processor decomposition */
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int isub_ref = myid % npx;
  int jsub_ref = ((myid - isub_ref) / npx) % npy;
  int ksub_ref = (myid - isub_ref - npx*jsub_ref) / (npx * npy);


  // create the kernel object based on nesting
  kernel = createKernel(input_vars->nesting, 3);

  // Create base quadrature set
  InitDirections(this, input_vars->num_dirsets_per_octant * input_vars->num_dirs_per_dirset);

  num_direction_sets = 8*input_vars->num_dirsets_per_octant;
  num_directions_per_set = input_vars->num_dirs_per_dirset;
  num_group_sets = input_vars->num_groupsets;
  num_groups_per_set = input_vars->num_groups_per_groupset;
  num_zone_sets = 1;



  int nx_g = input_vars->nx;
  int ny_g = input_vars->ny;
  int nz_g = input_vars->nz;


  num_legendre = input_vars->legendre_order;
  total_num_moments = (num_legendre+1)*(num_legendre+1);

  int num_subdomains = num_direction_sets*num_group_sets*num_zone_sets;

  Nesting_Order nest = input_vars->nesting;

  // Initialize Subdomains
  subdomains.resize(num_subdomains);
  for(int gs = 0;gs < num_group_sets;++ gs){
    for(int ds = 0;ds < num_direction_sets;++ ds){
      for(int zs = 0;zs < num_zone_sets;++ zs){
        // Compupte subdomain id
        int sdom_id = layout->setIdToSubdomainId(gs, ds, zs);

        // Setup the subdomain
        Subdomain &sdom = subdomains[sdom_id];
        sdom.setup(sdom_id, input_vars, gs, ds, zs, directions, kernel, layout);
      }
    }
  }

  /* Set ncalls */
  niter = input_vars->niter;

  // setup cross-sections
  sigma_tot.resize(num_group_sets*num_groups_per_set, 0.0);

  // Allocate moments variables
  int total_dirs = num_directions_per_set * num_direction_sets;
  int num_groups = num_groups_per_set * num_group_sets;

  phi = new SubTVec(nest, num_groups, total_num_moments, subdomains[0].num_zones);
  phi_out = new SubTVec(nest, num_groups, total_num_moments, subdomains[0].num_zones);

  ell = new SubTVec(kernel->nestingEll(), total_num_moments, total_dirs, 1);
  ell_plus = new SubTVec(kernel->nestingEllPlus(), total_num_moments, total_dirs, 1);


  delete layout;
}

Grid_Data::~Grid_Data(){
  delete kernel;
  delete phi;
  delete phi_out;
  delete ell;
  delete ell_plus;
}

/**
 * Randomizes all variables and matrices for testing suite.
 */
void Grid_Data::randomizeData(void){
  for(int i = 0;i < sigma_tot.size();++i){
    sigma_tot[i] = drand48();
  }

  for(int i = 0;i < directions.size();++i){
    directions[i].xcos = drand48();
    directions[i].ycos = drand48();
    directions[i].zcos = drand48();
  }


  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].randomizeData();
  }

  phi->randomizeData();
  phi_out->randomizeData();
  ell->randomizeData();
  ell_plus->randomizeData();
}

/**
 * Copies all variables and matrices for testing suite.
 * Correctly copies data from one nesting to another.
 */
void Grid_Data::copy(Grid_Data const &b){
  sigma_tot = b.sigma_tot;
  directions = b.directions;

  subdomains.resize(b.subdomains.size());
  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].copy(b.subdomains[s]);
  }
  phi->copy(*b.phi);
  phi_out->copy(*b.phi_out);
  ell->copy(*b.ell);
  ell_plus->copy(*b.ell_plus);
}

/**
 * Compares all variables and matrices for testing suite.
 * Correctly compares data from one nesting to another.
 */
bool Grid_Data::compare(Grid_Data const &b, double tol, bool verbose){
  bool is_diff = false;




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

  for(int s = 0;s < subdomains.size();++ s){
    is_diff |= subdomains[s].compare(
        b.subdomains[s], tol, verbose);

  }
  is_diff |= compareVector("sigma_tot", sigma_tot, b.sigma_tot, tol, verbose);
  is_diff |= phi->compare("phi", *b.phi, tol, verbose);
  is_diff |= phi_out->compare("phi_out", *b.phi_out, tol, verbose);
  is_diff |= ell->compare("ell", *b.ell, tol, verbose);
  is_diff |= ell_plus->compare("ell_plus", *b.ell_plus, tol, verbose);

  return is_diff;
}


