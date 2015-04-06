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
*/
Grid_Data::Grid_Data(Input_Variables *input_vars)
{
  // Create object to describe processor and subdomain layout in space
  // and their adjacencies
  Layout *layout = createLayout(input_vars);

  // create the kernel object based on nesting
  kernel = createKernel(input_vars->nesting, 3);

  // Create quadrature set (for all directions)
  int total_num_directions = input_vars->num_dirsets_per_octant * input_vars->num_dirs_per_dirset;
  InitDirections(this, total_num_directions);

  num_direction_sets = 8*input_vars->num_dirsets_per_octant;
  num_directions_per_set = input_vars->num_dirs_per_dirset;
  num_group_sets = input_vars->num_groupsets;
  num_groups_per_set = input_vars->num_groups_per_groupset;
  num_zone_sets = 1;
  for(int dim = 0;dim < 3;++ dim){
    num_zone_sets *= input_vars->num_zonesets_dim[dim];
  }

  legendre_order = input_vars->legendre_order;
  total_num_moments = (legendre_order+1)*(legendre_order+1);

  int num_subdomains = num_direction_sets*num_group_sets*num_zone_sets;

  Nesting_Order nest = input_vars->nesting;

  /* Set ncalls */
  niter = input_vars->niter;

  // setup cross-sections
  int total_num_groups = num_group_sets*num_groups_per_set;
  sigma_tot.resize(total_num_groups, 0.0);

  // Setup scattering transfer matrix for 2 materials
  sigs.resize(2);
  for(int mat = 0;mat < 2;++ mat){
    // allocate transfer matrix
    sigs[mat] = new SubTVec(kernel->nestingSigs(), total_num_groups, legendre_order+1, total_num_groups);

    // Set to identity for all moments
    sigs[mat]->clear(0.0);
    for(int l = 0;l < legendre_order+1;++ l){
      for(int g = 0;g < total_num_groups;++ g){
        (*sigs[mat])(g, l, g) = 1.0;
      }
    }
  }

  // just allocate pointer vectors, we will allocate them below
  ell.resize(num_direction_sets, NULL);
  ell_plus.resize(num_direction_sets, NULL);
  phi.resize(num_zone_sets, NULL);
  phi_out.resize(num_zone_sets, NULL);

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

        // Create ell and ell_plus, if this is the first of this ds
        if(ell[ds] == NULL){
          ell[ds] = new SubTVec(kernel->nestingEll(), total_num_moments, sdom.num_directions, 1);
          ell_plus[ds] = new SubTVec(kernel->nestingEllPlus(), total_num_moments, sdom.num_directions, 1);
        }

        // Create phi and phi_out, if this is the first of this zs
        if(phi[zs] == NULL){
          phi[zs] = new SubTVec(nest, total_num_groups, total_num_moments, sdom.num_zones);
          phi_out[zs] = new SubTVec(nest, total_num_groups, total_num_moments, sdom.num_zones);
        }

        // Set the variables for this subdomain
        sdom.setVars(ell[ds], ell_plus[ds], phi[zs], phi_out[zs]);
      }
    }
  }

  delete layout;

  // Now compute number of elements allocated globally
  long long vec_size[4] = {0,0,0,0};
  for(int sdom_id = 0;sdom_id < subdomains.size();++sdom_id){
    Subdomain &sdom = subdomains[sdom_id];
    vec_size[0] += sdom.psi->elements;
    vec_size[1] += sdom.psi->elements;
  }
  for(int zs = 0;zs < num_zone_sets;++ zs){
    vec_size[2] += phi[zs]->elements;
    vec_size[3] += phi_out[zs]->elements;
  }

  long long vec_global[4];
  MPI_Reduce(vec_size, vec_global, 4, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if(mpi_rank == 0){
    printf("Unknown counts: psi=%ld, rhs=%ld, phi=%ld, phi_out=%ld\n",
      (long)vec_global[0], (long)vec_global[1], (long)vec_global[2], (long)vec_global[3]);
  }
}

Grid_Data::~Grid_Data(){
  delete kernel;
  for(int zs = 0;zs < num_zone_sets;++ zs){
    delete phi[zs];
    delete phi_out[zs];
  }
  for(int ds = 0;ds < num_direction_sets;++ ds){
    delete ell[ds];
    delete ell_plus[ds];
  }
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

  for(int zs = 0;zs < num_zone_sets;++ zs){
    phi[zs]->randomizeData();
    phi_out[zs]->randomizeData();
  }

  for(int ds = 0;ds < num_direction_sets;++ ds){
    ell[ds]->randomizeData();
    ell_plus[ds]->randomizeData();
  }
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

  for(int zs = 0;zs < num_zone_sets;++ zs){
    phi[zs]->copy(*b.phi[zs]);
    phi_out[zs]->copy(*b.phi_out[zs]);
  }

  for(int ds = 0;ds < ell.size();++ ds){
    ell[ds]->copy(*b.ell[ds]);
    ell_plus[ds]->copy(*b.ell_plus[ds]);
  }
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

  for(int zs = 0;zs < num_zone_sets;++ zs){
    is_diff |= phi[zs]->compare("phi", *b.phi[zs], tol, verbose);
    is_diff |= phi_out[zs]->compare("phi_out", *b.phi_out[zs], tol, verbose);
  }

  for(int ds = 0;ds < ell.size();++ ds){
    is_diff |= ell[ds]->compare("ell", *b.ell[ds], tol, verbose);
    is_diff |= ell_plus[ds]->compare("ell_plus", *b.ell_plus[ds], tol, verbose);
  }

  return is_diff;
}


