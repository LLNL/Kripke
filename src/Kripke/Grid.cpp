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
  // Create object to describe processor and subdomain layout in space
  // and their adjacencies
  Layout *layout = createLayout(input_vars);

  // create the kernel object based on nesting
  kernel = createKernel(input_vars->nesting, 3);

  // Create quadrature set (for all directions)
  InitDirections(this, input_vars->num_dirsets_per_octant * input_vars->num_dirs_per_dirset);

  num_direction_sets = 8*input_vars->num_dirsets_per_octant;
  num_directions_per_set = input_vars->num_dirs_per_dirset;
  num_group_sets = input_vars->num_groupsets;
  num_groups_per_set = input_vars->num_groups_per_groupset;
  num_zone_sets = 1;
  for(int dim = 0;dim < 3;++ dim){
    num_zone_sets *= input_vars->num_zonesets_dim[dim];
  }

  num_legendre = input_vars->legendre_order;
  total_num_moments = (num_legendre+1)*(num_legendre+1);

  int num_subdomains = num_direction_sets*num_group_sets*num_zone_sets;

  Nesting_Order nest = input_vars->nesting;

  /* Set ncalls */
  niter = input_vars->niter;

  // setup cross-sections
  sigma_tot.resize(num_group_sets*num_groups_per_set, 0.0);

  // Allocate moments variables
  int total_num_groups = num_groups_per_set * num_group_sets;



  ell.resize(num_direction_sets);
  ell_plus.resize(num_direction_sets);
  for(int ds = 0;ds < num_direction_sets;++ ds){
    ell[ds] = new SubTVec(kernel->nestingEll(), total_num_moments, num_directions_per_set, 1);
    ell_plus[ds] = new SubTVec(kernel->nestingEllPlus(), total_num_moments, num_directions_per_set, 1);
  }

  // Initialize Subdomains
  subdomains.resize(num_subdomains);
  for(int gs = 0;gs < num_group_sets;++ gs){
    for(int ds = 0;ds < num_direction_sets;++ ds){
      for(int zs = 0;zs < num_zone_sets;++ zs){
        // Compupte subdomain id
        int sdom_id = layout->setIdToSubdomainId(gs, ds, zs);

        // Setup the subdomain
        Subdomain &sdom = subdomains[sdom_id];
        sdom.setup(sdom_id, input_vars, gs, ds, zs, directions, kernel, layout, ell[ds], ell_plus[ds]);
      }
    }
  }

  phi = new SubTVec(nest, total_num_groups, total_num_moments, subdomains[0].num_zones);
  phi_out = new SubTVec(nest, total_num_groups, total_num_moments, subdomains[0].num_zones);

  delete layout;
}

Grid_Data::~Grid_Data(){
  delete kernel;
  delete phi;
  delete phi_out;
  for(int ds = 0;ds < ell.size();++ ds){
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

  phi->randomizeData();
  phi_out->randomizeData();

  for(int ds = 0;ds < ell.size();++ ds){
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
  phi->copy(*b.phi);
  phi_out->copy(*b.phi_out);

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
  is_diff |= phi->compare("phi", *b.phi, tol, verbose);
  is_diff |= phi_out->compare("phi_out", *b.phi_out, tol, verbose);
  for(int ds = 0;ds < ell.size();++ ds){
    is_diff |= ell[ds]->compare("ell", *b.ell[ds], tol, verbose);
    is_diff |= ell_plus[ds]->compare("ell_plus", *b.ell_plus[ds], tol, verbose);
  }

  return is_diff;
}


