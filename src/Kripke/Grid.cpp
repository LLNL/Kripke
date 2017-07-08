/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

#include <Kripke/Grid.h>

#include <Kripke/Comm.h>
#include <Kripke/InputVariables.h>
#include <Kripke/Layout.h>
#include <Kripke/SubTVec.h>
#include <cmath>
#include <sstream>

#ifdef KRIPKE_USE_SILO
#include <sys/stat.h>
#include <silo.h>
#include <string.h>
#endif

/**
 * Grid_Data constructor
*/
Grid_Data::Grid_Data(InputVariables *input_vars)
{
  // Create object to describe processor and subdomain layout in space
  // and their adjacencies
  Layout *layout = createLayout(input_vars);

  // Create quadrature set (for all directions)
  int total_num_directions = input_vars->num_directions;
  InitDirections(this, input_vars);

  num_direction_sets = input_vars->num_dirsets;
  num_directions_per_set = total_num_directions/num_direction_sets;
  num_group_sets = input_vars->num_groupsets;
  num_groups_per_set = input_vars->num_groups/ num_group_sets;
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

  // setup mapping of moments to legendre coefficients
  moment_to_coeff.resize(total_num_moments);
  int nm = 0;
  for(int n = 0;n < legendre_order+1;++ n){
    for(int m = -n;m <= n; ++ m){
      moment_to_coeff[nm] = n;
      ++ nm;
    }
  }

  // setup cross-sections
  int total_num_groups = num_group_sets*num_groups_per_set;
  sigma_tot.resize(total_num_groups, 0.0);

  // Setup scattering transfer matrix for 3 materials  

  sigs = new SubTVec(NEST_DGZ, total_num_groups*total_num_groups, legendre_order+1, 3);

  // Set to isotropic scattering given user inputs
  sigs->clear(0.0);
  for(int mat = 0;mat < 3;++ mat){
    for(int g = 0;g < total_num_groups;++ g){
      int idx_g_gp = g*total_num_groups + g;
      (*sigs)(idx_g_gp, 0, mat) = input_vars->sigs[mat];
    }
  }

  // just allocate pointer vectors, we will allocate them below
  ell.resize(num_direction_sets, NULL);
  ell_plus.resize(num_direction_sets, NULL);

  // Initialize Subdomains
  zs_to_sdomid.resize(num_zone_sets);
  subdomains.resize(num_subdomains);
  for(int gs = 0;gs < num_group_sets;++ gs){
    for(int ds = 0;ds < num_direction_sets;++ ds){
      for(int zs = 0;zs < num_zone_sets;++ zs){
        // Compupte subdomain id
        int sdom_id = layout->setIdToSubdomainId(gs, ds, zs);

        // Setup the subdomain
        Subdomain &sdom = subdomains[sdom_id];
        sdom.setup(sdom_id, input_vars, gs, ds, zs, directions, layout);

        // Create ell and ell_plus, if this is the first of this ds
        bool compute_ell = false;
        if(ell[ds] == NULL){
          ell[ds] = new SubTVec(NEST_ZGD, total_num_moments, sdom.num_directions, 1);
          ell_plus[ds] = new SubTVec(NEST_ZDG, total_num_moments, sdom.num_directions, 1);

          compute_ell = true;
        }

        // setup zs to sdom mapping
        if(gs == 0 && ds == 0){
          zs_to_sdomid[zs] = sdom_id;
        }

        // Set the variables for this subdomain
        sdom.setVars(ell[ds], ell_plus[ds], nullptr, nullptr);

        if(compute_ell){
          // Compute the L and L+ matrices
          sdom.computeLLPlus(legendre_order);
        }
      }
    }
  }
  delete layout;



  // Now compute number of elements allocated globally,
  // and get each materials volume
  double vec_volume[3] = {0.0, 0.0, 0.0};
  for(int zs = 0;zs < num_zone_sets;++ zs){
    int sdom_id = zs_to_sdomid[zs];
    for(int mat = 0;mat < 3;++ mat){
      vec_volume[mat] += subdomains[sdom_id].reg_volume[mat];
    }
  }

  Kripke::Comm comm;

  double global_volume[3];
  for(size_t i = 0;i < 3;++ i){
    global_volume[i] = comm.allReduceSumDouble((long)vec_volume[i]);
  }

  int mpi_rank = comm.rank();
  if(mpi_rank == 0){
    printf("\n");
    printf("Region volumes:\n");
    printf("  Reg1 (source)          %e\n", global_volume[0]);
    printf("  Reg2                   %e\n", global_volume[1]);
    printf("  Reg3                   %e\n", global_volume[2]);
  }
}

Grid_Data::~Grid_Data(){
  for(int ds = 0;ds < num_direction_sets;++ ds){
    delete ell[ds];
    delete ell_plus[ds];
  }
  delete sigs;
}

/**
 * Randomizes all variables and matrices for testing suite.
 */
void Grid_Data::randomizeData(void){
  for(size_t i = 0;i < sigma_tot.size();++i){
    sigma_tot[i] = drand48();
  }

  for(size_t i = 0;i < directions.size();++i){
    directions[i].xcos = drand48();
    directions[i].ycos = drand48();
    directions[i].zcos = drand48();
  }


  for(size_t s = 0;s < subdomains.size();++ s){
    subdomains[s].randomizeData();
  }

  for(int ds = 0;ds < num_direction_sets;++ ds){
    ell[ds]->randomizeData();
    ell_plus[ds]->randomizeData();
  }

  sigs->randomizeData();
}



/**
 * Copies all variables and matrices for testing suite.
 * Correctly copies data from one nesting to another.
 */
void Grid_Data::copy(Grid_Data const &b){
  sigma_tot = b.sigma_tot;
  directions = b.directions;

  subdomains.resize(b.subdomains.size());
  for(size_t s = 0;s < subdomains.size();++ s){
    subdomains[s].copy(b.subdomains[s]);
  }

  for(size_t ds = 0;ds < ell.size();++ ds){
    ell[ds]->copy(*b.ell[ds]);
    ell_plus[ds]->copy(*b.ell_plus[ds]);
  }

  sigs->copy(*b.sigs);
}

/**
 * Compares all variables and matrices for testing suite.
 * Correctly compares data from one nesting to another.
 */
bool Grid_Data::compare(Grid_Data const &b, double tol, bool verbose){
  bool is_diff = false;

  for(size_t i = 0;i < directions.size();++i){
    std::stringstream dirname;
    dirname << "directions[" << i << "]";

    is_diff |= compareScalar(dirname.str()+".xcos",
        directions[i].xcos, b.directions[i].xcos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".ycos",
        directions[i].ycos, b.directions[i].ycos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".zcos",
        directions[i].zcos, b.directions[i].zcos, tol, verbose);
  }

  for(size_t s = 0;s < subdomains.size();++ s){
    is_diff |= subdomains[s].compare(
        b.subdomains[s], tol, verbose);

  }
  is_diff |= compareVector("sigma_tot", sigma_tot, b.sigma_tot, tol, verbose);

  for(size_t ds = 0;ds < ell.size();++ ds){
    is_diff |= ell[ds]->compare("ell", *b.ell[ds], tol, verbose);
    is_diff |= ell_plus[ds]->compare("ell_plus", *b.ell_plus[ds], tol, verbose);
  }

  is_diff |= sigs->compare("sigs", *b.sigs, tol, verbose);

  return is_diff;
}



